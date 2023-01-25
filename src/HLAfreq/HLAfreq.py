"""
HLAfreq

Download allele frequency data from
[allelefrequencies.net](www.allelefrequencies.net). Allele
frequencies from different populations can be combined to
estimate HLA frequencies of countries or other regions such as
global HLA frequencies.
"""

from bs4 import BeautifulSoup
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy as sp
import logging
import warnings
import matplotlib.colors as mcolors

def makeURL(
        country="",
        standard='s',
        locus="",
        resolution_pattern="bigger_equal_than",
        resolution=2,
        region="",
        ethnic="",
        study_type="",
        dataset_source="",
        sample_year="",
        sample_year_pattern="",
        sample_size="",
        sample_size_pattern=""
    ):
    """Create URL for search of allele frequency net database. All arguments are documented [here](http://www.allelefrequencies.net/extaccess.asp)
    Args:
        country (str, optional): Country name to retrieve records from. Defaults to "".
        standard (str, optional): Filter study quality standard to this or higher. {'g', 's', 'a'} Gold, silver, all. Defaults to 's'.
        locus (str, optional): The locus to return allele data for. Defaults to "".
        resolution_pattern (str, optional): Resolution comparitor {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Filter created using `resolution` and `resolution_pattern`. Defaults to "bigger_equal_than".
        resolution (int, optional): Number of fields of resolution of allele. Filter created using `resolution` and `resolution_pattern`. Defaults to 2.
        region (str, optional): Filter to geographic region. {Asia, Australia, Eastern Europe, ...}. All regions listed [here](http://www.allelefrequencies.net/pop6003a.asp). Defaults to "".
        ethnic (str, optional): Filter to ethnicity. {"Amerindian", "Black", "Caucasian", ...}. All ethnicities listed [here](http://www.allelefrequencies.net/pop6003a.asp). Defaults to "".
        study_type (str, optional): Type of study. {"Anthropology", "Blood+Donor", "Bone+Marrow+Registry", "Controls+for+Disease+Study", "Disease+Study+Patients", "Other", "Solid+Organd+Unrelated+Donors", "Stem+cell+donors"}. Defaults to "".
        dataset_source (str, optional): Source of data. {"Literature", "Proceedings+of+IHWs", "Unpublished"}. Defaults to "".
        sample_year (int, optional): Sample year to compare to. Filter created using sample_year and sample_year_pattern. Defaults to "".
        sample_year_pattern (str, optional): Pattern to compare sample year to. Filter created using sample_year and sample_year_pattern. {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Defaults to "".
        sample_size (int, optional): Sample size to compare to. Filter created using sample_size and sample_size_pattern. Defaults to "".
        sample_size_pattern (str, optional): Pattern to compare sample size to. Filter created using sample_size and sample_size_pattern. {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Defaults to "".

    Returns:
        str: URL to search allelefrequencies.net
    """
    base = "http://www.allelefrequencies.net/hla6006a.asp?"
    locus_type = "hla_locus_type=Classical&"
    hla_locus = "hla_locus=%s&" %(locus)
    country = "hla_country=%s&" %(country)
    region = "hla_region=%s&" %(region)
    ethnic = "hla_ethnic=%s&" %(ethnic)
    study_type = "hla_study=%s&" %(study_type)
    dataset_source = "hla_dataset_source=%s&" %(dataset_source)
    sample_year = "hla_sample_year=%s&" %(sample_year)
    sample_year_pattern = "hla_sample_year_pattern=%s&" %(sample_year_pattern)
    sample_size = "hla_sample_size=%s&" %(sample_size)
    sample_size_pattern = "hla_sample_size_pattern=%s&" %(sample_size_pattern)
    hla_level_pattern = "hla_level_pattern=%s&" %(resolution_pattern)
    hla_level = "hla_level=%s&" %(resolution)
    standard = "standard=%s&" %standard
    url = base + locus_type + hla_locus + country + hla_level_pattern + hla_level + standard + region + ethnic + study_type + dataset_source + sample_year + sample_year_pattern + sample_size + sample_size_pattern
    return url

def parseAF(bs):
    """Generate a dataframe from a given html page

    Args:
        bs (bs4.BeautifulSoup): BeautifulSoup object from allelefrequencies.net page

    Returns:
        pd.DataFrame: Table of allele, allele frequency, samplesize, and population
    """
    # Get the results table from the div `divGenDetail`
    tab = bs.find('div', {'id': 'divGenDetail'}).find('table', {'class': 'tblNormal'})
    # Get the column headers from the first row of the table
    columns = [
        'line', 'allele', 'flag', 'population', 'carriers%',
        'allele_freq', 'AF_graphic', 'sample_size', 'database',
        'distribution','haplotype_association', 'notes'
        ]
    rows =[]
    for row in tab.find_all('tr'):
        rows.append(
            [td.get_text(strip=True) for td in row.find_all('td')]
            )
    # Make dataframe of table rows
    # skip the first row as it's `th` headers
    df = pd.DataFrame(rows[1:], columns = columns)

    # Get HLA loci
    df['loci'] = df.allele.apply(lambda x: x.split("*")[0])

    # Drop unwanted columns
    df = df[['allele', 'loci', 'population', 'allele_freq', 'carriers%', 'sample_size']]
    return df
   
def Npages(bs):
    """How many pages of results are there?

    Args:
        bs (bs4.BeautifulSoup): BS object of allelefrequencies.net results page

    Returns:
        int: Total number of results pages
    """
    # Get the table with number of pages
    navtab = bs.find('div', {'id': 'divGenNavig'}).find('table', {'class': 'table10'})
    assert navtab, f"navtab does not evaluate to True. Check URL returns results in web browser."
    # Get cell with ' of ' in 
    pagesOfN = [td.get_text(strip=True) for td in navtab.find_all('td') if " of " in td.text]
    # Check single cell returned
    assert len(pagesOfN) == 1, "divGenNavig should contain 1 of not %s" %len(pagesOfN)
    # Get total number of pages
    N = pagesOfN[0].split("of ")[1]
    N = int(N)
    return N

def formatAF(AFtab, ignoreG=True):
    """Format allele frequency table.
    Convert sample_size and allele_freq to numeric data type.
    Removes commas from sample size. Removes "(*)" from allele frequency if 
    `ignoreG` is `True`. `formatAF()` is used internally by combineAF and getAFdata by default.

    Args:
        AFtab (pd.DataFrame): Allele frequency data downloaded from allelefrequency.net using `getAFdata()`.
        ignoreG (bool, optional): Treat G group alleles as normal. See http://hla.alleles.org/alleles/g_groups.html for details. Defaults to True.

    Returns:
        pd.DataFrame: The formatted allele frequency data.
    """
    df = AFtab.copy()
    if df.sample_size.dtype == "O":
        df.sample_size = pd.to_numeric(df.sample_size.str.replace(",", ""))
    if df.allele_freq.dtype == "O":
        if ignoreG:
            df.allele_freq = df.allele_freq.str.replace("(*)", "", regex=False)
        df.allele_freq = pd.to_numeric(df.allele_freq)
    return df

def getAFdata(base_url, format=True, ignoreG=True):
    """Get all allele frequency data from a search base url. Iterates over all
        pages regardless of which page is based.

    Args:
        base_url (str): URL for base search.
        format (bool): Format the downloaded data using `formatAF()`.
        ignoreG (bool): treat allele G groups as normal. See http://hla.alleles.org/alleles/g_groups.html for details. Default = True

    Returns:
        pd.DataFrame: allele frequency data parsed into a pandas dataframe
    """
    # Get BS object from base search
    bs = BeautifulSoup(requests.get(base_url).text, 'html.parser')
    # How many pages of results
    N = Npages(bs)
    print("%s pages of results" %N)
    # iterate over pages, parse and combine data from each
    tabs = []
    for i in range(N):
        # print (" Parsing page %s" %(i+1))
        print (" Parsing page %s" %(i+1), end="\r")
        url = base_url + "page=" + str(i+1)
        bs = BeautifulSoup(requests.get(url).text, 'html.parser')
        tab = parseAF(bs)
        tabs.append(tab)
    print("Download complete")
    tabs = pd.concat(tabs)
    if format:
        try:
            tabs = formatAF(tabs, ignoreG)
        except:
            print("Formatting failed, non-numeric datatypes may remain.")
    return tabs

def incomplete_studies(AFtab, llimit=0.95, ulimit=1.1, datasetID='population'):
    """Report any studies with allele freqs that don't sum to 1

    Args:
        AFtab (pd.DataFrame): Dataframe containing multiple studies
        llimit (float, optional): Lower allele_freq sum limit that counts as complete. Defaults to 0.95.
        ulimit (float, optional): Upper allele_freq sum limit that will not be reported. Defaults to 1.1.
        datasetID (str): Unique identifier column for study
    """
    poplocs = AFtab.groupby([datasetID, 'loci']).allele_freq.sum()
    lmask = poplocs < llimit
    if sum(lmask>0):
        print(poplocs[lmask])
        print(f"{sum(lmask)} studies have total allele frequency < {llimit}")
    umask = poplocs > ulimit
    if sum(umask>0):
        print(poplocs[umask])
        print(f"{sum(umask)} studies have total allele frequency > {ulimit}")
    incomplete = pd.concat([poplocs[lmask], poplocs[umask]])
    return incomplete

def only_complete(AFtab, llimit=0.95, ulimit=1.1, datasetID='population'):
    """Returns only complete studies. Studies are only dropped if their population and loci are in noncomplete together.
    This prevents throwing away data if another loci in the population is incomplete

    Args:
        AFtab (pd.DataFrame): Dataframe containing multiple studies
        llimit (float, optional): Lower allele_freq sum limit that counts as complete. Defaults to 0.95.
        ulimit (float, optional): Upper allele_freq sum limit that will not be reported. Defaults to 1.1.
        datasetID (str): Unique identifier column for study. Defaults to 'population'.

    Returns:
        pd.DataFrame: Allele frequency data of multiple studies, but only complete studies.
    """
    noncomplete = incomplete_studies(AFtab=AFtab, llimit=llimit, ulimit=ulimit, datasetID=datasetID)
    # Returns False if population AND loci are in the noncomplete.index
    # AS A PAIR
    # This is important so that we don't throw away all data on a population
    # just because one loci is incomplete.
    complete_mask = AFtab.apply(lambda x: (x[datasetID], x.loci) not in noncomplete.index, axis=1)
    df = AFtab[complete_mask]
    return df

def check_resolution(AFtab):
    """Check if all alleles in AFtab have the same resolution.
    Will print the number of records with each resolution.

    Args:
        AFtab (pd.DataFrame): Allele frequency data

    Returns:
        bool: True only if all alleles have the same resolution, else False.
    """
    resolution = 1 + AFtab.allele.str.count(":")
    resVC = resolution.value_counts()
    pass_check = len(resVC) == 1
    if not pass_check:
        print(resVC)
        print(f"Multiple resolutions in AFtab. Fix with decrease_resolution()")
    return pass_check

def decrease_resolution(AFtab, newres, datasetID='population'):
    """Decrease allele resolution to a specified value so all alleles have the same resolution.

    Args:
        AFtab (pd.DataFrame): Allele frequency data.
        newres (int): The desired number of fields for resolution.
        datasetID (str, optional): Column to use as stud identifier. Defaults to 'population'.

    Returns:
        pd.DataFrame: Allele frequency data with all alleles of requested resolution.
    """
    df = AFtab.copy()
    resolution = 1 + df.allele.str.count(":")
    assert all(resolution >= newres), f"Some alleles have resolution below {newres} fields"
    new_allele = df.allele.str.split(":").apply(lambda x: ":".join(x[:newres]))
    df.allele = new_allele
    collapsed = collapse_reduced_alleles(df, datasetID=datasetID)
    return collapsed

def collapse_reduced_alleles(AFtab, datasetID='population'):
    df = AFtab.copy()
    # Group by alleles within datasets
    grouped = df.groupby([datasetID,'allele'])
    # Sum allele freq but keep other columns
    collapsed = grouped.apply(
        lambda row: [
            sum(row.allele_freq),
            row.sample_size.unique()[0],
            row.loci.unique()[0],
            len(row.loci.unique()),
            len(row.sample_size.unique())
        ]
    )
    collapsed = pd.DataFrame(
        collapsed.tolist(),
        index=collapsed.index,
        columns = ['allele_freq', 'sample_size', 'loci', '#loci', '#sample_sizes']
    ).reset_index()
    # Within a study each all identical alleles should have the same loci and sample size
    assert all(collapsed['#loci'] == 1), "Multiple loci found for a single allele in a single population"
    assert all(collapsed['#sample_sizes'] == 1), "Multiple sample_sizes found for a single allele in a single population"
    collapsed = collapsed[['allele', 'loci','population','allele_freq','sample_size']]
    alleles_unique_in_study(collapsed)
    return collapsed

def unmeasured_alleles(AFtab, datasetID='population'):
    """When combining AF estimates, unreported alleles can inflate frequencies
        so AF sums to >1. Therefore we add unreported alleles with frequency zero.

    Args:
        AFtab (pd.DataFrame): Formatted allele frequency data
        datasetID (str): Unique identifier column for study

    Returns:
        pd.DataFrame: Allele frequency data with all locus alleles reported 
            for each dataset
    """
    df = AFtab.copy()
    loci = df.loci.unique()
    # Iterate over loci separately
    for locus in loci:
        # Iterate over each dataset reporting that locus
        datasets = df[df.loci == locus][datasetID].unique()
        for dataset in datasets:
            # Single locus, single dataset
            datasetAF = df[(df[datasetID] == dataset) & (df.loci == locus)]
            # What was the sample size for this data?
            dataset_sample_size = datasetAF.sample_size.unique()
            assert len(dataset_sample_size) == 1, "dataset_sample_size must be 1, not %s" %len(dataset_sample_size)
            dataset_sample_size = dataset_sample_size[0]
            # Get all alleles for this locus (across datasets)
            ualleles = df[df.loci == locus].allele.unique()
            # Which of these alleles are not in this dataset?
            missing_alleles = [allele for allele in ualleles if not allele in datasetAF.allele.values]
            missing_rows = [(al, locus, dataset, 0, 0, dataset_sample_size) for al in missing_alleles]
            missing_rows = pd.DataFrame(missing_rows, columns=['allele','loci',datasetID,'allele_freq','carriers%','sample_size'])
            # Add them in with zero frequency
            df = pd.concat([df, missing_rows], ignore_index=True)
    return df

def combineAF(AFtab, weights='2n', alpha = [], datasetID='population', format=True, ignoreG=True, add_unmeasured=True, complete=True, resolution=True, unique=True):
    """Combine allele frequencies from multiple studies. `datasetID` is the unique identifier for studies to combine.
    Allele frequencies combined using a Dirichlet distribution where each study's contribution to the concentration parameter is $2 * sample_size * allele_frequency$.
    Sample size is doubled to get `2n` due to diploidy. If an alternative `weights` is set it is not doubled.
    The total concentration parameter of the Dirichlet distribution is the contributions from all studies plus the prior `alpha`.
    If `alpha` is not set the prior defaults to 1 observation of each allele.

    Args:
        AFtab (pd.DataFrame): Table of Allele frequency data
        weights (str, optional): Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.
        alpha (list, optional): Prior to use for Dirichlet distribution. Defaults to [].
        datasetID (str, optional): Unique identifier column for study. Defaults to 'population'.
        format (bool, optional): Run `formatAF()`. Defaults to True.
        ignoreG (bool, optional): Treat allele G groups as normal, see `formatAF()`. Defaults to True.
        add_unmeasured (bool, optional): Add unmeasured alleles to each study. This is important to ensure combined allele frequencies sum to 1. See `add_unmeasured()`. Defaults to True.
        complete (bool, optional): Check study completeness. Uses default values for `incomplete_studies()`. If you are happy with your study completeness can be switched off with False. Defaults to True.
        resolution (bool, optional): Check that all alleles have the same resolution, see `check_resolution()`. Defaults to True.
        unique (bool, optional): Check that each allele appears no more than once per study. See `alleles_unique_in_study()`. Defaults to True.

    Returns:
        pd.DataFrame: Allele frequencies after combining estimates from all studies.
        *allele_freq* is the combined frequency estimate from the Dirichlet mean where the concentration is `alpha` + `c`.
        *alpha* is the prior used for the Dirichlet distribution.
        *c* is the observations used for the Dirichlet distribution.
        *sample_size* is the total sample size of all combined studies.
        *wav* is the weighted average.
    """
    df = AFtab.copy()
    single_loci(df)
    if unique:
        assert alleles_unique_in_study(df, datasetID=datasetID), "The same allele appears multiple times in a dataset"
    if complete:
        assert incomplete_studies(df, datasetID=datasetID).empty, "AFtab contains studies with AF that doesn't sum to 1. Check incomplete_studies(AFtab)"
    if resolution:
        assert check_resolution(df), "AFtab conains alleles at multiple resolutions, check check_resolution(AFtab)"
    if format:
        df = formatAF(df, ignoreG)
    if add_unmeasured:
        df = unmeasured_alleles(df, datasetID)
    try:
        df['2n'] = df.sample_size * 2
    except:
        print("column '2n' could not be created")
    df['c'] =  df.allele_freq * df[weights]
    grouped = df.groupby('allele', sort=True)
    combined = grouped.apply(
        lambda row: [
        row.name,
        row.loci.unique()[0],
        np.average(row.allele_freq, weights=row[weights]),
        row.c.sum(),
        row.sample_size.sum()
        ]
    )
    combined = pd.DataFrame(
        combined.tolist(),
        columns = ['allele', 'loci',
            'wav',
            'c', 'sample_size']
    )
    combined = combined.reset_index(drop=True)
    # Check that all alleles in a locus have the same sample size
    # after merging
    if duplicated_sample_size(combined):
        id_duplicated_allele(grouped)
    if not alpha:
        alpha = default_prior(len(combined.allele))
    combined['alpha'] = alpha
    # Calculate Dirichlet mean for each allele
    combined['allele_freq'] = sp.stats.dirichlet(combined.alpha + combined.c).mean()

    return combined

def default_prior(k):
    """Calculate a default prior, 1 observation of each class.

    Args:
        k (int): Number of classes in the Dirichlet distribution.

    Returns:
        list: List of k 1s to use as prior.
    """
    alpha = [1] * k
    return alpha

def single_loci(AFtab):
    """Check that allele frequency data is only of one locus

    Args:
        AFtab (pd.DataFrame): Allele frequency data
    """
    assert len(AFtab.loci.unique()) == 1, f"'AFtab' must conatain only 1 loci"

def alleles_unique_in_study(AFtab, datasetID='population'):
    """Are all alleles unique in each study? Checks that no alleles are reported more than once in a single study. Study is defined by `datasetID`.

    Args:
        AFtab (pd.DataFrame): Allele frequency data
        datasetID (str, optional): Unique identifier column to define study. Defaults to 'population'.

    Returns:
        bool: `True` on if no alleles occur more than once in any study, otherwise `False`.
    """
    df = AFtab.copy()
    grouped = df.groupby([datasetID,'allele'])
    # Are allele alleles unique? i.e. do any occur multiple times in grouping?
    unique = grouped.size()[grouped.size()>1].empty
    if not unique:
        print(f"Non unique alleles in study, is datasetID correct? {datasetID}")
        print(grouped.size()[grouped.size()>1])
    return unique

def duplicated_sample_size(AFtab):
    """Returns True if any loci has more than 1 unique sample size"""
    locus_sample_sizes = AFtab.groupby('loci').sample_size.apply(lambda x: len(x.unique()))
    return any(locus_sample_sizes != 1)

def id_duplicated_allele(grouped):
    """ Reports the allele that has mupltiple sample sizes """
    duplicated_population = grouped.population.apply(lambda x: any(x.duplicated()))
    assert all(~duplicated_population), "duplicated population within allele %s" %duplicated_population[duplicated_population].index.tolist()

def population_coverage(p):
    """Calculate the proportion of people with at least 1 copy of this allele
        assuming HWE.

    Args:
        p (float): Allele frequency

    Returns:
        float: Sum of homozygotes and heterozygotes for this allele
    """
    q = 1-p
    homo = p**2
    hetero = 2*p*q
    return homo + hetero

def betaAB(alpha):
    """Given the alpha vector defining a Dirichlet distribution calculate the a b values for all composite beta distributions.

    Args:
        alpha (list): Values defining a Dirichlet distribution. This will be the prior (for a naive distribution) or the prior + caf.c for a posterior distribution.

    Returns:
        list: List of a b values defining beta values, i.e. for each allele it is the number of times it was and wasn't observed.
    """
    ab = [(a,sum(alpha)-a) for a in alpha]
    return ab

# def betaCI(a,b,credible_interval=0.95):
#     """Calculat the central credible interval of a beta distribution

#     Args:
#         a (float): Beta shape parameter `a`, i.e. the number of times the allele was observed.
#         b (float): Beta shape parameter `b`, i.e. the number of times the allele was not observed.
#         credible_interval (float, optional): The size of the credible interval requested. Defaults to 0.95.

#     Returns:
#         tuple: Lower and upper credible interval of beta distribution.
#     """
#     bd = sp.stats.beta(a,b)
#     lower_quantile = (1-credible_interval)/2
#     upper_quantile = 1-lower_quantile
#     lower_interval = bd.ppf(lower_quantile)
#     upper_interval = bd.ppf(upper_quantile)
#     return lower_interval, upper_interval

# def AFci(caf, credible_interval=0.95):
#     """Calculate credible interval for combined allele frequency table.
#     Note that this ignores sampling error so confidence interval is too tight.
#     Use HLAhdi.AFhdi() instead.

#     Args:
#         caf (pd.DataFrame): Table produced by combineAF()
#         credible_interval (float, optional): The desired confidence interval. Defaults to 0.95.

#     Returns:
#         list: Lower and upper credible intervals as a list of tuples
#     """
#     ab = betaAB(
#         caf.alpha + caf.c,
#     )
#     ci = [betaCI(a, b, credible_interval) for a,b in ab]
#     return ci

def plot_prior(concentration, ncol=2, psteps=1000, labels=""):
    """Plot probability density function for prior values.

    Args:
        concentration (list): Vector of the prior Dirichlet concentration values.
        ncol (int, optional): Number of columns. Defaults to 2.
        labels (list, optional): Labels for elements of concentration in the same order. Defaults to "".
    """
    ab = betaAB(concentration)
    pline = np.linspace(0, 1, psteps)
    nrow = math.ceil(len(concentration)/ncol)
    fig, ax = plt.subplots(nrow, ncol, sharex=True)
    fig.suptitle("Probability density")
    # If labels is a list nothing happens,
    # But if it's a series it converts to a list
    labels = list(labels)
    if not labels:
        labels = [""] * len(concentration)
    assert len(concentration) == len(labels), "concentration must be same length as labels"
    for i,alpha in enumerate(concentration):
        a,b = ab[i]
        bd = sp.stats.beta(a,b)
        pdf = [bd.pdf(p) for p in pline]
        ax.flat[i].plot(pline, pdf)
        ax.flat[i].set_title(labels[i])
    for axi in ax[-1,:]:
        axi.set(xlabel='Allele freq')
    for axi in ax[:,0]:
        axi.set(ylabel='PDF')
    plt.show()

def plotAF(
    caf=pd.DataFrame(),
    AFtab=pd.DataFrame(),
    cols = list(mcolors.TABLEAU_COLORS.keys()),
    datasetID="population",
    weights = "2n",
    credible_interval=None,
    ):
    """Plot combined allele frequencies, individual allele frequencies,
    and credible intervals on combined allele frequency estimates.
    Credible interval is only plotted if a value is given for `credible_interval`.
    Credible interval is calculated using HLAfreq_pymc.AFhdi() and PyMc.

    Args:
        caf (pd.DataFrame, optional): Combined allele frequency estimates from HLAfreq.combineAF. Defaults to pd.DataFrame().
        AFtab (pd.DataFrame, optional): Table of allele frequency data. Defaults to pd.DataFrame().
        cols (list, optional): List of colours to use for each individual dataset. Defaults to list(mcolors.TABLEAU_COLORS.keys()).
        datasetID (str, optional): Column used to define separate datasets. Defaults to "population".
        weights (str, optional): Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.
        credible_interval (float, optional): The credible interval to plot. Defaults to None.
    """
    # Plot allele frequency for each dataset
    if not AFtab.empty:
        # Cols must be longer than the list of alleles
        # If not, repeat the list of cols
        repeat_cols = np.ceil(len(AFtab[datasetID])/len(cols))
        repeat_cols = int(repeat_cols)
        cols = cols * repeat_cols
        # Make a dictionary mapping datasetID to colours
        cmap = dict(zip(AFtab[datasetID].unique(), cols))
        plt.scatter(
            x = AFtab.allele_freq,
            y = AFtab.allele,
            c = [cmap[x] for x in AFtab[datasetID]],
            alpha = 0.7,
            zorder=2
            )
    # Plot combined allele frequency
    if not caf.empty:
        plt.scatter(
            x = caf.allele_freq,
            y = caf.allele,
            edgecolors='black',
            facecolors='none',
            zorder=3
            )
    # Plot high density interval
    if credible_interval:
        assert not AFtab.empty, "AFtab is needed to calculate credible interval"
        from HLAfreq import HLAfreq_pymc as HLAhdi
        print("Fitting model with PyMC, make take a few seconds")
        hdi = HLAhdi.AFhdi(AFtab=AFtab, weights=weights, datasetID=datasetID, credible_interval=credible_interval)
        for interval in hdi.iterrows():
            # .iterrows returns a index and data as a tuple for each row
            plt.hlines(
                y = interval[1]['allele'],
                xmin = interval[1]['lo'],
                xmax = interval[1]['hi'],
                color="black"
            )
    plt.xlabel("Allele frequency")
    plt.grid(zorder=0)
    plt.show()
