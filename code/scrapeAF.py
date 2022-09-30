"""
Scrape Allele [allelefrequencies.net](www.allelefrequencies.net)
for HLA frequencies by population or larger region to calculate
the regional HLA frequencies e.g. global.

The database can be searched based on url as described in
[automated access](http://www.allelefrequencies.net/extaccess.asp).
"""

from bs4 import BeautifulSoup
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import dirichlet, beta
import logging

def makeURL(country, standard='s', locus="", resolution_pattern="bigger_equal_than", resolution=2):
    base = "http://www.allelefrequencies.net/hla6006a.asp?"
    locus_type = "hla_locus_type=Classical&"
    hla_locus = "hla_locus=%s&" %(locus)
    country = "hla_country=%s&" %(country)
    hla_level_pattern = "hla_level_pattern=%s&" %(resolution_pattern)
    hla_level = "hla_level=%s&" %(resolution)
    standard = "standard=%s&" %standard
    url = base + locus_type + hla_locus + country + hla_level_pattern + hla_level + standard
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

def getAFdata(base_url):
    """Get all allele frequency data from a search base url. Iterates over all
        pages regardless of which page is based.

    Args:
        base_url (str): URL for base search

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
    tabs = pd.concat(tabs)
    return tabs

def formatAF(AFtab):
    df = AFtab.copy()
    if df.sample_size.dtype == "O":
        df.sample_size = pd.to_numeric(df.sample_size.str.replace(",", ""))
    return df

def incomplete_studies(AFtab, llimit=0.95, ulimit=1.1, datasetID='population'):
    """Report any studies with allele freqs that don't sum to 1

    Args:
        df (pd.DataFrame): Dataframe containing multiple studies
        llimit (float, optional): Lower allele_freq sum limit that counts as complete. Defaults to 0.95.
        ulimit (float, optional): Upper allele_freq sum limit that will not be reported. Defaults to 1.1.
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

def check_resolution(AFtab):
    resolution = 1 + AFtab.allele.str.count(":")
    resVC = resolution.value_counts()
    pass_check = len(resVC) == 1
    if not pass_check:
        print(resVC)
        print(f"Multiple resolutions in AFtab. Fix with decrease_resolution()")
    return pass_check

def decrease_resolution(AFtab, newres, datasetID='population'):
    df = AFtab.copy()
    resolution = 1 + df.allele.str.count(":")
    assert all(resolution >= newres), f"Some alleles have resolution below {newres} fields"
    new_allele = df.allele.str.split(":").apply(lambda x: ":".join(x[:newres]))
    df.allele = new_allele
    collapsed = collapse_reduced_alleles(df, datasetID=datasetID)
    return collapsed

def collapse_reduced_alleles(AFtab, datasetID='population'):
    df = AFtab.copy()
    # Group by alleles withing datasets
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

def combineAF(AFtab, weights='2n', alpha = [], datasetID='population', format=True, add_unmeasured=True, complete=True, resolution=True, unique=True):
    """Combine allele frequencies at multiple levels

    Args:
        df (pd.DataFrame): Table of Allele frequency data
        weights (str): Column to weight averages by. Default 'sample_size'
        format (boolean): run formatAF(), Default = True
        add_unmeasured (boolean): run unmeasured_alleles(), Default = True
        datasetID (str): Unique identifier column for study used by unmeasured_alleles()

    Returns:
        pd.DataFrame: Table of allele frequency data with alleles
        grouped to make a weighted average based on weights.
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
        df = formatAF(df)
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
    combined['allele_freq'] = dirichlet(combined.alpha + combined.c).mean()

    return combined

def default_prior(k):
    alpha = [1] * k
    return alpha

def single_loci(AFtab):
    assert len(AFtab.loci.unique()) == 1, f"'AFtab' must conatain only 1 loci"

def alleles_unique_in_study(AFtab, datasetID='population'):
    df = AFtab.copy()
    grouped = df.groupby([datasetID,'allele'])
    # Are allele alleles unique? i.e. do any occur multiple times in grouping?
    unique = grouped.size()[grouped.size()>1].empty
    if not unique:
        print("Non unique alleles in population")
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

def plotAFprob(alpha, caf=pd.DataFrame(), AFtab=pd.DataFrame(), datasetID="population", log=False, psteps=1000, ncol=2, ci=0.95):
    """Plot the (log) posterior density function of all frequencies
        for all alleles based on supplied alpha for Dirichlet
        distribution. Options for adding empirical values
        and plotting the log pdf.

    Args:
        alpha (list): Dirichlet concentration parameters, can be prior + caf.c
        caf (pd.DataFrame, optional): Combined allele frequence data produced by scrapeAF.combineAF(). Defaults to pd.DataFrame().
        AFtab (pd.DataFrame, optional): The uncombined allele frequency data used by scrapeAF.combinedAF(). You must use the same dataframe as this function doesn't have the error checking that scrapeAF.combineAF() has. Defaults to pd.DataFrame().
        datasetID (str, optional): The column used to define datasets. Defaults to "population".
        log (bool, optional): Plot log pdf instead of pdf? Defaults to False.
        psteps (int, optional): Number of increments in pdf calculation, higher values make smoother plots. Defaults to 1000.
        ncol (int, optional): How many columns to arrange subplots in. Defaults to 2.
        ci (float, optional): Central credible interval to plot. Set as 0 to hide. Defaults to 0.95.
    """
    # Get beta parameters for each k in Dirichlet
    ab = betaAB(alpha)
    pline = np.linspace(0,1,psteps)
    if not AFtab.empty:
        # Format the population allele frequencies
        df = AFtab.copy()
        df = unmeasured_alleles(df, datasetID=datasetID)
        df = df.sort_values('allele')
        assert all(df.groupby('population').allele.apply(list).apply(lambda x: x == caf.allele.tolist())), "Alleles not matching between AFtab and caf"
    fig, axs = plt.subplots(math.ceil(len(alpha)/ncol), ncol)
    for i,x in enumerate(alpha):
        subplotselector = i//ncol, i%ncol
        a,b = ab[i]
        bd = beta(a,b)
        if log:
            pdf = [bd.logpdf(p) for p in pline]
        else:
            pdf = [bd.pdf(p) for p in pline]
        ax = axs[subplotselector]
        if not AFtab.empty:
            # Add the empirical allele frequency for each population
            for af in df.groupby('population').allele_freq.apply(list).apply(lambda x: x[i]):
                # x is the reported allele freq, y is it's pdf with scatter
                if log:
                    ax.scatter(af, bd.logpdf(af)*(0.95+np.random.random()/5))
                else:
                    ax.scatter(af, bd.pdf(af)*(0.95+np.random.random()/5))
        ax.plot(pline, pdf)
        # Annotate with alpha
        ax.text(0.5, max(pdf)/2, f"alpha {round(alpha[i])}")
        if not caf.empty:
            # Plot the combined average
            ax.axvline(caf.allele_freq[i], color="black", ls="--")
        if ci:
            cl,cu = betaCI(a,b, ci)
            ax.axvline(cl, color="black", linestyle="dotted")
            ax.axvline(cu, color="black", linestyle="dotted")
        # fig.tight_layout()
    plt.show()

def betaCI(a,b,credible_interval=0.95):
    """Calculat the central credible interval of a beta distribution

    Args:
        a (float): Beta shape parameter `a`, i.e. the number of times the allele was observed.
        b (float): Beta shape parameter `b`, i.e. the number of times the allele was not observed.
        credible_interval (float, optional): The size of the credible interval requested. Defaults to 0.95.

    Returns:
        _type_: _description_
    """
    bd = beta(a,b)
    lower_quantile = (1-credible_interval)/2
    upper_quantile = 1-lower_quantile
    lower_interval = bd.ppf(lower_quantile)
    upper_interval = bd.ppf(upper_quantile)
    return lower_interval, upper_interval

# url = "http://www.allelefrequencies.net/hla6006a.asp?hla_selection=A*01%3A01&hla_region=South+Asia"
# base_url = makeURL("Philippines")

# aftab = getAFdata(base_url)
# aftab

# import numpy as np

# aftab.groupby('allele').apply(lambda x: np.average(x.allele_freq, weights=x.sample_size))
# aftab.groupby('allele').size()
