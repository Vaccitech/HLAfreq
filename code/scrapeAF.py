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
import logging

def makeURL(country, standard='s'):
    base = "http://www.allelefrequencies.net/hla6006a.asp?"
    locus_type = "hla_locus_type=Classical&"
    country = "hla_country=%s&" %(country)
    resolution = "hla_level=2&"
    standard = "standard=%s&" %standard
    url = base + locus_type + country + resolution + standard
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

def formatAF(df):
    if df.sample_size.dtype == "O":
        df.sample_size = pd.to_numeric(df.sample_size.str.replace(",", ""))
    return df

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
    loci = AFtab.loci.unique()
    # Iterate over loci separately
    for locus in loci:
        # Iterate over each dataset reporting that locus
        datasets = AFtab[AFtab.loci == locus][datasetID].unique()
        for dataset in datasets:
            # Single locus, single dataset
            datasetAF = AFtab[(AFtab[datasetID] == dataset) & (AFtab.loci == locus)]
            # What was the sample size for this data?
            dataset_sample_size = datasetAF.sample_size.unique()
            assert len(dataset_sample_size) == 1, "dataset_sample_size must be 1, not %s" %len(dataset_sample_size)
            dataset_sample_size = dataset_sample_size[0]
            # Get all alleles for this locus (across datasets)
            ualleles = AFtab[AFtab.loci == locus].allele.unique()
            # Which of these alleles are not in this dataset?
            missing_alleles = [allele for allele in ualleles if not allele in datasetAF.allele.values]
            missing_rows = [(al, locus, dataset, 0, 0, dataset_sample_size) for al in missing_alleles]
            missing_rows = pd.DataFrame(missing_rows, columns=['allele','loci',datasetID,'allele_freq','carriers%','sample_size'])
            # Add them in with zero frequency
            AFtab = pd.concat([AFtab, missing_rows], ignore_index=True)
    return AFtab

def combineAF(df, weights='sample_size', format=True, add_unmeasured=True, datasetID='population'):
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
    if format:
        df = formatAF(df)
    if add_unmeasured:
        df = unmeasured_alleles(df, datasetID)
    grouped = df.groupby('allele')
    combined = grouped.apply(
        lambda row: [
        row.name,
        row.loci.unique()[0],
        np.average(row.allele_freq, weights=row[weights]),
        row.sample_size.sum()
        ]
    )
    combined = pd.DataFrame(
        combined.tolist(),
        columns = ['allele', 'loci', 'allele_freq', 'sample_size']
    )
    combined = combined.reset_index(drop=True)
    # Check that all alleles in a locus have the same sample size
    # after merging
    if duplicated_sample_size(combined):
        id_duplicated_allele(grouped)

    return combined

def duplicated_sample_size(df):
    """Returns True if any loci has more than 1 unique sample size"""
    locus_sample_sizes = df.groupby('loci').sample_size.apply(lambda x: len(x.unique()))
    return any(locus_sample_sizes != 1)

def id_duplicated_allele(grouped):
    """ Reports the allele that has mupltiple sample sizes """
    duplicated_population = grouped.population.apply(lambda x: any(x.duplicated()))
    assert all(~duplicated_population), "duplicated population within allele %s" %duplicated_population[duplicated_population]


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

# url = "http://www.allelefrequencies.net/hla6006a.asp?hla_selection=A*01%3A01&hla_region=South+Asia"
# base_url = makeURL("Philippines")

# aftab = getAFdata(base_url)
# aftab

# import numpy as np

# aftab.groupby('allele').apply(lambda x: np.average(x.allele_freq, weights=x.sample_size))
# aftab.groupby('allele').size()
