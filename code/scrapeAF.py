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
import logging

def makeURL(country):
    url = "http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=Classical&hla_country=%s&hla_level=2&" %(country)
    return url

def parseAF(bs):
    """Generate a dataframe from a given url

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
    tabs = formatAF(tabs)
    validateAF(tabs)
    return tabs

def formatAF(df):
    # Force convert columns to numeric
    df.allele_freq = pd.to_numeric(df.allele_freq)
    df.sample_size = pd.to_numeric(df.sample_size)
    return df

def validateAF(df):
    if not all(df.allele_freq <= 1):
        print("Allele frequency greater 1 reported")


url = "http://www.allelefrequencies.net/hla6006a.asp?hla_selection=A*01%3A01&hla_region=South+Asia"
base_url = makeURL("Philippines")

# aftab = getAFdata(base_url)
aftab

import numpy as np

aftab.groupby('allele').apply(lambda x: np.average(x.allele_freq, weights=x.sample_size))
aftab.groupby('allele').size()
