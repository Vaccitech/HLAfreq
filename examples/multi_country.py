"""
Calculate multi-country average

We can calculate the average allele frequency of multiple countries
by combining studies within countries and then between countries.
This second step is very similar except we set `datasetID` as the country
rather than the population which is the default. 
"""

import code.HLAfreq as HLAfreq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

countries = ['Cameroon','Cape+Verde','Ghana','Guinea',
    'Guinea-Bissau', 'Kenya','Sao+Tome+and+Principe','Senegal',
    'South+Africa','Uganda','Zimbabwe']

# Download HLA allele frequencies
for country in countries:
    print(country)
    base_url = HLAfreq.makeURL(
        country, standard='s', locus="A",
        resolution_pattern="bigger_equal_than", resolution=2)
    aftab = HLAfreq.getAFdata(base_url)
    aftab.to_csv("data/example/multi_country/%s_raw.csv" %country, index=False)

# Combine alleles frequencies within country
cafs = []
for country in countries:
    # Load raw county data
    aftab = pd.read_csv("data/example/multi_country/%s_raw.csv" %country)
    # Drop any incomplete studies
    aftab = HLAfreq.only_complete(aftab)
    # Ensure all alleles have the same resolution
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    # Combine studies within country
    caf = HLAfreq.combineAF(aftab)
    # Add country name to dataset, this is used as `datasetID` going forward
    caf['country'] = country
    cafs.append(caf)

# Concatenate all combined single country allele frequency data
# into a single dataframe
cafs = pd.concat(cafs, ignore_index=True)
# Combine allele frequency data of all countries
international = HLAfreq.combineAF(cafs, datasetID='country')

# Plot international averages as bar plot
international.plot.barh('allele', 'allele_freq')
plt.show()

# Plot national averages as grouped bar plot
cafs.pivot(index='allele', columns='country', values='allele_freq').plot.bar()
plt.show()

# Plot international allele frequencies estimates and individual countries
HLAfreq.plotAFprob(international, cafs, datasetID='country', ncol=4)

# Plot specific alleles and zoom in on frequencies
# Select alleles to plot
hifreq = international[international.allele_freq > 0.01].allele
# Must be a list
hifreq = hifreq.tolist()
# Plot only selected alleles
HLAfreq.plotAFprob(
    international,
    cafs,
    datasetID='country',
    ncol=4,
    xmax=cafs.allele_freq.max(),
    alleles=hifreq)

"""
From these plots we can clearly see that one country has much
higher frequencies of A*24:02, A*11:01, and A*34:01.
It is also clear that the international average allele
frequency for A*24:02 is being skewed by this country. We can view the
`cafs` dataset to see which country this is, Guinea.

This method can
also be applied to individual countries to identify studies which
differ from the majority, possibly because it focuses on a specific
ethnic group.
"""

cafs[cafs.allele == "A*24:02"].sort_values('allele_freq')

# How to weight countries differently?

"""
When calculating the average allele frequency the "concentration"
parameter of the Dirichlet distribution is calculated as `weights`
multiplied by the dataset allele frequency (where dataset is usually
a single study but can be a country etc, this is what `datasetID`
refers to). By default `weights` is double the sample size as each
sample provides 2n due to diploidy.

If estimating allele frequency for a region, it may be important to
account for the size of national populations. Below we calculate an
individual weight for each country, individuals from large countries
contibute more to the total sample but the total sample size is
unchanged. A county's populaion, as a proportion of the sum of
countries' populations, multiplied by the number of countries is
used as the weight for each individual in that country.
"""

population_sizes = {'Cameroon':24348251,
    'Cape+Verde':563198,
    'Ghana':30832019,
    'Guinea':12907395,
    'Guinea-Bissau':1646077,
    'Kenya':47564296,
    'Sao+Tome+and+Principe':214610,
    'Senegal':17223497,
    'South+Africa':60604992,
    'Uganda':42885900,
    'Zimbabwe':15178979
}

country_data = pd.DataFrame(
    {'country':population_sizes.keys(),
    'population':population_sizes.values()}
)

# What proportion of the regional population does each country account for
country_data['proportion'] = country_data.population/country_data.population.sum()
# How much will each individual in the country count towards the sample size?
country_data['individual_weight'] = country_data['proportion'] * len(country_data.country)
# Add country data to Combined Allele Frequency data
cafs = pd.merge(cafs, country_data, how="left", on='country')
# Sample size is multiplied by this individual weight and doubled
# this accounts for diploid samples from each individual
cafs['weighted_sample_size'] = cafs.sample_size * 2 * cafs.individual_weight

# Combine allele frequency data of all countries
# Using the custom column 'weighted_sample_size'
# The interpretation of the column used for `weights` is that n
# samples were observed with allele x. If unspecified this defaults
# to double the sample size (to accound for diploid samples from each
# individual)
winternational = HLAfreq.combineAF(cafs, datasetID='country', weights='weighted_sample_size')

HLAfreq.plotAFprob(international, cafs, datasetID='country', ncol=4)

hifreq = winternational[winternational.allele_freq > 0.01].allele
# Must be a list
hifreq = hifreq.tolist()
# Plot only selected alleles
HLAfreq.plotAFprob(
    winternational,
    cafs,
    datasetID='country',
    ncol=4,
    xmax=cafs.allele_freq.max(),
    alleles=hifreq)
