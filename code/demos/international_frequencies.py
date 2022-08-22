"""
Calculate International average allele frequencies
using population estimates.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

countries = [
    'United+Kingdom',
    'Thailand',
    'Uganda'
]
population_sizes = {
    'United+Kingdom': 67081234,
    'Thailand': 66813717,
    'Uganda': 42885900
}

# Get weighted average for each county
wavs = []
for country in countries:
    # Load raw county data
    AFtab = pd.read_csv("data/example/%s_raw.csv" %country)
    # Weighted average alleles
    wav = scrapeAF.combineAF(AFtab)
    # Add country to dataset
    wav['country'] = country
    wavs.append(wav)

international = pd.concat(wavs, ignore_index=True)
international = scrapeAF.unmeasured_alleles(international, 'country')
# Add population size to weight averages
international['population_size'] = international.country.apply(lambda x: population_sizes[x])

#International Weighted Average
iwav = scrapeAF.combineAF(international, 'population_size')

intA = international[international.loci == 'A']
iwava = iwav[iwav.loci == 'A']

plt.scatter(
    intA.allele,
    intA.allele_freq,
    c=[intA.country.unique().tolist().index(x) for x in intA.country],
    s=intA.population_size/10000000
    )
plt.scatter(
    iwava.allele,
    iwava.allele_freq,
    c='black',
    marker='x'
)
plt.xticks(rotation=90)
plt.show()
