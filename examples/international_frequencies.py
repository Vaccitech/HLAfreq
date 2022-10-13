"""
Calculate International average allele frequencies
using population estimates.
"""

import code.HLAfreq as HLAfreq
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
    AFtab = AFtab[AFtab.loci == "A"]
    # Weighted average alleles
    wav = HLAfreq.combineAF(AFtab)
    # Add country to dataset
    wav['country'] = country
    wavs.append(wav)

AFtab[AFtab.allele=="A*66:02"]
AFtab[AFtab.allele=="A*66:03"]


international = pd.concat(wavs, ignore_index=True)
# Add population size to weight averages
international['population_size'] = international.country.apply(lambda x: population_sizes[x])

#International Weighted Average
iwav = HLAfreq.combineAF(international, weights='population_size', datasetID='country')

plt.scatter(
    international.allele,
    international.allele_freq,
    c=[international.country.unique().tolist().index(x) for x in international.country],
    s=international.population_size/10000000
    )
plt.scatter(
    iwav.allele,
    iwav.allele_freq,
    c='black',
    marker='x'
)
plt.xticks(rotation=90)
plt.show()
