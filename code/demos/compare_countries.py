"""
Compare HLA frequencies of different countries.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

countries = [
    'United+Kingdom',
    'Thailand',
    'Uganda'
]

# Download HLA allele frequencies
for country in countries:
    print(country)
    base_url = scrapeAF.makeURL(country)
    aftab = scrapeAF.getAFdata(base_url)
    aftab.to_csv("data/example/%s_raw.csv" %country, index=False)

wavs = []
for country in countries:
    df = pd.read_csv("data/example/%s_raw.csv" %country)
    df = df[df.loci=="A"]
    wav = scrapeAF.combineAF(df)
    wav['country'] = country
    wavs.append(wav)

wavs = pd.concat(wavs, axis=0).reset_index()

# Keep only alleles with frequency over threshold
threshold = 0.02
common_alleles = wavs[wavs.allele_freq > threshold].allele.unique()
wavs = wavs[wavs.allele.isin(common_alleles)]

fig, ax = plt.subplots(2)

sns.barplot(
    data = wavs,
    x = 'allele',
    y = 'allele_freq',
    hue = 'country',
    ax=ax[0]
)

sns.barplot(
    data = wavs,
    x = 'allele',
    y = 'wav',
    hue = 'country',
    ax=ax[1]
)

plt.xticks(rotation=90)
plt.show()
