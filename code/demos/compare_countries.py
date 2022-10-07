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
    'Uganda'
]

# Download HLA allele frequencies
for country in countries:
    print(country)
    base_url = scrapeAF.makeURL(country)
    aftab = scrapeAF.getAFdata(base_url)
    aftab.to_csv("data/example/%s_raw.csv" %country, index=False)
# If you use the data directly without saving and loading it
# pandas treats the numeric columns as strings (allele_freq and sample size)
# which throws errors later when doing maths

wavs = []
for country in countries:
    df = pd.read_csv("data/example/%s_raw.csv" %country)
    df = df[df.loci=="A"]
    wav = scrapeAF.combineAF(df)
    # Add credible intervals
    wav['cl','cu'] = scrapeAF.AFci(wav)
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

# Similarity of weighted allele freq and Diriclet mean
x = max(max(wavs.allele_freq),max(wavs.wav))
plt.plot((0,x), (0,x), 'k--')
plt.scatter(
    wavs.allele_freq,
    wavs.wav
)
plt.show()
