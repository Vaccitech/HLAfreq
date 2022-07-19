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
    # Handle commmas in sample size
    if df.sample_size.dtype == 'O':
        df.sample_size = pd.to_numeric(df.sample_size.str.replace(",",""))
    # Weighted AVerage
    wav = df.groupby('allele').apply(lambda x: np.average(x.allele_freq, weights=x.sample_size))
    wav.name = "allele_freq"
    wav = pd.DataFrame(wav)
    wav['country'] = country
    wavs.append(wav)

wavs = pd.concat(wavs, axis=0).reset_index()

# Keep only alleles with frequency over threshold
threshold = 0.02
common_alleles = wavs[wavs.allele_freq > threshold].allele.unique()
wavs = wavs[wavs.allele.isin(common_alleles)]

# Focus on locus A and locus B
wavsa = wavs[wavs.allele.str.contains("^A\*")]
wavsb = wavs[wavs.allele.str.contains("^B\*")]

fig, ax = plt.subplots(2)

sns.barplot(
    data = wavsa,
    x = 'allele',
    y = 'allele_freq',
    hue = 'country',
    ax=ax[0]
)
sns.barplot(
    data = wavsb,
    x = 'allele',
    y = 'allele_freq',
    hue = 'country',
    ax=ax[1]
)

plt.xticks(rotation=90)
plt.show()