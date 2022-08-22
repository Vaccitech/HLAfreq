"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA, UMAP, t-SNE). Then plot
those points coloured by geographic area.

For this each country must have an estimate for the same set
of alleles.
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

# # Download HLA allele frequencies
# for country in countries:
#     print(country)
#     base_url = scrapeAF.makeURL(country)
#     aftab = scrapeAF.getAFdata(base_url)
#     aftab.to_csv("data/example/%s_raw.csv" %country, index=False)

wavs = []
for country in countries:
    df = pd.read_csv("data/example/%s_raw.csv" %country)
    wav = scrapeAF.combineAF(df)
    wav['country'] = country
    wavs.append(wav)

wavs = pd.concat(wavs, axis=0).reset_index()

scrapeAF.unmeasured_alleles(wavs, 'country')