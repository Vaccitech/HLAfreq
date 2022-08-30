"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA, UMAP, t-SNE). Then plot
those points coloured by geographic area.

For this each country must have an estimate for the same set
of alleles.
"""

import os
import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

regions = pd.read_csv("data/example/globalPCA/countries.csv")

# countries = [
#     'United+Kingdom',
#     'Thailand',
#     'Uganda',
#     'Greece',
#     'Japan',
#     'Germany',
#     'Brazil',
#     'Iraq',
#     'Iran',
#     'Israel',
#     'Chile',
#     'Kenya',
#     'Spain',
# ]

countries = regions.Country.tolist()

# Download HLA allele frequencies
for country in countries:
    print()
    print(country)
    if not os.path.exists("data/example/globalPCA/%s_raw.csv" %country):
        base_url = scrapeAF.makeURL(country, standard="g", locus="A")
        try:
            aftab = scrapeAF.getAFdata(base_url)
            aftab.to_csv("data/example/globalPCA/%s_raw.csv" %country, index=False)
        except:
            print(f"Failed to get data for {country}: {base_url}")

wavs = []
for country in countries:
    try:
        df = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
        wav = scrapeAF.combineAF(df)
        wav['country'] = country
        wavs.append(wav)
    except:
        pass

wavs = pd.concat(wavs, axis=0).reset_index(drop=True)

# Give record for all alleles to all countries
wavs = scrapeAF.unmeasured_alleles(wavs, 'country')

# Filter to single loci
wavsa = wavs[wavs.loci == 'A']

# Check all countries have the same number of alleles
wavsa.groupby('country').allele.unique().apply(len).unique()

# Sort by allele then select a single country,
# check that its alleles match a specified list
# then use the frequencies as features
wavsa = wavsa.sort_values('allele')

sorted_alleles = sorted(list(wavsa.allele.unique()))

AFeatures = []
for country in countries:
    mask = wavsa['country'] == country
    # Only generate features for countries with records
    # e.g. Thailand has no loci A records
    if any(mask):
        assert all(wavsa[mask].allele.unique() == sorted_alleles)
        features = [country] + wavsa[mask].allele_freq.tolist()
        AFeatures.append(features)

# Dataframe of allele frequencies read for dimension reduction
AFeatures = pd.DataFrame(AFeatures, columns = ['country'] + sorted_alleles)
AFeatures.to_csv("data/example/globalPCA/AF_features.csv", index=False)
