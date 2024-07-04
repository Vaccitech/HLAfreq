"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA). Then plot
those points coloured by geographic area.

This script requires scikit-learn which is not included
with HLAFreq.

Takes ~30 minutes to download data
"""

import os
import requests
import HLAfreq as HLAfreq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Create folder for data
try:
    os.makedirs("data/example/globalPCA")
except FileExistsError:
    pass

# Download countries in regions as defined on 
# http://www.allelefrequencies.net/datasets.asp#tag_4
r = requests.get("https://raw.githubusercontent.com/Vaccitech/HLAfreq/main/data/example/countries.csv")
with open("data/example/countries.csv", "w") as f:
    f.write(r.text)

regions = pd.read_csv("data/example/countries.csv")

countries = regions.Country.tolist()

############
# Download data
############
# Download HLA allele frequencies
# Not all countries have data
# those without will print "Failed to get data for..."
for country in countries:
    print()
    print(country)
    if not os.path.exists("data/example/globalPCA/%s_raw.csv" %country):
        base_url = HLAfreq.makeURL(country, standard="g", locus="A")
        try:
            aftab = HLAfreq.getAFdata(base_url)
            aftab.to_csv("data/example/globalPCA/%s_raw.csv" %country, index=False)
        except:
            pass

############
# Using prior for Venezuela
############
# Load Venezuela data
country = "Venezuela"
venezuelaAF = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
venezuelaAF = HLAfreq.only_complete(venezuelaAF)
venezuelaAF = HLAfreq.decrease_resolution(venezuelaAF, 2)
cafV = HLAfreq.combineAF(venezuelaAF)

# Load Colombia data
country = "Colombia"
colombiaAF = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
colombiaAF = HLAfreq.only_complete(colombiaAF)
colombiaAF = HLAfreq.decrease_resolution(colombiaAF, 2)
cafC = HLAfreq.combineAF(colombiaAF)
cafC['country'] = country
cafC['population'] = country

# Create a prior from Colombia
study_prior = cafC.copy()
study_prior.sample_size = study_prior.sample_size / 30
# Combine Allele Frequency with Study Prior
cafSP = HLAfreq.combineAF(
    pd.concat([venezuelaAF, study_prior], join="inner")
    )
cafSP['country'] = "Venezuela with prior"

############
# Load each country data
############
cafs = []
for country in countries:
    try:
        df = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
        df = HLAfreq.only_complete(df)
        df = HLAfreq.decrease_resolution(df, 2)
        caf = HLAfreq.combineAF(df)
        caf['country'] = country
        cafs.append(caf)
    except:
        pass

# Add Venezuela with prior to dataset
cafs.append(cafSP)

cafs = pd.concat(cafs, axis=0).reset_index(drop=True)

# Give record for all alleles to all countries
cafs = HLAfreq.unmeasured_alleles(cafs, 'country')

# Check all countries have the same number of alleles
cafs.groupby('country').allele.unique().apply(len).unique()

# Sort by allele then select a single country,
# check that its alleles match a specified list
# then use the frequencies as features
cafs = cafs.sort_values('allele')

sorted_alleles = sorted(list(cafs.allele.unique()))

AFeatures = []
for country in cafs.country.unique():
    mask = cafs['country'] == country
    # Only generate features for countries with records
    # e.g. Thailand has no loci A records
    if any(mask):
        assert all(cafs[mask].allele.unique() == sorted_alleles)
        features = [country] + cafs[mask].allele_freq.tolist()
        AFeatures.append(features)

# Dataframe of allele frequencies read for dimension reduction
AFeatures = pd.DataFrame(AFeatures, columns = ['country'] + sorted_alleles)
AFeatures.to_csv("data/example/globalPCA/AF_features.csv", index=False)

###################
#   Prep plots    #
###################
AFeatures = pd.read_csv("data/example/globalPCA/AF_features.csv")
regions = pd.read_csv("data/example/countries.csv")

# Remove Guinea because it includes data from
# Guinea-Bissau and Papua New Guinea
AFeatures = AFeatures[AFeatures.country != "Guinea"]

pca = PCA(n_components=2)
pca_embedding = pca.fit_transform(AFeatures.loc[:,'A*01:01':].values)
AFeatures['pca0'] = pca_embedding[:,0]
AFeatures['pca1'] = pca_embedding[:,1]

# Add region to countries
AFeatures = pd.merge(AFeatures, regions, how="left", left_on="country", right_on="Country")

# Select Venezuela With Prior
vwp = AFeatures[AFeatures.country == "Venezuela with prior"]
# Drop countries without largeRegion, this is Venezuela with a prior
AFeatures = AFeatures.dropna(subset=['largeRegion'])

# Define offset values for each country label
offset = dict(zip(AFeatures.country, [(0,0)]*AFeatures.shape[0]))
# Customise offsets
offset['Zimbabwe'] = (0,5)
offset['United+States'] = (0,5)
offset['Georgia'] = (0,5)
offset['Germany'] = (0,5)
offset['Italy'] = (0,-5)
offset['Czech+Republic'] = (0,-4)
offset['Poland'] = (0,-7)
offset['France'] = (0,-2)
offset['Bosnia+and+Herzegovina'] = (0,4)
offset['Russia'] = (0,-5)
offset['Ireland'] = (0,-5)
offset['Netherlands'] = (0,-6)
offset['Austria'] = (0,-7)
offset['Romania'] = (0,4)
offset['Wales'] = (0,-5)
offset['Cameroon'] = (0,2)
offset['Ghana'] = (0,-3)
offset['United+Kingdom'] = (0,-4)
offset['Sao+Tome+and+Principe'] = (0,-5)

# Plot PCA
regions = AFeatures.largeRegion.unique()
regions.sort()
for region in regions:
    mask = AFeatures.largeRegion == region
    plt.scatter(AFeatures.pca0[mask], AFeatures.pca1[mask], label=region)
    for i in AFeatures[mask].index:
        # Label each country
        plt.annotate(
            AFeatures.loc[i].country.replace("+", " "),
            AFeatures.loc[i][['pca0','pca1']],
            offset[AFeatures.loc[i].country],
            textcoords="offset points",
            fontsize=8)

plt.legend()
# Arrow shows change of Venezuela with prior
plt.annotate(
    "",
    xy=(vwp.pca0, vwp.pca1),
    xytext=(AFeatures[AFeatures.country == "Venezuela"].pca0, AFeatures[AFeatures.country == "Venezuela"].pca1),
    arrowprops=dict(arrowstyle="->")
)
# Open point to show updated Venezuela
plt.scatter(vwp.pca0, vwp.pca1, color="tab:blue", facecolors="none")
plt.xlabel('1st principal component'); plt.ylabel('2nd principal component')
plt.gcf().set_size_inches(10,10)
plt.show()
