"""
Global PCA

Calculate average AF for ~all countries and perform dimension
reduction on that data to 2D (PCA, UMAP, t-SNE). Then plot
those points coloured by geographic area.

For this each country must have an estimate for the same set
of alleles.
"""

# Use umap conda env

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import umap

AFeatures = pd.read_csv("data/example/globalPCA/AF_features.csv")
regions = pd.read_csv("data/example/countries.csv")

# Remove Guinea because it includes data from
# Guinea-Bissau and Papua New Guinea
AFeatures = AFeatures[AFeatures.country != "Guinea"]

pca = PCA(n_components=2)
pca_embedding = pca.fit_transform(AFeatures.loc[:,'A*01:01':].values)
AFeatures['pca0'] = pca_embedding[:,0]
AFeatures['pca1'] = pca_embedding[:,1]

# umapper = umap.UMAP(transform_seed=1234)
# umap_embedding = umapper.fit_transform(AFeatures.loc[:,'A*01:01':].values)
# AFeatures['umap0'] = umap_embedding[:,0]
# AFeatures['umap1'] = umap_embedding[:,1]

# Add region to countries
AFeatures = pd.merge(AFeatures, regions, how="left", left_on="country", right_on="Country")

# Select Venezuela With Prior
vwp = AFeatures[AFeatures.country == "Venezuela with prior"]
# Drop countries without largeRegion, this is Venezuela with a prior
AFeatures = AFeatures.dropna(subset=['largeRegion'])

# plt.style.use('fivethirtyeight')

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
plt.savefig('writeup/globalPCA.png')
plt.close()

# for region in AFeatures.largeRegion.unique():
#     mask = AFeatures.largeRegion == region
#     plt.scatter(AFeatures.umap0[mask], AFeatures.umap1[mask], label=region)
#     for i in AFeatures[mask].index:
#         plt.annotate(AFeatures.loc[i].country.replace("+", " "), AFeatures.loc[i][['umap0','umap1']])
# plt.legend()
# plt.show()
