"""
Entropy of country average allele frequencies
"""

# use maps conda env

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import geoplot as gplt
from scipy.stats import dirichlet

regions = pd.read_csv("data/example/countries.csv")
regions.index = regions.Country
regions['country'] = regions.Country.str.replace("+"," ")

regions['entropy'] = np.nan
regions['total_sample_size'] = np.nan
regions['populations'] = np.nan

countries = regions.Country.tolist()

for country in countries:
    print(country)
    try:
        df = pd.read_csv("data/example/globalPCA/%s_raw.csv" %country)
        df = scrapeAF.decrease_resolution(df, 2)
        caf = scrapeAF.combineAF(df)
        # Recreate the Dirichlet distribution
        alpha = scrapeAF.default_prior(len(caf.allele))
        dd = dirichlet(alpha + caf.c)
        # Calculate entropy
        regions.loc[country, 'entropy'] = dd.entropy()
        # How many populations in the country?
        regions.loc[country, 'populations'] = len(df.population.unique())
        # Calculate total sample size
        regions.loc[country, 'total_sample_size'] = caf.sample_size.unique()[0]
    except:
        pass

######## World plot
world = gpd.read_file(
    gpd.datasets.get_path('naturalearth_lowres')
)

# Not all countries are in both datasets,
# Some of this is due to alternative naming
# Will require manual editting.
regions[~regions.country.isin(world.name)].country
world.name[~world.name.isin(regions.country)]

world_regions = pd.merge(world, regions, how="left", left_on="name", right_on="country")

fig, ax = plt.subplots(2,2)
gplt.choropleth(world_regions, hue="entropy", legend=True, ax=ax[0,0])
ax[0,0].set_title('Entropy')
gplt.choropleth(world_regions, hue="populations", legend=True, ax=ax[0,1])
ax[0,1].set_title('Populations')
gplt.choropleth(world_regions, hue="total_sample_size", legend=True, ax=ax[1,0])
ax[1,0].set_title('Total sample size')
plt.show()

world_regions[world_regions.name.str.contains('ahara')]
world_regions.sort_values('entropy')