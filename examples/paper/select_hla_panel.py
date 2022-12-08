"""
Some countries are not well covered by the IEDB reference set.
We can select a panel of alleles to better represent these
populations, being sure to select a range of supertypes.
"""

import HLAfreq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# The IEDB reference panel
iedb_ref = ['B*35:01', 'A*30:01', 'A*02:03', 'B*15:01', 'B*57:01', 'A*68:02', 'A*23:01', 'B*53:01', 'A*03:01', 'B*40:01', 'B*44:02', 'B*51:01', 'A*33:01', 'A*01:01', 'A*68:01', 'A*24:02', 'B*07:02', 'B*08:01', 'A*02:06', 'A*26:01', 'A*11:01', 'A*30:02', 'A*31:01', 'B*58:01', 'B*44:03', 'A*32:01', 'A*02:01']

regions = pd.read_csv("data/example/countries.csv")
countries = regions.Country.tolist()

# Which supertypes are represented
super1 = pd.read_csv("data/HLA1supertypes_Sidney2008.csv")
super1.allele = super1.allele.apply(lambda x: x[:4]+':'+x[4:])

##################
# Download A and B data for all countries to identify 
# Countries poorly covered by the IEDB reference set
# for locus in ['A', 'B']:
#     # for country in ['Indonesia', 'Philippines']:
#     for country in countries:
#         print(country)
#         try:
#             base_url = HLAfreq.makeURL(country, standard='s', locus=locus)
#             aftab = HLAfreq.getAFdata(base_url)
#             aftab.to_csv(f"data/example/population_coverage/{country}_{locus}_raw.csv", index=False)
#         except:
#             pass

# Calculate total population coverage of IEDB reference
# set (locus A and B) for each country

iedb_cov = []
iedb_country = []

# Load and combine within country data for HLA-A and HLA-B separately
# Calculate the proportion of each population with no alleles
# in the IEDB reference set, assuming no linkage disequilibrium
# between HLA-A and HLA-B
for country in countries:
    try:
        afa = pd.read_csv(f"data/example/population_coverage/{country}_A_raw.csv")
        afa = HLAfreq.only_complete(afa)
        afa = HLAfreq.decrease_resolution(afa, 2)
        cafa = HLAfreq.combineAF(afa)
        maska = cafa.allele.isin(iedb_ref)
        p = cafa[maska].allele_freq.sum()

        afb = pd.read_csv(f"data/example/population_coverage/{country}_B_raw.csv")
        afb = HLAfreq.only_complete(afb)
        afb = HLAfreq.decrease_resolution(afb, 2)
        cafb = HLAfreq.combineAF(afb)
        maskb = cafb.allele.isin(iedb_ref)
        m = cafb[maskb].allele_freq.sum()
        # Population proportion with no IEDB ref alleles
        # at HLA A or HLA B
        coverage = 1 - ((1-p)**2 * (1-m)**2)
        iedb_cov.append(coverage)
        iedb_country.append(country)
    except:
        pass

df = pd.DataFrame({'country':iedb_country, 'coverage':iedb_cov})
df = df.sort_values('coverage')
df = pd.merge(df, regions, how="left", left_on="country", right_on="Country")

# Plot the proportion of each population with no alleles
# in the IEDB reference set, assuming no linkage disequilibrium
# between HLA-A and HLA-B
plt.barh(df.country, df.coverage, color="black")
regions = df.largeRegion.unique()
regions.sort()
for region in regions:
    mask = df.largeRegion == region
    plt.barh(df[mask].country, df[mask].coverage, label=region, zorder=3)
plt.axvline(0.97, linestyle="--", c="black", zorder=5)
plt.tight_layout()
plt.legend(loc="upper left")
plt.grid(zorder=0)
plt.show()

#############
# Plot the cumulative allele frequency of HLA-A and HLA-B
# a single country, with alleles coloured by supertype
country = "Indonesia"
# country="Philippines"
afa = pd.read_csv(f"data/example/population_coverage/{country}_A_raw.csv")
afa = HLAfreq.only_complete(afa)
afa = HLAfreq.decrease_resolution(afa, 2)
cafa = HLAfreq.combineAF(afa)
maska = cafa.allele.isin(iedb_ref)
p = cafa[maska].allele_freq.sum()

afb = pd.read_csv(f"data/example/population_coverage/{country}_B_raw.csv")
afb = HLAfreq.only_complete(afb)
afb = HLAfreq.decrease_resolution(afb, 2)
cafb = HLAfreq.combineAF(afb)
maskb = cafb.allele.isin(iedb_ref)
m = cafb[maskb].allele_freq.sum()
# Population proportion with no IEDB ref alleles
# at HLA A or HLA B
coverage = 1 - ((1-p)**2 * (1-m)**2)
print(coverage)

cafa = cafa.sort_values("allele_freq", ascending=False, ignore_index=True)
cafb = cafb.sort_values("allele_freq", ascending=False, ignore_index=True)
# Add HLA1 supertypes
cafa = pd.merge(cafa, super1, how="left", on="allele")
cafb = pd.merge(cafb, super1, how="left", on="allele")

fig, ax = plt.subplots(1,2, gridspec_kw={'width_ratios': [1,2]})
ax[0].plot(cafa.allele,
    cafa.allele_freq.cumsum().apply(HLAfreq.population_coverage),
    c="black",
    zorder=0)
for supertype in cafa.supertype.unique():
    mask = cafa.supertype == supertype
    ax[0].scatter(
        cafa[mask].allele,
        cafa.allele_freq.cumsum().apply(HLAfreq.population_coverage)[mask],
        label=supertype,
        zorder=10)
ax[0].set_title("HLA A")
ax[0].set_ylabel("Cumulative population coverage")
ax[0].grid(True)
ax[0].legend()
ax[0].tick_params(labelrotation=90)

ax[1].plot(cafb.allele,
    cafb.allele_freq.cumsum().apply(HLAfreq.population_coverage),
    c="black",
    zorder=0)
for supertype in cafb.supertype.unique():
    mask = cafb.supertype == supertype
    ax[1].scatter(
        cafb[mask].allele,
        cafb.allele_freq.cumsum().apply(HLAfreq.population_coverage)[mask],
        label=supertype,
        zorder=10)
ax[1].set_title("HLA B")
ax[1].grid(True)
ax[1].legend()
ax[1].tick_params(labelrotation=90)

plt.tight_layout()
plt.show()

#########################
#                       #
#   Design a new panel  #
#                       #
#########################
# When selecting a panel of HLA alleles to represent a country,
# consider different genes,
# and different supertypes,
# not only the most common alleles.