"""
Download HLA data for all countries and calculate
population coverage of IEDB set of 27 alleles.

Plot population coverage in detail for Indonesia.

Takes ~1hour with data downloads
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import HLAfreq
import requests

# Create folder for data
try:
    os.makedirs("data/example/population_coverage")
except FileExistsError:
    pass

# Population coverage of IEDB 27 alleles
# The IEDB reference panel
iedb_ref = ['B*35:01', 'A*30:01', 'A*02:03', 'B*15:01', 'B*57:01', 'A*68:02', 'A*23:01', 'B*53:01', 'A*03:01', 'B*40:01', 'B*44:02', 'B*51:01', 'A*33:01', 'A*01:01', 'A*68:01', 'A*24:02', 'B*07:02', 'B*08:01', 'A*02:06', 'A*26:01', 'A*11:01', 'A*30:02', 'A*31:01', 'B*58:01', 'B*44:03', 'A*32:01', 'A*02:01']
iedb_refa = [i for i in iedb_ref if "A" in i]
iedb_refb = [i for i in iedb_ref if "B" in i]

# Download countries in regions as defined on 
# http://www.allelefrequencies.net/datasets.asp#tag_4
r = requests.get("https://raw.githubusercontent.com/Vaccitech/HLAfreq/main/data/example/countries.csv")
with open("data/example/countries.csv", "w") as f:
    f.write(r.text)

regions = pd.read_csv("data/example/countries.csv")
countries = regions.Country.tolist()


# Download A and B data for all countries to identify 
# Countries poorly covered by the IEDB reference set
for locus in ['A', 'B']:
    # for country in ['Indonesia', 'Philippines']:
    for country in countries:
        print(country)
        try:
            base_url = HLAfreq.makeURL(country, standard='s', locus=locus)
            aftab = HLAfreq.getAFdata(base_url)
            aftab.to_csv(f"data/example/population_coverage/{country}_{locus}_raw.csv", index=False)
        except:
            pass


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

##############
# Plot cumulative frequency of IEDB panel
# vs
# country specific alleles
#A
country = "Indonesia"
afa = pd.read_csv(f"data/example/population_coverage/{country}_A_raw.csv")
afa = HLAfreq.only_complete(afa)
afa = HLAfreq.decrease_resolution(afa, 2)
cafa = HLAfreq.combineAF(afa)

cafa = cafa.sort_values('allele_freq', ascending=False, ignore_index=True)
plt.scatter(cafa.allele, cafa.allele_freq.cumsum().apply(HLAfreq.population_coverage))
plt.grid(True); plt.xticks(rotation=90); plt.tight_layout()
plt.xlabel('Allele'); plt.ylabel('Cumulative population coverage')
# plt.tight_layout(); plt.show()

maska = cafa.allele.isin(iedb_refa)
cafai = cafa[maska].sort_values('allele_freq', ascending=False, ignore_index=True)
# Add unobserved IEDB ref alleles
missing_iedb = pd.DataFrame([[i, 'A',0,0,473,1,0] for i in iedb_refa if not i in cafai.allele.values], columns=cafai.columns)
cafai = pd.concat([cafai, missing_iedb])
plt.scatter(cafai.allele, cafai.allele_freq.cumsum().apply(HLAfreq.population_coverage))
plt.grid(True); plt.xticks(rotation=90); plt.tight_layout()
plt.xlabel('Allele'); plt.ylabel('Cumulative population coverage')
plt.ylim(0,1.05)
plt.tight_layout(); plt.show()
