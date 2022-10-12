"""
Some countries are not well covered by the IEDB reference set.
We can select a panel of alleles to better represent these
populations, being sure to select a range of supertypes.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
#             base_url = scrapeAF.makeURL(country, standard='s', locus=locus)
#             aftab = scrapeAF.getAFdata(base_url)
#             aftab.to_csv(f"data/example/population_coverage/{country}_{locus}_raw.csv", index=False)
#         except:
#             pass

# Calculate total population coverage of IEDB reference
# set (locus A and B) for each country

iedb_cov = []
iedb_country = []

for country in countries:
    try:
        afa = pd.read_csv(f"data/example/population_coverage/{country}_A_raw.csv")
        afa = scrapeAF.only_complete(afa)
        afa = scrapeAF.decrease_resolution(afa, 2)
        cafa = scrapeAF.combineAF(afa)
        maska = cafa.allele.isin(iedb_ref)
        p = cafa[maska].allele_freq.sum()

        afb = pd.read_csv(f"data/example/population_coverage/{country}_B_raw.csv")
        afb = scrapeAF.only_complete(afb)
        afb = scrapeAF.decrease_resolution(afb, 2)
        cafb = scrapeAF.combineAF(afb)
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

plt.barh(df.country, df.coverage, color="black")
plt.barh()
for region in df.largeRegion.unique():
    mask = df.largeRegion == region
    plt.barh(df[mask].country, df[mask].coverage, label=region)
plt.axvline(0.95, linestyle="--", c="black")
plt.legend()
plt.show()

#############
country = "Indonesia"
# country="Philippines"
afa = pd.read_csv(f"data/example/population_coverage/{country}_A_raw.csv")
afa = scrapeAF.only_complete(afa)
afa = scrapeAF.decrease_resolution(afa, 2)
cafa = scrapeAF.combineAF(afa)
maska = cafa.allele.isin(iedb_ref)
p = cafa[maska].allele_freq.sum()

afb = pd.read_csv(f"data/example/population_coverage/{country}_B_raw.csv")
afb = scrapeAF.only_complete(afb)
afb = scrapeAF.decrease_resolution(afb, 2)
cafb = scrapeAF.combineAF(afb)
maskb = cafb.allele.isin(iedb_ref)
m = cafb[maskb].allele_freq.sum()
# Population proportion with no IEDB ref alleles
# at HLA A or HLA B
coverage = 1 - ((1-p)**2 * (1-m)**2)

cafa = cafa.sort_values("allele_freq", ascending=False, ignore_index=True)
cafb = cafb.sort_values("allele_freq", ascending=False, ignore_index=True)
# Add HLA1 supertypes
cafa = pd.merge(cafa, super1, how="left", on="allele")
cafb = pd.merge(cafb, super1, how="left", on="allele")

fig, ax = plt.subplots(1,2, gridspec_kw={'width_ratios': [1,2]})
ax[0].plot(cafa.allele,
    cafa.allele_freq.cumsum().apply(scrapeAF.population_coverage),
    c="black",
    zorder=0)
for supertype in cafa.supertype.unique():
    mask = cafa.supertype == supertype
    ax[0].scatter(
        cafa[mask].allele,
        cafa.allele_freq.cumsum().apply(scrapeAF.population_coverage)[mask],
        label=supertype,
        zorder=10)
ax[0].set_title("HLA A")
ax[0].set_ylabel("Cumulative population coverage")
ax[0].grid(True)
ax[0].legend()
ax[0].tick_params(labelrotation=90)

ax[1].plot(cafb.allele,
    cafb.allele_freq.cumsum().apply(scrapeAF.population_coverage),
    c="black",
    zorder=0)
for supertype in cafb.supertype.unique():
    mask = cafb.supertype == supertype
    ax[1].scatter(
        cafb[mask].allele,
        cafb.allele_freq.cumsum().apply(scrapeAF.population_coverage)[mask],
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

def multilocus_homozygous_coverage(AFs):
    coverage = np.prod([p**2 for p in AFs])
    return coverage

def remainingAF(x, cumAF, is_locus):
    return [1-c-x*modifier for c,modifier in zip(cumAF, is_locus)]
    

def population_uncoverage(row, covered_dict):
    loci = np.array(list(covered_dict.keys()))
    is_locus = loci == row.loci
    remains = remainingAF(row.allele_freq, covered_dict.values(), is_locus)
    uncovered = multilocus_homozygous_coverage(remains)
    return uncovered

def sum_panel_AF(panel, locus):
    return sum([x.allele_freq for x in panel if x.loci == locus])

all_loci = pd.concat([cafa, cafb], ignore_index=True)

covered_dict = {
    'A':0,
    'B':0
}
panel = []
mask = [True] *  all_loci.shape[0]

for n in range(4):
    # Select allele
    i = all_loci[mask].apply(population_uncoverage, axis=1, args=[covered_dict]).idxmin()
    # Add selected allele to panel
    panel.append(all_loci.loc[i])
    # Update mask to omited selected loci from choice
    mask = ~all_loci.allele.isin([x.allele for x in panel])
    # Update coverd_dict
    covered_dict = {key: sum_panel_AF(panel, key) for key, value in covered_dict.items()}


panel = pd.concat(panel, axis=1, ignore_index=True).T
panel

# To my surprise it generally favours adding more alleles to which
# ever locus already has the highest cumulative frequency

def f(p,m):
    return (1-p)**2*(1-m)**2
    # return (1-p)*(1-m)

# def f(q,n):
    # return (q**2) * (n**2)
    # return q*n

x = np.linspace(1, 0, 21)
y = np.linspace(1, 0, 21)
z = np.array([f(j,i) for j in y for i in x])
X, Y = np.meshgrid(x, y)
Z = z.reshape(21, 21)

from matplotlib import cm
# Plot the surface.
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()

plt.plot(f(0.2, x))
plt.plot(f(0.8, x))
plt.show()