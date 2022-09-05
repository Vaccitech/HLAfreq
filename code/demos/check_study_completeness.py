"""
Check that study-loci allele frequency sums to ~1.
Some studies report alleles that do not sum to 1.
Also download errors can lead to incomplete datasets.
Some download settings can also lead to unexpected results.
E.g. some studies report alleles to mixed resolution, therefore;
if you search for a specific allele resolution you will download
an incomplete dataset.

Having the complete dataset is important, especially when
calculating the Dirichlet distribution as reported alleles
will be over weighted if some alleles are unreported.

How to identify and resolve allele resolution. Because allele
frequency estimates can only be combined if the allele name is
identical all alleles in a study must have the same resolution
i.e. A*01 is 2 digit resolution (or 1 field) and A*01:01 is
4 digit (or 2 field). These cannot be combined as we cannot
tell whether A*01 is A*01:01, A*01:02 etc. For the same reason
we can decrease the resolution of an allele but we cannot
increase it.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

country = 'Thailand'

# Because resolution equals 2, an incomplete dataset is returned
base_url = scrapeAF.makeURL(country, resolution_pattern="equal", resolution=2, standard="s")
aftab = scrapeAF.getAFdata(base_url)
aftab.to_csv("data/example/incomplete_raw.csv", index=False)

df = pd.read_csv("data/example/incomplete_raw.csv")

# How to filter incomplete studies

# Identify incomplete studies
noncomplete = scrapeAF.incomplete_studies(df)

# Returns False if population AND loci are in the noncomplete.index
# AS A PAIR
# This is important so that we don't throw away all data on a population
# just because one loci is incomplete.
complete_mask = df.apply(lambda x: (x.population, x.loci) not in noncomplete.index, axis=1)

df = df[complete_mask]

# Check filtering success
scrapeAF.incomplete_studies(df)

#################################
#                               #
#  Differing allele resolution  #
#                               #
#################################

base_url = scrapeAF.makeURL(
    'Thailand', locus="A",
    resolution_pattern="bigger_equal_than", resolution=1
)
aftab = scrapeAF.getAFdata(base_url)
aftab.to_csv("data/example/multiresolution_raw.csv", index=False)

df = pd.read_csv("data/example/multiresolution_raw.csv")
scrapeAF.check_resolution(df)

# How many fields are reported for each allele?
1 + df.allele.str.count(":")
# How many alleles report to exactly n fields?
(1 + df.allele.str.count(":")).value_counts()

# You cant increase resolution 
scrapeAF.decrease_resolution(df, 2)

# But you can decrease
df = scrapeAF.decrease_resolution(df, 1)
scrapeAF.check_resolution(df)

# In principal you could use 2 digit resolution
# to estimate 4 digit resolution. E.g if A*01 has
# frequency p we could increase the resolution by
# spreading this across all possible A*01:xx alleles
# More formally for k A*:01:xx alleles they would have
# frequency p/k. This would effectively update the
# Dirichlet distribution with a flat prior about which
# allele A*01 was actually measured.

# If this is of interest to you please open an issue or
# pull request on github.

