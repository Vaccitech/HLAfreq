"""
Identify poulations that are outliers relative to the country.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import dirichlet

# Dirichlet throws an error if frequencies do not sum to 1.
# Massage the frequencies to sum to 2
def sumto1(af):
    # Compress/Inflate all values by a factor so they sum to 1
    ciFactor = 1/np.sum(af)
    af1 = [x*ciFactor for x in af]
    return af1

country = "Thailand"
# Download dataset
# base_url = scrapeAF.makeURL(country, 'g', locus='DQB1')
# AFtab = scrapeAF.getAFdata(base_url)
# AFtab.to_csv("data/example/outlier_raw.csv", index=False)

AFtab = pd.read_csv(f"data/example/outlier_raw.csv")
# Reduce resolution to 2 fields and collapse alleles
AFtab = scrapeAF.decrease_resolution(AFtab, 2)
# Combined allle frequency
caf = scrapeAF.combineAF(AFtab)

# Recreate the Dirichlet distribution
alpha = scrapeAF.default_prior(len(caf.allele))
dd = dirichlet(alpha + caf.c)

# Get allele freqs from a single population
# it must be in the same order as the dirichlet
datasetID = 'population'
df = scrapeAF.unmeasured_alleles(AFtab, datasetID=datasetID)
df = df.sort_values('allele')

df.groupby('population').allele.apply(list).apply(lambda x: x == caf.allele.tolist())
df.groupby('population').allele_freq.apply(list)

# The probability density function of actual samples is almost always zero
# This may be because some alleles are never seen in some samples
# This could be solved by adding pseudo counts...
for i in range(5):
    # print(i)
    af = df.groupby('population').allele_freq.apply(list)[i]
    # pseudoAF = [x if x>0 else x+0.1 for x in af]
    # af1 = sumto1(pseudoAF)
    af1 = sumto1(af)
    # [round(x,3) for x in af1]
    dd.logpdf(af1)

[round(x,3) for x in caf.allele_freq.tolist()]

# Random draws from the Dirichlet distribution
# do not have zero pdf
draws = pd.DataFrame(np.random.dirichlet(alpha+caf.c, 1000))
pdfs = draws.apply(lambda x: dd.pdf(list(x)), axis=1)
plt.hist(pdfs)
plt.show()



# Should we have a function/class that creates
# and holds the dirichlet distribution?