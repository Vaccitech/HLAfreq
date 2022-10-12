"""
Single country

Calculate average
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#########################
#                       #
#   Download AF data    #
#                       #
#########################
"""
Pick a country to download AF data for. Spaces in country names
are replaced with '+', if in doubt use the allelefrequencies.net
website to search for a country and then check the url for the
appropriate spelling.
"""

country = "Thailand"

base_url = scrapeAF.makeURL(country)
aftab = scrapeAF.getAFdata(base_url)
aftab.to_csv("data/example/%s_raw.csv" %country, index=False)
#aftab = pd.read_csv("data/example/%s_raw.csv" %country)

#########################
#                       #
#   Check completeness  #
#                       #
#########################
"""
We can't yet combine allele frequency estimates because not all
studies sum to an allele frequency of 1. In this case it is because some studies
reported allele frequencies at 1 field of resolution and we only downloaded
alleles with >2.

To drop incomplete studies use `only_complete()` which will print out the studies
that were incomplete. To view incomplete studies without dropping them use
`incomplete_studies()` before dropping studies.
"""

# Drop any incomplete studies
aftab = scrapeAF.only_complete(aftab)

#########################
#                       #
#   Check resolution    #
#                       #
#########################
"""
Because allele frequency estimates can only be combined if the
allele name is identical, all alleles in a set of studies must have
the same resolution i.e. A*01 is 2 digit resolution (or 1 field) and
A*01:01 is 4 digit (or 2 field). These cannot be combined as we
cannot tell whether A*01 is A*01:01, A*01:02 etc. For the same reason
we can decrease the resolution of an allele but we cannot
increase it.

Before combining allele frequency estimates we may have to reduce the
resolution of some alleles. `check_resolution()` will report if all
alleles are at the same resolution and the number of estimates at
each resolution.

`decrease_resolution()` reduces the resolution of all alleles to a
specified level so they can be combined.
"""

scrapeAF.check_resolution(aftab)

aftab = scrapeAF.decrease_resolution(aftab, 2)


#########################
#                       #
#     Combine AFs       #
#                       #
#########################

"""
We can only combine allele frequency estimates for a single loci
at a time. So we filter the data to a single loci and then combine
allele frequency estimates using `combineAF()`
"""

scrapeAF.combineAF(aftab)

afloc = aftab[aftab.loci=="DRB1"]
caf = scrapeAF.combineAF(afloc)

# Calculate credible interval
caf[['cl','cu']] = scrapeAF.AFci(caf)

caf.plot.barh('allele', 'allele_freq')
plt.tight_layout()
plt.show()

scrapeAF.plotAFprob(caf, afloc)
scrapeAF.plotAFprob(caf, afloc, alleles=['DRB1*09:01', 'DRB1*15:02', 'DRB1*07:01', 'DRB1*12:02',])

#############################
#                           #
#    Population coverage    #
#                           #
#############################

# Sort by allele frequency
caf = caf.sort_values('allele_freq', ascending=False, ignore_index=True)

plt.scatter(caf.allele, caf.allele_freq.cumsum().apply(scrapeAF.population_coverage))
plt.plot(caf.allele_freq.cumsum().apply(scrapeAF.population_coverage), label="Cumulative coverage")
plt.grid(True)
plt.xticks(rotation=90)
plt.legend()
plt.tight_layout()
plt.show()
