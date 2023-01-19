"""
CI empirical data
"""

import HLAfreq
from HLAfreq import HLAfreq_pymc as HLAhdi
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import pymc as pm
import numpy as np
import arviz as az
import pandas as pd
import scipy as sp

#                   #
#   Download data   #
#                   #

country = "Thailand"
locus = "DRB1"

country = "Mongolia"
locus = "DQB1"

base_url = HLAfreq.makeURL(country, locus=locus)
aftab = HLAfreq.getAFdata(base_url)
aftab = HLAfreq.only_complete(aftab)
HLAfreq.check_resolution(aftab)
aftab = HLAfreq.decrease_resolution(aftab, 2)
caf = HLAfreq.combineAF(aftab)

HLAfreq.plotAFprob(caf, aftab, xmax=1.1*caf.allele_freq.max())

#                                   #
#   Calculate credible intervals    #
#                                   #

# pymc hdi
hdi = HLAhdi.AFhdi(aftab)
# Dirichlet only ci, too narrow
Dci = HLAfreq.AFci(caf)

# To access the model you have to call the internal functions like so
c_array, allele_names = HLAhdi._make_c_array(aftab)
idata = HLAhdi._fit_Dirichlet_Multinomial(c_array)
hdi = az.hdi(idata).frac.values

# Compare dirichlet and pymc estimates and intervals
shift = .3
plt.scatter(caf.allele_freq, range(caf.shape[0]))
plt.scatter(az.summary(idata)['mean'][:-1], shift+pd.Series(range(caf.shape[0])))
for i in range(caf.shape[0]):
    plt.hlines(
        y = (i),
        xmin = Dci[i][0],
        xmax = Dci[i][1]
    )
    plt.hlines(
        y = i+shift,
        xmin = float(hdi[i][0]),
        xmax = float(hdi[i][1]),
        color="tab:orange"
    )
plt.show()

# Compare sub populations to country
aftab['allele_i'] = ""
for i,allele in enumerate(caf.allele):
    aftab.loc[aftab.allele==allele,'allele_i'] = i

for i in range(caf.shape[0]):
    plt.hlines(
        y = i,
        xmin = float(hdi[i][0]),
        xmax = float(hdi[i][1]),
        color="black"
    )
    plt.hlines(
        y = i,
        xmin = float(Dci[i][0]),
        xmax = float(Dci[i][1]),
        color="red",
    )
    # mask = aftab.allele==caf.allele[i]
    # plt.scatter(
    #     aftab[mask].allele_freq,
    #     [i] * mask.sum()
    # )
for pop in aftab.population.unique():
    pop
    mask = aftab.population == pop
    plt.scatter(
        aftab[mask].allele_freq,
        aftab[mask].allele_i
    )
plt.scatter(caf.allele_freq, range(caf.shape[0]), color='black')
plt.show()
