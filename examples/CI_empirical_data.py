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

#                   #
#   Calculate CI    #
#                   #

hdi = HLAhdi.AFhdi(aftab)
# Add to caf
caf = pd.merge(caf, hdi, 'left', on='allele')
caf[['allele', 'allele_freq', 'post_mean', 'lo', 'hi']]

#           #
#   Plot CI #
#           #

HLAfreq.plotAF(caf, AFtab=aftab, hdi=hdi, compound_mean=hdi)

#                           #
#   Access model directly   #
#                           #

c_array, allele_names = HLAhdi._make_c_array(aftab)
idata = HLAhdi._fit_Dirichlet_Multinomial(c_array)
idata.posterior.frac

az.plot_trace(idata); plt.show()