"""
Visualise high k Dirichlet
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import dirichlet, beta
import pandas as pd
import code.scrapeAF as scrapeAF

#####################
#                   #
#   Define alpha    #
#                   #
#####################

# Load allele frequency data
AFtab = pd.read_csv(f"data/example/outlier_raw.csv")
# Reduce resolution to 2 fields and collapse alleles
AFtab = scrapeAF.decrease_resolution(AFtab, 2)
# Combined allle frequency
# noncomplete = scrapeAF.incomplete_studies(AFtab)
# complete_mask = AFtab.apply(lambda x: (x.population, x.loci) not in noncomplete.index, axis=1)
# AFtab = AFtab[complete_mask]
caf = scrapeAF.combineAF(AFtab)
# Recreate the Dirichlet distribution
naivealpha = scrapeAF.default_prior(len(caf.allele))

# Dirichlet distribution concentration parameters 
# based on allele frequency data
alpha = naivealpha + caf.c

######################
#                    #
# View k=3 Dirichlet #
#                    #
######################
if len(alpha) == 3:
    dd = dirichlet(alpha)
    # Make data for k=3.
    k1 = np.arange(0, 1, 0.05)
    k2 = np.arange(0, 1, 0.05)
    k1, k2 = np.meshgrid(k1, k2)
    k3 = 1 - k1 - k2
    k3[k3<0] = np.nan

    pdf = np.zeros(k1.shape)
    for i in range(k1.shape[0]):
        for j in range(k1.shape[1]):
            pdf[i,j] =  dd.pdf([k1[i,j], k2[i,j], k3[i,j]])

    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(k1, k2, pdf, cmap=cm.plasma,
        linewidth=0, antialiased=False)
    ax.set_title('PDF')
else:
    print( "Surface plot only works with k==3")

######################
#                    #
#   View Dirichlet   #
#       as beta      #
#                    #
######################

# View pdf only
scrapeAF.plotAFprob(concentration=alpha, ci=0)
# View estimated AF
scrapeAF.plotAFprob(caf)
# View estimated AF with individual studies
scrapeAF.plotAFprob(caf, AFtab)
# View subset of alleles
scrapeAF.plotAFprob(caf, AFtab, alleles=["DQB1*03:01","DQB1*05:01","DQB1*05:02"])

# CI may be misleading

# Use cafs from multi_country.py
x = cafs[cafs.allele == 'A*23:01']

clist = 1 + x.allele_freq * x.sample_size
long_c = []
for c in clist:
    long_c.append(c)
    long_c.append(0)
scrapeAF.plotAFprob(concentration=pd.Series(long_c), ncol=2)

pline = np.linspace(0,1,500)
fig, ax = plt.subplots()
for a,b in zip(x.c, x.sample_size):
    bd = beta(a,b)
    pdf = [bd.pdf(p) for p in pline]
    ax.plot(pline, pdf)
    ax.set_xlim(0,0.4)
plt.show()
plt.subplots()