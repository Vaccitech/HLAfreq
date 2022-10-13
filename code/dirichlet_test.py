"""
Tests with the dirichlet distribution
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import dirichlet, beta
import pandas as pd
import code.HLAfreq as HLAfreq

# define a Dirichlet Distributions with parameters alpha
alpha = np.array([0.4, 5, 15])  # specify concentration parameters
dd = dirichlet(alpha)

dd.mean()
dd.var()
dd.entropy()

dd.pdf([.6,.3,.1])

# Make data for k=3.
k1 = np.arange(0, 1, 0.05)
k2 = np.arange(0, 1, 0.05)
k1, k2 = np.meshgrid(k1, k2)
k3 = 1 - k1 - k2
k3[k3<0] = np.nan

################################
#                              #
#   Compare pdf and log pdf    #
#                              #
################################
alpha = (3,3,3)

pdf = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        pdf[i,j] =  dirichlet(alpha).pdf([k1[i,j], k2[i,j], k3[i,j]])

logpdf = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        logpdf[i,j] =  dirichlet(alpha).logpdf([k1[i,j], k2[i,j], k3[i,j]])

# Plot the surface.
fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={"projection": "3d"})
ax[0].plot_surface(k1, k2, pdf, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[0].set_title('PDF')
ax[1].plot_surface(k1, k2, logpdf, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[1].set_title('Log pdf')
plt.show()

################################
#                              #
#       Update dirichlet       #
#                              #
################################

# update a dirichlet distribution with new info
# plot pdf before and after
alpha = np.array([1,1,1])
c = np.array([2,0,1])
alpha2 = alpha+c

pdf = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        pdf[i,j] =  dirichlet(alpha).pdf([k1[i,j], k2[i,j], k3[i,j]])

pdf2 = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        pdf2[i,j] =  dirichlet(alpha2).pdf([k1[i,j], k2[i,j], k3[i,j]])

# Plot the surface.
fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={"projection": "3d"})
ax[0].plot_surface(k1, k2, pdf, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[0].set_title(f'PDF {alpha}')
ax[1].plot_surface(k1, k2, pdf2, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[1].set_title(f'PDF {alpha2}')
plt.show()

############################
#                          #
#   Study error ignored    #
#                          #
############################

# Does not account for study level error
alpha = np.array([1,1,1])
# Two studies in agreement
agree = alpha + np.array([3,3,0]) + np.array([3,3,0])
#Two studies disagreeing
disagree = alpha + np.array([6,0,0]) + np.array([0,6,0])

all(agree == disagree)

pdf = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        pdf[i,j] =  dirichlet(agree).pdf([k1[i,j], k2[i,j], k3[i,j]])

pdf2 = np.zeros(k1.shape)
for i in range(k1.shape[0]):
    for j in range(k1.shape[1]):
        pdf2[i,j] =  dirichlet(disagree).pdf([k1[i,j], k2[i,j], k3[i,j]])

# Plot the surface.
fig, ax = plt.subplots(nrows=1, ncols=2, subplot_kw={"projection": "3d"})
ax[0].plot_surface(k1, k2, pdf, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[0].set_title(f'PDF agree')
ax[1].plot_surface(k1, k2, pdf2, cmap=cm.plasma,
    linewidth=0, antialiased=False)
ax[1].set_title(f'PDF disagree')
plt.show()

#############################
#                           #
#   Dirichlet allele freq   #
#                           #
#############################

country = "Thailand"
df = pd.read_csv("data/example/%s_raw.csv" %country)
loci = "DPB1"
dfa = df[df.loci == loci]
caf = HLAfreq.combineAF(dfa)

plt.scatter(
    caf.wav,
    caf.allele_freq
)
plt.show()

plt.scatter(
    caf.wav,
    caf.wav-caf.allele_freq
)
plt.show()

########################
#                      #
#   Central credible   #
#   region from beta   #
#                      #
########################

pline = np.linspace(0, 1, 1000)
a,b = 2,5
bd = beta(a,b)
pdf = bd.pdf(pline)
cl,cu = HLAfreq.betaCI(a,b)

plt.plot(pline, pdf)
plt.fill_between([cl,cu], [bd.pdf(cl), bd.pdf(cu)], alpha=0.5)
plt.show()
