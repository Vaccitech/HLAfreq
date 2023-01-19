"""
Estimate confidence intervals on simulated data
using pymc and pymc_env
"""

import pymc as pm
import matplotlib.pyplot as plt
import numpy as np
import arviz as az
import pandas as pd
import scipy as sp

def simulateAF(k, concentration, n, samplesize_lo=5, samplesize_hi=200):
    """ Simulate allele frequency dataset for k alleles and n populations.
    The global allele frequency is a random draw from a uniform Dirichlet.
    This is used with the concentration parameter to define a new dirichlet.
    Each of n populations is then an independent draw from this new Dirichlet.
    The observed allele_counts then use these population allele frequencies to
    define a multinomial distribution and take a random number of samples from it.
    The range of posible sample size is between samplesize_lo and samplesize_hi.
    """
    # Global AF
    alpha = sp.stats.dirichlet(np.ones(k)).rvs()
    alpha = alpha.reshape(-1)
    # The Dirichlet distribution which populations are drawn from
    dd = sp.stats.dirichlet(alpha*concentration)
    # Population AF
    popAF =  dd.rvs(n)
    sample_size = np.random.randint(samplesize_lo, samplesize_hi, size=n)
    allele_counts = np.vstack([sp.stats.multinomial(n=samples, p=af).rvs() for samples,af in zip(sample_size, popAF)])
    return alpha, popAF, allele_counts

results = []
for i in range(100):
    print(i)
    alpha, popAF, allele_counts = simulateAF(k=15, concentration=10, n=5)

    # Key data derived from simulated data
    # Number of alleles
    k = len(alpha)
    # Number of populations
    n = popAF.shape[0]
    # Sample size for each population
    sample_size = np.apply_along_axis(sum, 1, allele_counts)

    # Define and fit model
    with pm.Model() as m1:
        frac = pm.Dirichlet('frac', a=np.ones(k))
        conc = pm.Lognormal('conc', mu=1, sigma=1)
        y = pm.DirichletMultinomial('y', n=sample_size, a=frac*conc, shape=(n,k), observed=allele_counts)

    with m1:
        idata1 = pm.sample()

    # Store results
    result = {
        'lo': az.summary(idata1)[:k]['hdi_3%'],
        'hi':az.summary(idata1)[:k]['hdi_97%'],
        'alpha': alpha,
        'popAF': popAF
    }
    results.append(result)

# View summary results
def alpha_in_hdi(result):
    """What proportion of true global allele frequencies are in the estimated hdi?"""
    in_hdi = (result['lo'] < result['alpha']) & (result['alpha'] < result['hi'])
    aih = in_hdi.mean()
    return aih

def popAF_in_hdi(result):
    """What proportion of true population allele frequencies are in the estimated hdi?
    Note that the hdi applies to the global allele frequency."""
    in_hdi = np.apply_along_axis(lambda x: (result['lo']<x)&(x<result['hi']), 1, result['popAF'])
    pih = in_hdi.mean()
    return pih

aih = [alpha_in_hdi(result) for result in results]
pih = [popAF_in_hdi(result) for result in results]

# 95-97% simulated global AF are within the estimated hdi intervals
np.mean(aih)
# 50-80% simulated population AF are within the estimated hdi intervals
# This is expected to be <95% because the hdi is on the global AF
np.mean(pih)

fig, ax = plt.subplots(2)
ax[0].hist(aih)
ax[0].axvline(np.mean(aih), color="black")
ax[1].hist(pih)
ax[1].axvline(np.mean(pih), color="black")
plt.show()
