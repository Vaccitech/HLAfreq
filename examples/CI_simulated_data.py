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
import HLAfreq
from HLAfreq import HLAfreq_pymc as HLAhdi

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
    # Flat prior
    # alpha = sp.stats.dirichlet(np.ones(k)).rvs()
    
    # prior with one high freq allele
    alpha = np.ones(k)
    # alpha[0] = 5*k
    
    alpha = sp.stats.dirichlet(alpha).rvs()

    alpha = alpha.reshape(-1)
    # The Dirichlet distribution which populations are drawn from
    dd = sp.stats.dirichlet(alpha*concentration)
    # Population AF
    popAF =  dd.rvs(n)
    sample_size = np.random.randint(samplesize_lo, samplesize_hi, size=n)
    allele_counts = np.vstack([sp.stats.multinomial(n=samples, p=af).rvs() for samples,af in zip(sample_size, popAF)])
    return alpha, popAF, allele_counts

def serialise_aftab(allele_counts):
    aftab = pd.DataFrame(allele_counts)
    aftab['sample_size'] = aftab.apply(sum, axis=1)
    aftab['population'] = aftab.index
    aftab = pd.melt(aftab,
            id_vars=['population', 'sample_size'],
            value_vars=range(allele_counts.shape[1]),
            value_name='allele_count',
            var_name="allele")
    aftab['loci'] = "X"
    aftab['allele_freq'] = aftab.allele_count/aftab.sample_size
    aftab.allele = aftab.allele+1000
    aftab.population = aftab.population.astype(str)
    aftab.allele = aftab.allele.astype(str)
    return aftab


results = []
for k in range(10, 110, 10):
    for i in range(1):
        print(f"{i} for k={k}")
        alpha, popAF, allele_counts = simulateAF(k=k, concentration=100, n=5)

        # # Key data derived from simulated data
        # # Number of alleles
        # k = len(alpha)
        # # Number of populations
        # n = popAF.shape[0]
        # # Sample size for each population
        # sample_size = np.apply_along_axis(sum, 1, allele_counts)

        # # Define and fit model
        # with pm.Model() as m1:
        #     frac = pm.Dirichlet('frac', a=np.ones(k))
        #     conc = pm.Lognormal('conc', mu=1, sigma=1)
        #     y = pm.DirichletMultinomial('y', n=sample_size, a=frac*conc, shape=(n,k), observed=allele_counts)

        # with m1:
        #     idata1 = pm.sample()

        # Store results
        # result = {
        #     'lo': az.summary(idata1)[:k]['hdi_3%'],
        #     'hi':az.summary(idata1)[:k]['hdi_97%'],
        #     'alpha': alpha,
        #     'popAF': popAF
        # }

        aftab = serialise_aftab(allele_counts)
        caf = HLAfreq.combineAF(aftab)
        hdi = HLAhdi.AFhdi(aftab, prior=[1/k]*k)
        caf = pd.merge(caf, hdi, how="left", on="allele")

        # Store results
        result = {
            'lo': caf.lo,
            'hi': caf.hi,
            'post_mean': caf.post_mean,
            'allele_freq': caf.allele_freq,
            'alpha': alpha,
            'popAF': popAF,
            'k':k,
            'total_c': caf.c.sum()
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

#               #
#   k/c_total   #
#               #

plt.plot([r['k']/r['total_c'] for r in results]);plt.show()

def reduction_vs_kc(estimate, results, ax):
    # post mean reduction vs k/c
    for result in results:
        ax.scatter(
            [result['k']/result['total_c']] * result['k'],
            result['alpha'] - result[estimate]
        )
    ax.set(xlabel='k/total_c', ylabel=f'{estimate} reduction')

def est_vs_true(estimate, results, ax):
    # post_mean vs alpha coloured by k/c
    ax.scatter(
        pd.concat([pd.Series(r['alpha']) for r in results]),
        pd.concat([r[estimate] for r in results]),
        c=pd.concat([pd.Series([r['k']/result['total_c']]*r['k']) for r in results])
    )
    ax.plot([0,1], [0,1], c='black', linestyle="--")
    ax.set(xlabel='alpha', ylabel=estimate)

def kc_effect(results):
    fig, ax = plt.subplots(2,2)
    reduction_vs_kc('post_mean', results, ax[0,0])
    est_vs_true('post_mean', results, ax[1,0])
    reduction_vs_kc('allele_freq', results, ax[0,1])
    est_vs_true('allele_freq', results, ax[1,1])
    plt.show()

kc_effect(results)
kc_effect(results[:100])
kc_effect(results[100:110])
kc_effect(results[110:120])
kc_effect(results[120:130])

fig,ax = plt.subplots()
# post_mean vs alpha coloured by k/c
ax.scatter(
    pd.concat([pd.Series(r['alpha']) for r in results]),
    pd.concat([r['post_mean'] for r in results]),
    c=pd.concat([pd.Series([r['k']]*r['k']) for r in results])
)
ax.plot([0,1], [0,1], c='black', linestyle="--")
ax.set(xlabel='alpha', ylabel='post_mean')