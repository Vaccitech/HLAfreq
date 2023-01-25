import math
import pymc as pm
import numpy as np
import arviz as az
import pandas as pd
import HLAfreq

def _make_c_array(AFtab, weights="2n", datasetID="population"):
    df = AFtab.copy()
    df = HLAfreq.unmeasured_alleles(df, datasetID)
    try:
        df['2n'] = df.sample_size * 2
    except:
        print("column '2n' could not be created")
    df['c'] =  df.allele_freq * df[weights]

    # Sort by alleles so it matches the combined alleles
    df = df.sort_values('allele')
    c_array = np.array(df.groupby(datasetID).c.apply(list).tolist())
    allele_names = sorted(df.allele.unique())
    # Imperfect check that allele order matches between caf and c_array.
    # caf is sorted automatically so should match sorted AFloc
    # Therefore we check that sorted AFloc matches c_array
    # The check is that the sum of allele i is the same
    for a,b in zip(
        np.apply_along_axis(sum, 0, c_array),
        df.groupby('allele').c.sum()
        ):
        assert math.isclose(a,b), "Error making c_array sum of single allele frequency differs between c_array and AFloc"
    return c_array, allele_names

def _fit_Dirichlet_Multinomial(c_array, prior=[]):
    # Number of populations
    n = c_array.shape[0]
    # number of alleles
    k = c_array.shape[1]

    # Round c array so that it and effective ssamples are whole numbers
    # for the multinomial
    c_array = np.round(c_array)
    effective_samples = np.apply_along_axis(sum, 1, c_array)
    if len(prior)==0:
        prior = np.ones(k)
    assert len(prior) == k, "For k alleles, prior must be length k"
    with pm.Model() as mod:
        frac = pm.Dirichlet('frac', a=prior)
        conc = pm.Lognormal('conc', mu=1, sigma=1)
        y = pm.DirichletMultinomial('y', n=effective_samples, a=frac*conc, shape=(n,k), observed=c_array)

    with mod:
        idata = pm.sample()
    return idata

def AFhdi(AFtab, weights="2n", datasetID="population", credible_interval=0.95, prior=[]):
    """Calculate mean and high posterior density interval on combined allele frequency.
    Fits a Marginalized Dirichlet-Multinomial Model in PyMc as described [here](https://docs.pymc.io/en/v3/pymc-examples/examples/mixture_models/dirichlet_mixture_of_multinomials.html).
    
    In brief, the global allele frequency is modelled as a Dirichlet distribution,
    and each population (defined by `datasetID`) is a Dirichlet distribution draw from
    the global Dirichlet distribution, and the observed allele count data of that
    population is multinomial count data drawn from the population Dirichlet distribution.

    The observed allele frequencies are transformed into allele counts using `weights`.
    The variability of population allele frequencies around the global mean is defined
    by a latent, lognormal variable `conc`.

    Args:
        AFtab (pd.DataFrame): Table of allele frequency data
        weights (str, optional): Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.
        datasetID (str, optional): Unique identifier column for study. Defaults to 'population'.
        credible_interval (float, optional): The size of the credible interval requested. Defaults to 0.95.
        prior (list, optional): Prior vector for global allele frequency. Order should match alphabetical alleles, i.e. the first value is used for the alphabetically first allele.

    Returns:
        np.array: Pairs of high density interval limits, allele name, and posterior mean.
            as a 4 by n array.
            In alphabetical order of alleles, regardless of input order.
            This way it matches the output of combineAF().
    """
    c_array, allele_names = _make_c_array(AFtab, weights, datasetID)
    idata = _fit_Dirichlet_Multinomial(c_array, prior)
    hdi = az.hdi(idata, hdi_prob=credible_interval).frac.values
    post_mean = az.summary(idata, var_names='frac')['mean']
    post = pd.DataFrame([hdi[:,0], hdi[:,1], allele_names, post_mean]).T
    post.columns = ['lo','hi','allele','post_mean']
    return post