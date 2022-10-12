"""
Working with priors to improve estimates of undersampled populations.
Priors can be based on other related, populations, or specified by
hand.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Venezuela is one of the least sampled countries on allelefrequency.net
We found only 55 gold standard samples for HLA A.
"""

venezuelaAF = pd.read_csv("data/example/globalPCA/Venezuela_raw.csv")
venezuelaAF = scrapeAF.only_complete(venezuelaAF)
venezuelaAF = scrapeAF.decrease_resolution(venezuelaAF, 2)
cafV = scrapeAF.combineAF(venezuelaAF)
cafV['study'] =  'Venezuela'

scrapeAF.plotAFprob(cafV, venezuelaAF)

"""
We can specify priors manually. For simplicity we use
a prior of 1 for each observed allele but this could be any
n length list of positive real numbers.

The interpretation of the prior is that n copies were found for each
allele. Note that `combineAF()` uses the first value in the prior for
the first allele alphabetically and so on.

By viewing the posterior we can see that the prior has moved the
estimate away from the observed study.
"""
manual_prior = [10, 1, 1, 1, 1]

# View prior distribution
scrapeAF.plotAFprob(concentration=manual_prior)

# Combine Allele Freq of study with Manual Prior
cafMP = scrapeAF.combineAF(venezuelaAF, alpha=manual_prior)

# View posterior distribution
scrapeAF.plotAFprob(cafMP, AFtab=venezuelaAF)

"""
Alternatively we can use another population as our prior.
For example we will use data for Colombia, a neighbouring country
with many more samples for HLA A.

When using other studies as a prior we can downweight their samples
so they count less than the samples in our actual population.
The amount of downweighting depends on how relevent you believe the
other study is. Here we use downweight by a factor 0.01 so that the
Venezuela data is not drowned out by the Colombia data.

When recalculating the prior to plot, sample size is doubled
as each person in the study is diploid and so provides two alleles.
"""

colombiaAF = pd.read_csv("data/example/globalPCA/Colombia_raw.csv")
colombiaAF = scrapeAF.only_complete(colombiaAF)
colombiaAF = scrapeAF.decrease_resolution(colombiaAF, 2)
cafC = scrapeAF.combineAF(colombiaAF)
cafC['population'] =  'Colombia'

study_prior = cafC.copy()
study_prior.sample_size = study_prior.sample_size * 0.01

# View prior
scrapeAF.plotAFprob(
    concentration = (2 * study_prior.sample_size * study_prior.allele_freq).tolist()
)

# (2 * study_prior.sample_size * study_prior.allele_freq)-(study_prior.c * 0.01)
# plt.scatter(
# (2 * study_prior.sample_size * study_prior.allele_freq).tolist(),
# (study_prior.c * 0.01).tolist()
# )
# plt.show()

# Combine Allele Frequency with Study Prior
cafSP = scrapeAF.combineAF(
    pd.concat([venezuelaAF, study_prior], join="inner")
    )

# View posterior
scrapeAF.plotAFprob(cafSP, pd.concat([venezuelaAF, study_prior], join="inner"))
scrapeAF.plotAFprob(cafSP, pd.concat([venezuelaAF, study_prior], join="inner"), alleles=venezuelaAF.allele.tolist())
