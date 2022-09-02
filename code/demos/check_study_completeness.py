"""
Check that study-loci allele frequency sums to ~1.
Some studies report alleles that do not sum to 1.
Also download errors can lead to incomplete datasets.
Some download settings can also lead to unexpected results.
E.g. some studies report alleles to mixed resolution, therefore;
if you search for a specific allele resolution you will download
an incomplete dataset.

Having the complete dataset is important, especially when
calculating the Dirichlet distribution as reported alleles
will be over weighted if some alleles are unreported.
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

country = 'Thailand'
# Because resolution equals 2, an incomplete dataset is returned
base_url = scrapeAF.makeURL(country, standard="s")
aftab = scrapeAF.getAFdata(base_url)
aftab.to_csv("data/example/incomplete_raw.csv", index=False)

df = pd.read_csv("data/example/incomplete_raw.csv")

scrapeAF.study_completeness(df)

# How to filter incomplete studies...