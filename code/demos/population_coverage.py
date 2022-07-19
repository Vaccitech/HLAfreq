"""
How much of the population is covered by the top n alleles?
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import matplotlib.pyplot as plt

alleles = pd.read_csv("data/example/Uganda_raw.csv")
alleles = scrapeAF.formatAF(alleles)

# Weighted average alleles
wav = scrapeAF.combineAF(alleles)

wav = wav.sort_values('allele_freq', ascending=False, ignore_index = True)

plt.plot(wav.allele_freq)
plt.show()

plt.plot(wav.allele_freq.cumsum())
plt.show()

plt.plot(wav.allele_freq.cumsum().apply(scrapeAF.population_coverage))
plt.show()