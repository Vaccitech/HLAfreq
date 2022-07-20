"""
How much of the population is covered by the top n alleles?
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import matplotlib.pyplot as plt

AFtab = pd.read_csv("data/example/Uganda_raw.csv")
AFtab = scrapeAF.formatAF(AFtab)

# Add unreported alleles so AF sums to 1
AFtab = scrapeAF.unmeasured_alleles(AFtab)

# Weighted average alleles
wav = scrapeAF.combineAF(AFtab)

wava = wav[wav.loci == "A"].sort_values('allele_freq', ascending=False, ignore_index = True)

plt.bar(wava.allele, wava.allele_freq)
plt.xticks(rotation=90)
plt.show()

plt.plot(wava.allele_freq.cumsum(), label="Cumulative AF")
plt.plot(wava.allele_freq.cumsum().apply(scrapeAF.population_coverage), label="Cumulative coverage")
plt.axhline(y=1, color="black", linestyle="--")
plt.legend()
plt.show()
