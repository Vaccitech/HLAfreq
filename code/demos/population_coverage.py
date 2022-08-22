"""
How much of the population is covered by the top n alleles?
What is the population coverage of different supertypes?
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Which supertypes are represented
super1 = pd.read_csv("data/HLA1supertypes_Sidney2008.csv")
super1.allele = super1.allele.apply(lambda x: x[:4]+':'+x[4:])

AFtab = pd.read_csv("data/example/Uganda_raw.csv")
# Weighted average alleles
wav = scrapeAF.combineAF(AFtab)
# Add HLA1 supertypes
wav = pd.merge(wav, super1)
wava = wav[wav.loci == "A"].sort_values('allele_freq', ascending=False, ignore_index = True)

# Plot allele frequencies coloured by supertype
sns.barplot(
    y=wava.allele,
    x=wava.allele_freq,
    hue=wava.supertype,
    dodge=False
    )
plt.show()

plt.plot(wava.allele_freq.cumsum(), label="Cumulative AF")
plt.plot(wava.allele_freq.cumsum().apply(scrapeAF.population_coverage), label="Cumulative coverage")
plt.axhline(y=1, color="black", linestyle="--")
plt.legend()
plt.show()

for supertype in wava.supertype.unique():
    print(supertype)
    supertypedf = wava[wava.supertype == supertype].reset_index()
    plt.plot(supertypedf.allele_freq.cumsum().apply(scrapeAF.population_coverage), label=supertype)
plt.legend()
plt.show()

mask = wava.allele_freq.cumsum() > 0.95
wava[~mask]
wava[~mask].supertype.unique()

