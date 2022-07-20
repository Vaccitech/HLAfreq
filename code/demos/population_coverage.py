"""
How much of the population is covered by the top n alleles?
"""

import code.scrapeAF as scrapeAF
import pandas as pd
import matplotlib.pyplot as plt

AFtab = pd.read_csv("data/example/Uganda_raw.csv")
AFtab = scrapeAF.formatAF(AFtab)

loci = AFtab.loci.unique()
for locus in loci:
    populations = AFtab[AFtab.loci == locus].population.unique()
    for population in populations:
# locus = "A"
# population = 'Uganda Kampala'
        popAF = AFtab[(AFtab.population == population) & (AFtab.loci == locus)]
        pop_sample_size = popAF.sample_size.unique()
        assert len(pop_sample_size) == 1, "pop_sample_size must be 1, not %s" %len(pop_sample_size)
        pop_sample_size = pop_sample_size[0]
        ualleles = AFtab[AFtab.loci == locus].allele.unique()
        missing_alleles = [allele for allele in ualleles if not allele in popAF.allele.values]
        missing_rows = [(al, locus, population, 0, 0, pop_sample_size) for al in missing_alleles]
        missing_rows = pd.DataFrame(missing_rows, columns=['allele','loci','population','allele_freq','carriers%','sample_size'])
        AFtab = pd.concat([AFtab, missing_rows], ignore_index=True)

# Weighted average alleles
wav = scrapeAF.combineAF(AFtab)

wav = wav.sort_values('allele_freq', ascending=False, ignore_index = True)

wava = wav[wav.loci == "A"]

plt.bar(wava.allele, wava.allele_freq)
plt.xticks(rotation=90)
plt.show()

# cumsum exceeds 1, possibly because rare alleles are not found
# in all populations within a country

plt.plot(wava.allele_freq.cumsum())
plt.plot(wava.allele_freq.cumsum().apply(scrapeAF.population_coverage))
plt.axhline(y=1, color="black", linestyle="--")
plt.show()