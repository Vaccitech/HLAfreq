import matplotlib.pyplot as plt
from HLAfreq import makeURL, getAFdata, combineAF, plotAF, population_coverage

base_url = makeURL("Mongolia", standard="g", locus="DQB1")
aftab = getAFdata(base_url)
aftab.population.unique()
caf = combineAF(aftab)
plotAF(caf, aftab, credible_interval=0.95)

caf = caf.sort_values('allele_freq', ascending=False, ignore_index=True)
plt.scatter(caf.allele, caf.allele_freq.cumsum().apply(population_coverage))
plt.grid(True); plt.xticks(rotation=90); plt.tight_layout()
plt.xlabel('Allele'); plt.ylabel('Cumulative population coverage')
plt.tight_layout(); plt.show()
