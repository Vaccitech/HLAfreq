import matplotlib.pyplot as plt
from HLAfreq import makeURL, getAFdata, combineAF, plotAF, population_coverage
from HLAfreq import HLAfreq_pymc

base_url = makeURL("Mongolia", standard="g", locus="DQB1")
aftab = getAFdata(base_url)
caf = combineAF(aftab)
hdi = HLAfreq_pymc.AFhdi(aftab, credible_interval=0.95)
plotAF(caf, aftab, hdi=hdi, compound_mean=hdi)

caf = caf.sort_values('allele_freq', ascending=False, ignore_index=True)
plt.scatter(caf.allele, caf.allele_freq.cumsum().apply(population_coverage))
plt.grid(True); plt.xticks(rotation=90); plt.tight_layout()
plt.xlabel('Allele'); plt.ylabel('Cumulative population coverage')
plt.tight_layout(); plt.show()
