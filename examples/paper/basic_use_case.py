import matplotlib.pyplot as plt
import HLAfreq
from HLAfreq import HLAfreq_pymc

base_url = HLAfreq.makeURL(country="Mongolia", standard="g", locus="DQB1")
aftab = HLAfreq.getAFdata(base_url)
caf = HLAfreq.combineAF(aftab)
hdi = HLAfreq_pymc.AFhdi(aftab, credible_interval=0.95)
HLAfreq.plotAF(caf, aftab, hdi=hdi, compound_mean=hdi)

caf = caf.sort_values('allele_freq', ascending=False, ignore_index=True)
plt.scatter(caf.allele, caf.allele_freq.cumsum().apply(HLAfreq.population_coverage))
plt.grid(True); plt.xticks(rotation=90); plt.tight_layout()
plt.xlabel('Allele'); plt.ylabel('Cumulative population coverage')
plt.tight_layout(); plt.show()
