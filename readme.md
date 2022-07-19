# scrapeAF

The aim of this tool is to be able to estimate the HLA allele
frequencies at large multi population scale, e.g. globally.
Scrape allele frequency data from
[allele frequencies.net](http://www.allelefrequencies.net/).
This gives frequencies by population, which can be averaged
across studies, weighted by sample size. Frequencies can then
be averaged across populations if a relative weighting can
found, i.e. population size.

```
import sys
sys.path.append("/home/dwells/work/tools/scrapeAF/code")
import scrapeAF
```

When selecting alleles to represent a population, all supertypes
should be represented.

When combining extimates of allele frequency larger samples should
weigh more heavily. However, when combining estiamtes that disagree
the uncertainty on the estimate should be larger than when combining
equivalent samples that agree.

When combining studies we must weight averages by sample size but
when combining countries are we not better to weight by national population?
Should sample size also have an effect?