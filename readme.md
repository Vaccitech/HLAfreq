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