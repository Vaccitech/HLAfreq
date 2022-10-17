# HLAfreq

`HLAfreq` allows you to estimate the HLA allele
frequencies at large multi population scale, e.g. combine data from
several studies within a country or combine countries.

Automated download of allele frequency data download from 
[allele frequencies.net](http://www.allelefrequencies.net/).

## Details
Estimates are combined by modelling allele frequency as a 
Dirichlet distribution which defines the probability of drawing each
allele. When combining studies their estimates are weighted as 2x sample size by
default. Sample size is doubled as each person in the study
contributes two alleles. Alternative weightings can be used
for example population size when averaging across countries.

A note on credible intervals...

When selecting a panel of HLA alleles to represent a population,
allele frequency is not the only thing to consider. Depending on
the purpose of the panel, you should include a range of loci and
supertypes (groups alleles sharing binding specificies).

## Install

## Minimal example

## Detailed examples

## Citation
