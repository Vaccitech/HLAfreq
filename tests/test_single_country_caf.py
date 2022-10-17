"""
Tests on a toy allele frequency data set

Tests dropping incomplete studies, merging different resolution
alleles, and combing allele frequency estimates.
"""
import pytest
import HLAfreq

import pandas as pd
from scipy.stats import dirichlet

dfa =  pd.DataFrame({
    'allele':['A*02:03','A*02:05','A*02:07'],
    'loci':['A','A','A'],
    'population':['test1','test1','test1'],
    'allele_freq':[0.1, 0.3, 0.6],
    'sample_size':[10,10,10]
})
dfb =  pd.DataFrame({
    'allele':['A*02:03','A*02:05:01','A*02:05:02'],
    'loci':['A','A', 'A'],
    'population':['test2','test2', 'test2'],
    'allele_freq':[0.5, 0.25, 0.25],
    'sample_size':[5, 5, 5]
})
dfc =  pd.DataFrame({
    'allele':['A*02:03','A*02:05'],
    'loci':['A','A'],
    'population':['test3','test3'],
    'allele_freq':[0.5, 0.1],
    'sample_size':[5, 5]
})
aftab = pd.concat([dfa, dfb, dfc])

# Drop incomplete test3
aftab = HLAfreq.only_complete(aftab)
def test_drop_incomplete_test():
    assert not any(aftab.population == 'test3')

# Collapse 3 field alleles
aftab = HLAfreq.decrease_resolution(aftab, 2)
def test_resolution_decreased_correctly():
    mask = (aftab.allele=="A*02:05") & (aftab.population=='test2')
    assert sum(mask) == 1
    collapsed = aftab[mask]
    assert all(collapsed.allele_freq == 0.5)
    assert all(collapsed.sample_size == 5)

# Combine allele frequencies
caf = HLAfreq.combineAF(aftab)
def test_combined_allele_freq():
    assert all(caf.allele_freq == dirichlet([8,12,13]).mean())