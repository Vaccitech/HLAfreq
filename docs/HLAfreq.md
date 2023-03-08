---
description: |
    API documentation for modules: HLAfreq, HLAfreq.HLAfreq, HLAfreq.HLAfreq_data, HLAfreq.HLAfreq_pymc.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...


    
# Module `HLAfreq` {#HLAfreq}




    
## Sub-modules

* [HLAfreq.HLAfreq](#HLAfreq.HLAfreq)
* [HLAfreq.HLAfreq_data](#HLAfreq.HLAfreq_data)
* [HLAfreq.HLAfreq_pymc](#HLAfreq.HLAfreq_pymc)






    
# Module `HLAfreq.HLAfreq` {#HLAfreq.HLAfreq}

HLAfreq

Download allele frequency data from
[allelefrequencies.net](www.allelefrequencies.net). Allele
frequencies from different populations can be combined to
estimate HLA frequencies of countries or other regions such as
global HLA frequencies.




    
## Functions


    
### Function `Npages` {#HLAfreq.HLAfreq.Npages}



    
> `def Npages(bs)`


How many pages of results are there?


###### Args

**```bs```** :&ensp;<code>bs4.BeautifulSoup</code>
:   BS object of allelefrequencies.net results page



###### Returns

<code>int</code>
:   Total number of results pages



    
### Function `alleles_unique_in_study` {#HLAfreq.HLAfreq.alleles_unique_in_study}



    
> `def alleles_unique_in_study(AFtab, datasetID='population')`


Are all alleles unique in each study? Checks that no alleles are reported more than once in a single study. Study is defined by <code>datasetID</code>.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Allele frequency data


**```datasetID```** :&ensp;<code>str</code>, optional
:   Unique identifier column to define study. Defaults to 'population'.



###### Returns

<code>bool</code>
:   <code>True</code> on if no alleles occur more than once in any study, otherwise <code>False</code>.



    
### Function `betaAB` {#HLAfreq.HLAfreq.betaAB}



    
> `def betaAB(alpha)`


Given the alpha vector defining a Dirichlet distribution calculate the a b values for all composite beta distributions.


###### Args

**```alpha```** :&ensp;<code>list</code>
:   Values defining a Dirichlet distribution. This will be the prior (for a naive distribution) or the prior + caf.c for a posterior distribution.



###### Returns

<code>list</code>
:   List of a b values defining beta values, i.e. for each allele it is the number of times it was and wasn't observed.



    
### Function `check_resolution` {#HLAfreq.HLAfreq.check_resolution}



    
> `def check_resolution(AFtab)`


Check if all alleles in AFtab have the same resolution.
Will print the number of records with each resolution.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Allele frequency data



###### Returns

<code>bool</code>
:   True only if all alleles have the same resolution, else False.



    
### Function `collapse_reduced_alleles` {#HLAfreq.HLAfreq.collapse_reduced_alleles}



    
> `def collapse_reduced_alleles(AFtab, datasetID='population')`




    
### Function `combineAF` {#HLAfreq.HLAfreq.combineAF}



    
> `def combineAF(AFtab, weights='2n', alpha=[], datasetID='population', format=True, ignoreG=True, add_unmeasured=True, complete=True, resolution=True, unique=True)`


Combine allele frequencies from multiple studies. <code>datasetID</code> is the unique identifier for studies to combine.
Allele frequencies combined using a Dirichlet distribution where each study's contribution to the concentration parameter is $2 * sample_size * allele_frequency$.
Sample size is doubled to get <code>2n</code> due to diploidy. If an alternative <code>weights</code> is set it is not doubled.
The total concentration parameter of the Dirichlet distribution is the contributions from all studies plus the prior <code>alpha</code>.
If <code>alpha</code> is not set the prior defaults to 1 observation of each allele.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Table of Allele frequency data


**```weights```** :&ensp;<code>str</code>, optional
:   Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.


**```alpha```** :&ensp;<code>list</code>, optional
:   Prior to use for Dirichlet distribution. Defaults to [].


**```datasetID```** :&ensp;<code>str</code>, optional
:   Unique identifier column for study. Defaults to 'population'.


**```format```** :&ensp;<code>bool</code>, optional
:   Run <code>[formatAF()](#HLAfreq.HLAfreq.formatAF "HLAfreq.HLAfreq.formatAF")</code>. Defaults to True.


**```ignoreG```** :&ensp;<code>bool</code>, optional
:   Treat allele G groups as normal, see <code>[formatAF()](#HLAfreq.HLAfreq.formatAF "HLAfreq.HLAfreq.formatAF")</code>. Defaults to True.


**```add_unmeasured```** :&ensp;<code>bool</code>, optional
:   Add unmeasured alleles to each study. This is important to ensure combined allele frequencies sum to 1. See <code>add\_unmeasured()</code>. Defaults to True.


**```complete```** :&ensp;<code>bool</code>, optional
:   Check study completeness. Uses default values for <code>[incomplete\_studies()](#HLAfreq.HLAfreq.incomplete\_studies "HLAfreq.HLAfreq.incomplete\_studies")</code>. If you are happy with your study completeness can be switched off with False. Defaults to True.


**```resolution```** :&ensp;<code>bool</code>, optional
:   Check that all alleles have the same resolution, see <code>[check\_resolution()](#HLAfreq.HLAfreq.check\_resolution "HLAfreq.HLAfreq.check\_resolution")</code>. Defaults to True.


**```unique```** :&ensp;<code>bool</code>, optional
:   Check that each allele appears no more than once per study. See <code>[alleles\_unique\_in\_study()](#HLAfreq.HLAfreq.alleles\_unique\_in\_study "HLAfreq.HLAfreq.alleles\_unique\_in\_study")</code>. Defaults to True.



###### Returns

<code>pd.DataFrame</code>
:   Allele frequencies after combining estimates from all studies.


*allele_freq* is the combined frequency estimate from the Dirichlet mean where the concentration is <code>alpha</code> + <code>c</code>.
*alpha* is the prior used for the Dirichlet distribution.
*c* is the observations used for the Dirichlet distribution.
*sample_size* is the total sample size of all combined studies.
*wav* is the weighted average.

    
### Function `decrease_resolution` {#HLAfreq.HLAfreq.decrease_resolution}



    
> `def decrease_resolution(AFtab, newres, datasetID='population')`


Decrease allele resolution to a specified value so all alleles have the same resolution.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Allele frequency data.


**```newres```** :&ensp;<code>int</code>
:   The desired number of fields for resolution.


**```datasetID```** :&ensp;<code>str</code>, optional
:   Column to use as stud identifier. Defaults to 'population'.



###### Returns

<code>pd.DataFrame</code>
:   Allele frequency data with all alleles of requested resolution.



    
### Function `default_prior` {#HLAfreq.HLAfreq.default_prior}



    
> `def default_prior(k)`


Calculate a default prior, 1 observation of each class.


###### Args

**```k```** :&ensp;<code>int</code>
:   Number of classes in the Dirichlet distribution.



###### Returns

<code>list</code>
:   List of k 1s to use as prior.



    
### Function `duplicated_sample_size` {#HLAfreq.HLAfreq.duplicated_sample_size}



    
> `def duplicated_sample_size(AFtab)`


Returns True if any loci has more than 1 unique sample size

    
### Function `formatAF` {#HLAfreq.HLAfreq.formatAF}



    
> `def formatAF(AFtab, ignoreG=True)`


Format allele frequency table.
Convert sample_size and allele_freq to numeric data type.
Removes commas from sample size. Removes "(*)" from allele frequency if 
<code>ignoreG</code> is <code>True</code>. <code>[formatAF()](#HLAfreq.HLAfreq.formatAF "HLAfreq.HLAfreq.formatAF")</code> is used internally by combineAF and getAFdata by default.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Allele frequency data downloaded from allelefrequency.net using <code>[getAFdata()](#HLAfreq.HLAfreq.getAFdata "HLAfreq.HLAfreq.getAFdata")</code>.


**```ignoreG```** :&ensp;<code>bool</code>, optional
:   Treat G group alleles as normal. See <http://hla.alleles.org/alleles/g_groups.html> for details. Defaults to True.



###### Returns

<code>pd.DataFrame</code>
:   The formatted allele frequency data.



    
### Function `getAFdata` {#HLAfreq.HLAfreq.getAFdata}



    
> `def getAFdata(base_url, format=True, ignoreG=True)`


Get all allele frequency data from a search base url. Iterates over all
    pages regardless of which page is based.


###### Args

**```base_url```** :&ensp;<code>str</code>
:   URL for base search.


**```format```** :&ensp;<code>bool</code>
:   Format the downloaded data using <code>[formatAF()](#HLAfreq.HLAfreq.formatAF "HLAfreq.HLAfreq.formatAF")</code>.


**```ignoreG```** :&ensp;<code>bool</code>
:   treat allele G groups as normal. See <http://hla.alleles.org/alleles/g_groups.html> for details. Default = True



###### Returns

<code>pd.DataFrame</code>
:   allele frequency data parsed into a pandas dataframe



    
### Function `id_duplicated_allele` {#HLAfreq.HLAfreq.id_duplicated_allele}



    
> `def id_duplicated_allele(grouped)`


Reports the allele that has mupltiple sample sizes

    
### Function `incomplete_studies` {#HLAfreq.HLAfreq.incomplete_studies}



    
> `def incomplete_studies(AFtab, llimit=0.95, ulimit=1.1, datasetID='population')`


Report any studies with allele freqs that don't sum to 1


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Dataframe containing multiple studies


**```llimit```** :&ensp;<code>float</code>, optional
:   Lower allele_freq sum limit that counts as complete. Defaults to 0.95.


**```ulimit```** :&ensp;<code>float</code>, optional
:   Upper allele_freq sum limit that will not be reported. Defaults to 1.1.


**```datasetID```** :&ensp;<code>str</code>
:   Unique identifier column for study



    
### Function `makeURL` {#HLAfreq.HLAfreq.makeURL}



    
> `def makeURL(country='', standard='s', locus='', resolution_pattern='bigger_equal_than', resolution=2, region='', ethnic='', study_type='', dataset_source='', sample_year='', sample_year_pattern='', sample_size='', sample_size_pattern='')`


Create URL for search of allele frequency net database. All arguments are documented [here](<http://www.allelefrequencies.net/extaccess.asp>)

###### Args

**```country```** :&ensp;<code>str</code>, optional
:   Country name to retrieve records from. Defaults to "".


**```standard```** :&ensp;<code>str</code>, optional
:   Filter study quality standard to this or higher. {'g', 's', 'a'} Gold, silver, all. Defaults to 's'.


**```locus```** :&ensp;<code>str</code>, optional
:   The locus to return allele data for. Defaults to "".


**```resolution_pattern```** :&ensp;<code>str</code>, optional
:   Resolution comparitor {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Filter created using <code>resolution</code> and <code>resolution\_pattern</code>. Defaults to "bigger_equal_than".


**```resolution```** :&ensp;<code>int</code>, optional
:   Number of fields of resolution of allele. Filter created using <code>resolution</code> and <code>resolution\_pattern</code>. Defaults to 2.


**```region```** :&ensp;<code>str</code>, optional
:   Filter to geographic region. {Asia, Australia, Eastern Europe, ...}. All regions listed [here](<http://www.allelefrequencies.net/pop6003a.asp>). Defaults to "".


**```ethnic```** :&ensp;<code>str</code>, optional
:   Filter to ethnicity. {"Amerindian", "Black", "Caucasian", ...}. All ethnicities listed [here](<http://www.allelefrequencies.net/pop6003a.asp>). Defaults to "".


**```study_type```** :&ensp;<code>str</code>, optional
:   Type of study. {"Anthropology", "Blood+Donor", "Bone+Marrow+Registry", "Controls+for+Disease+Study", "Disease+Study+Patients", "Other", "Solid+Organd+Unrelated+Donors", "Stem+cell+donors"}. Defaults to "".


**```dataset_source```** :&ensp;<code>str</code>, optional
:   Source of data. {"Literature", "Proceedings+of+IHWs", "Unpublished"}. Defaults to "".


**```sample_year```** :&ensp;<code>int</code>, optional
:   Sample year to compare to. Filter created using sample_year and sample_year_pattern. Defaults to "".


**```sample_year_pattern```** :&ensp;<code>str</code>, optional
:   Pattern to compare sample year to. Filter created using sample_year and sample_year_pattern. {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Defaults to "".


**```sample_size```** :&ensp;<code>int</code>, optional
:   Sample size to compare to. Filter created using sample_size and sample_size_pattern. Defaults to "".


**```sample_size_pattern```** :&ensp;<code>str</code>, optional
:   Pattern to compare sample size to. Filter created using sample_size and sample_size_pattern. {'equal', 'different', 'less_than', 'bigger_than', 'less_equal_than', 'bigger_equal_than'}. Defaults to "".



###### Returns

<code>str</code>
:   URL to search allelefrequencies.net



    
### Function `only_complete` {#HLAfreq.HLAfreq.only_complete}



    
> `def only_complete(AFtab, llimit=0.95, ulimit=1.1, datasetID='population')`


Returns only complete studies. Studies are only dropped if their population and loci are in noncomplete together.
This prevents throwing away data if another loci in the population is incomplete


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Dataframe containing multiple studies


**```llimit```** :&ensp;<code>float</code>, optional
:   Lower allele_freq sum limit that counts as complete. Defaults to 0.95.


**```ulimit```** :&ensp;<code>float</code>, optional
:   Upper allele_freq sum limit that will not be reported. Defaults to 1.1.


**```datasetID```** :&ensp;<code>str</code>
:   Unique identifier column for study. Defaults to 'population'.



###### Returns

<code>pd.DataFrame</code>
:   Allele frequency data of multiple studies, but only complete studies.



    
### Function `parseAF` {#HLAfreq.HLAfreq.parseAF}



    
> `def parseAF(bs)`


Generate a dataframe from a given html page


###### Args

**```bs```** :&ensp;<code>bs4.BeautifulSoup</code>
:   BeautifulSoup object from allelefrequencies.net page



###### Returns

<code>pd.DataFrame</code>
:   Table of allele, allele frequency, samplesize, and population



    
### Function `plotAF` {#HLAfreq.HLAfreq.plotAF}



    
> `def plotAF(caf=Empty DataFrame
Columns: []
Index: [], AFtab=Empty DataFrame
Columns: []
Index: [], cols=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'], datasetID='population', hdi=Empty DataFrame
Columns: []
Index: [], compound_mean=Empty DataFrame
Columns: []
Index: [])`


Plot combined allele frequencies, individual allele frequencies,
and credible intervals on combined allele frequency estimates.
Credible interval is only plotted if a value is given for <code>hdi</code>.
The plotted Credible interval is whatever was passed to HLAfreq_pymc.AFhdi() when calculating hdi.


###### Args

**```caf```** :&ensp;<code>pd.DataFrame</code>, optional
:   Combined allele frequency estimates from HLAfreq.combineAF. Defaults to pd.DataFrame().


**```AFtab```** :&ensp;<code>pd.DataFrame</code>, optional
:   Table of allele frequency data. Defaults to pd.DataFrame().


**```cols```** :&ensp;<code>list</code>, optional
:   List of colours to use for each individual dataset. Defaults to list(mcolors.TABLEAU_COLORS.keys()).


**```datasetID```** :&ensp;<code>str</code>, optional
:   Column used to define separate datasets. Defaults to "population".


**```weights```** :&ensp;<code>str</code>, optional
:   Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.


**```hdi```** :&ensp;<code>pd.DataFrame</code>, optional
:   The high density interval object to plot credible intervals. Produced by HLAfreq.HLA_pymc.AFhdi(). Defaults to pd.DataFrame().


**```compound_mean```** :&ensp;<code>pd.DataFrame</code>, optional
:   The high density interval object to plot post_mean. Produced by HLAfreq.HLA_pymc.AFhdi(). Defaults to pd.DataFrame().



    
### Function `plot_prior` {#HLAfreq.HLAfreq.plot_prior}



    
> `def plot_prior(concentration, ncol=2, psteps=1000, labels='')`


Plot probability density function for prior values.


###### Args

**```concentration```** :&ensp;<code>list</code>
:   Vector of the prior Dirichlet concentration values.


**```ncol```** :&ensp;<code>int</code>, optional
:   Number of columns. Defaults to 2.


**```labels```** :&ensp;<code>list</code>, optional
:   Labels for elements of concentration in the same order. Defaults to "".



    
### Function `population_coverage` {#HLAfreq.HLAfreq.population_coverage}



    
> `def population_coverage(p)`


Calculate the proportion of people with at least 1 copy of this allele
    assuming HWE.


###### Args

**```p```** :&ensp;<code>float</code>
:   Allele frequency



###### Returns

<code>float</code>
:   Sum of homozygotes and heterozygotes for this allele



    
### Function `single_loci` {#HLAfreq.HLAfreq.single_loci}



    
> `def single_loci(AFtab)`


Check that allele frequency data is only of one locus


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Allele frequency data



    
### Function `unmeasured_alleles` {#HLAfreq.HLAfreq.unmeasured_alleles}



    
> `def unmeasured_alleles(AFtab, datasetID='population')`


When combining AF estimates, unreported alleles can inflate frequencies
    so AF sums to >1. Therefore we add unreported alleles with frequency zero.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Formatted allele frequency data


**```datasetID```** :&ensp;<code>str</code>
:   Unique identifier column for study



###### Returns

<code>pd.DataFrame</code>
:   Allele frequency data with all locus alleles reported 
    for each dataset






    
# Module `HLAfreq.HLAfreq_data` {#HLAfreq.HLAfreq_data}

Data loaders




    
## Functions


    
### Function `load_HLA1supertypes_Sidney2008` {#HLAfreq.HLAfreq_data.load_HLA1supertypes_Sidney2008}



    
> `def load_HLA1supertypes_Sidney2008()`




    
### Function `load_countries` {#HLAfreq.HLAfreq_data.load_countries}



    
> `def load_countries()`







    
# Module `HLAfreq.HLAfreq_pymc` {#HLAfreq.HLAfreq_pymc}






    
## Functions


    
### Function `AFhdi` {#HLAfreq.HLAfreq_pymc.AFhdi}



    
> `def AFhdi(AFtab, weights='2n', datasetID='population', credible_interval=0.95, prior=[], conc_mu=1, conc_sigma=1, compare_models=True)`


Calculate mean and high posterior density interval on combined allele frequency.
Fits a Marginalized Dirichlet-Multinomial Model in PyMc as described [here](<https://docs.pymc.io/en/v3/pymc-examples/examples/mixture_models/dirichlet_mixture_of_multinomials.html>).

In brief, the global allele frequency is modelled as a Dirichlet distribution,
and each population (defined by <code>datasetID</code>) is a Dirichlet distribution draw from
the global Dirichlet distribution, and the observed allele count data of that
population is multinomial count data drawn from the population Dirichlet distribution.

The observed allele frequencies are transformed into allele counts using <code>weights</code>.
The variability of population allele frequencies around the global mean is defined
by a latent, lognormal variable <code>conc</code>.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Table of allele frequency data


**```weights```** :&ensp;<code>str</code>, optional
:   Column to be weighted by allele frequency to generate concentration parameter of Dirichlet distribution. Defaults to '2n'.


**```datasetID```** :&ensp;<code>str</code>, optional
:   Unique identifier column for study. Defaults to 'population'.


**```credible_interval```** :&ensp;<code>float</code>, optional
:   The size of the credible interval requested. Defaults to 0.95.


**```prior```** :&ensp;<code>list</code>, optional
:   Prior vector for global allele frequency. Order should match alphabetical alleles, i.e. the first value is used for the alphabetically first allele.


**```conc_mu```** :&ensp;<code>float</code>, optional
:   Mean to parameterise lognormal distribution of <code>conc</code> prior. Defaults to 1.


**```conc_sigma```** :&ensp;<code>float</code>, optional
:   Standard deviation to parameterise lognormal distribution of <code>conc</code> prior. Defaults to 1.


**```compare_models```** :&ensp;<code>bool</code>, optional
:   Check that default estimated allele_freq is within compound model estimated credible intervals. Defaults to True.



###### Returns

<code>np.array</code>
:   Pairs of high density interval limits, allele name, and posterior mean.
    as a 4 by n array.
    In alphabetical order of alleles, regardless of input order.
    This way it matches the output of combineAF().



    
### Function `compare_estimates` {#HLAfreq.HLAfreq_pymc.compare_estimates}



    
> `def compare_estimates(AFtab, hdi, datasetID)`


Does the defaul estimate of <code>allele\_freq</code> sit within the compound
model's estimated credible intervals? If not, print warnings.


###### Args

**```AFtab```** :&ensp;<code>pd.DataFrame</code>
:   Table of allele frequency data


**```hdi```** :&ensp;<code>np.array</code>
:   Pairs of high density interval limits, allele name, and posterior mean from compound model.


**```datasetID```** :&ensp;<code>str</code>, optional
:   Unique identifier column for study.





-----
Generated by *pdoc* 0.8.1 (<https://pdoc3.github.io>).
