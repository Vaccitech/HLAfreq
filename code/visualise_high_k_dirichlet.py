"""
Visualise high k Dirichlet
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.stats import dirichlet, beta
import pandas as pd
import math
import code.scrapeAF as scrapeAF

#####################
#                   #
#   Define alpha    #
#                   #
#####################

# Load allele frequency data
AFtab = pd.read_csv(f"data/example/outlier_raw.csv")
# Reduce resolution to 2 fields and collapse alleles
AFtab = scrapeAF.decrease_resolution(AFtab, 2)
# Combined allle frequency
# noncomplete = scrapeAF.incomplete_studies(AFtab)
# complete_mask = AFtab.apply(lambda x: (x.population, x.loci) not in noncomplete.index, axis=1)
# AFtab = AFtab[complete_mask]
caf = scrapeAF.combineAF(AFtab)
# Recreate the Dirichlet distribution
naivealpha = scrapeAF.default_prior(len(caf.allele))

# Dirichlet distribution concentration parameters 
# based on allele frequency data
alpha = naivealpha + caf.c

######################
#                    #
# View k=3 Dirichlet #
#                    #
######################
if len(alpha) == 3:
    dd = dirichlet(alpha)
    # Make data for k=3.
    k1 = np.arange(0, 1, 0.05)
    k2 = np.arange(0, 1, 0.05)
    k1, k2 = np.meshgrid(k1, k2)
    k3 = 1 - k1 - k2
    k3[k3<0] = np.nan

    pdf = np.zeros(k1.shape)
    for i in range(k1.shape[0]):
        for j in range(k1.shape[1]):
            pdf[i,j] =  dd.pdf([k1[i,j], k2[i,j], k3[i,j]])

    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(k1, k2, pdf, cmap=cm.plasma,
        linewidth=0, antialiased=False)
    ax.set_title('PDF')
else:
    print( "Surface plot only works with k==3")

######################
#                    #
#   View Dirichlet   #
#       as beta      #
#                    #
######################
    
def plotAFprob(alpha, caf=pd.DataFrame(), AFtab=pd.DataFrame(), datasetID="population", log=False, psteps=1000, ncol=2):
    """Plot the (log) posterior density function of all frequencies
        for all alleles based on supplied alpha for Dirichlet
        distribution. Options for adding empirical values
        and plotting the log pdf.

    Args:
        alpha (list): Dirichlet concentration parameters, can be prior + caf.c
        caf (pd.DataFrame, optional): Combined allele frequence data produced by scrapeAF.combineAF(). Defaults to pd.DataFrame().
        AFtab (pd.DataFrame, optional): The uncombined allele frequency data used by scrapeAF.combinedAF(). You must use the same dataframe as this function doesn't have the error checking that scrapeAF.combineAF() has. Defaults to pd.DataFrame().
        datasetID (str, optional): The column used to define datasets. Defaults to "population".
        log (bool, optional): Plot log pdf instead of pdf? Defaults to False.
        psteps (int, optional): Number of increments in pdf calculation, higher values make smoother plots. Defaults to 1000.
        ncol (int, optional): How many columns to arrange subplots in. Defaults to 2.
    """
    # Get beta parameters for each k in Dirichlet
    ab = [(a,sum(alpha)-a) for a in alpha]
    pline = np.linspace(0,1,psteps)
    if not AFtab.empty:
        # Format the populatoin allele frequencies
        df = AFtab.copy()
        df = scrapeAF.unmeasured_alleles(df, datasetID=datasetID)
        df = df.sort_values('allele')
        assert all(df.groupby('population').allele.apply(list).apply(lambda x: x == caf.allele.tolist())), "Alleles not matching between AFtab and caf"
    fig, axs = plt.subplots(math.ceil(len(alpha)/ncol), ncol)
    for i,x in enumerate(alpha):
        subplotselector = i//ncol, i%ncol
        a,b = ab[i]
        bd = beta(a,b)
        if log:
            pdf = [bd.logpdf(p) for p in pline]
        else:
            pdf = [bd.pdf(p) for p in pline]
        ax = axs[subplotselector]
        if not AFtab.empty:
            # Add the empirical allele frequency for each population
            for af in df.groupby('population').allele_freq.apply(list).apply(lambda x: x[i]):
                # x is the reported allele freq, y is it's pdf with scatter
                if log:
                    ax.scatter(af, bd.logpdf(af)*(0.95+np.random.random()/5))
                else:
                    ax.scatter(af, bd.pdf(af)*(0.95+np.random.random()/5))
        ax.plot(pline, pdf)
        # Annotate with alpha
        ax.text(0.5, max(pdf)/2, f"alpha {round(alpha[i])}")
        if not caf.empty:
            # Plot the combined average
            ax.axvline(caf.allele_freq[i], color="black", ls="--")
        # fig.tight_layout()
    plt.show()

# plotAFprob(alpha)
# plotAFprob(alpha, caf)
plotAFprob(alpha, caf, AFtab)
plotAFprob(alpha, caf, AFtab, log=True)

