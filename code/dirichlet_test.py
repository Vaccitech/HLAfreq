"""
Tests with the dirichlet distribution
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import dirichlet

# define a Dirichlet Distributions with parameters alpha
alpha = np.array([0.4, 5, 15])  # specify concentration parameters
dd = dirichlet(alpha)

dd.mean()
dd.var()
dd.entropy()

dd.pdf([.6,.3,.1])

# Visualise a dirichlet distribution as k normal distriubtions?
# or let k=3 and plot a 3d surface like https://en.wikipedia.org/wiki/Dirichlet_distribution#/media/File:LogDirichletDensity-alpha_0.3_to_alpha_2.0.gif

# update a dirichlet distribution with new info

# show the reduction in variance

# Compare addition of n agreeing samples and n disagreeing samples
