"""
Contribution Correlation Plots
==============================
This example demonstrates the use of the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function to generate a contribution 
correlation plot.
"""

# %%
# First, we import the necessary modules, and generate teh uncertainty contributions for the desired application and experiment
# pair (as described in :ref:`../sphx_glr_auto_examples_uncertainty_contributions.py`).

from tsunami_ip_utils.viz.viz import correlation_plot
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-001.sdf' ]
experiment_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-004-009.sdf' ]

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(application_files, experiment_files)

# %%
# Nuclide-Wise Contributions
# --------------------------
# Next, we generate the correlation plot using the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function. First we start
# with the nuclide-wise contributions.

application_contributions = contributions_nuclide['application'][0]
experiment_contributions = contributions_nuclide['experiment'][0]
fig, axs = correlation_plot(application_contributions, experiment_contributions)
fig.show()