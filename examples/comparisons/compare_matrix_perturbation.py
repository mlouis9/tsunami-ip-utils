"""
Matrix Perturbation Correlation Comparison
==========================================
This example demonstrates how to generate a matrix comparison plot for the perturbation correlation method (described in 
:ref:`sphx_glr_auto_examples_matrix_plot_matrix_of_perturbation_plots.py` and 
:ref:`sphx_glr_auto_examples_correlation_plot_perturbation_correlation_plot.py`). This example compares the calculated correlation
correlation coefficient to the TSUNAMI-IP calculated :math:`c_k` value for a set of SDFs computed from the HMF series of critical
experiments [Bess2019]_. **Note** before running this example, please ensure that the necessary cross section libraries are cached
by running the :ref:`sphx_glr_auto_examples_correlation_plot_perturbation_correlation_plot.py` example.)
"""

# %%
# Getting TSUNAMI-IP :math:`c_k` Values
# -------------------------------------
# First, we need to get the TSUNAMI-IP :math:`c_k` values for the experiment and application SDFs we want to compare. We can do this
# by using the :func:`tsunami_ip_utils.integral_indices.get_integral_indices` function, as shown in 
# :ref:`sphx_glr_auto_examples_integral_index_reader.py`.

from tsunami_ip_utils.integral_indices import get_integral_indices
from paths import EXAMPLES

application_sdfs = [ EXAMPLES / 'data' / 'example_sdfs' / 'u235-dummy' / f'sphere_model_{i}.sdf' for i in range(1, 4) ]
experiment_sdfs = application_sdfs

# Get the TSUNAMI-IP integral indices
coverx_library = '56groupcov7.1'
integral_indices = get_integral_indices(application_sdfs, experiment_sdfs, coverx_library=coverx_library)
c_k = integral_indices['c_k']

# %%
# Generating the Comparison
# -------------------------
# A matrix comparison plot comparing the computed correlation coefficient to the TSUNAMI-IP :math:`c_k` values can be generated using
# the :func:`tsunami_ip_utils.comparisons.correlation_comparison` function. Note this function is used to compare any method for computing
# any of the available TSUNAMI-IP integral indices (:math:`E` and :math:`c_k`). 

from pathlib import Path
from tsunami_ip_utils.comparisons import correlation_comparison

# Define paths used for generating cross section perturbations
multigroup_library = EXAMPLES / 'data' / 'dummy_56_v7.1'
perturbation_factors = Path("~/codes/SCALE-6.3.1/data/perturb/56n.v7.1")

# %%
# 50 Points
# ---------
# First, we'll generate a comparison for 50 points

num_perturbations = 50
comparisons, fig = correlation_comparison(
    integral_index_matrix=c_k,
    integral_index_name='c_k',
    application_files=application_sdfs, 
    experiment_files=experiment_sdfs, 
    method='perturbation',
    base_library=multigroup_library,
    perturbation_factors=perturbation_factors,
    num_perturbations=num_perturbations,
)

fig.show()
print(comparisons)

# %%
# From this we can see an interactive matrix plot showing the perturbation correlation plots with the computed pearson and spearman
# correlation coefficients and the TSUNAMI-IP computed :math:`c_k` values. Note that any plots disagreeing more than 5% are highlighted
# in red. In addition, we get a multi-index pandas dataframe showing the comparisons (which can be written to excel as in
# :ref:`sphx_glr_auto_examples_comparisons_compare_E.py`).
# Note that all of the plot objects are generated using the defaults, and the matrix plot is also generated using the defaults.
# If you want to customize the plot, just pass in ``plot_object_kwargs`` or ``matrix_plot_kwargs`` to the function.

# %%
# 200 Points
# ----------
# Now, using a slightly larger number of points: 200 (remember 1000 is the maximum), we can generate another comparison, but this time
# using static matplotlib plots and including labels

num_perturbations = 200
comparisons, fig = correlation_comparison(
    integral_index_matrix=c_k,
    integral_index_name='c_k',
    application_files=application_sdfs, 
    experiment_files=experiment_sdfs, 
    method='perturbation',
    base_library=multigroup_library,
    perturbation_factors=perturbation_factors,
    num_perturbations=num_perturbations,
)

fig.show()
print(comparisons)

# %%
# Note that for these examples, adding more poins hardly changes the agreement, because the TSUNAMI_IP :math:`c_k` values are ~1, only
# a small number of points are needed to get a small uncertainty on the correlation coefficient (see :ref:`the technical manual <sec-pearson-coefficient>`
# for details).