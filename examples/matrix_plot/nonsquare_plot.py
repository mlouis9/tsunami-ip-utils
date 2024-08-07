"""
Non-Square Matrix Plots
=======================
This example demonstrates how to create a non-square matrix of correlation plots. This is useful for comparing a single application
(or perhaps a few) with a suite of experiments (like those from the ICSBEP handbook [Bess2019]_).
"""

# %%
# Generating the Plot
# -------------------
# The steps for creating the matrix plot are just as described in :ref:`sphx_glr_auto_examples_matrix_plot_matrix_of_contribution_plots.py`
# except that the application files and experiment files lists are not the same length. Recall, first the uncertainty contributions
# must be generated from the given SDF profiles using :func:`tsunami_ip_utils.integral_indices.get_uncertainty_contributions`. In this
# example, we will choose a single application from the MCT directory, and compare it to all other SDFS in the MCT directory.

from tsunami_ip_utils.viz.viz import matrix_plot
from tsunami_ip_utils.viz.plot_utils import generate_plot_objects_array_from_contributions
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES
import os

# Get the filenames of all of the SDF profiles in the MCT directory
all_mct_sdfs = os.listdir(EXAMPLES / 'data' / 'example_sdfs' / 'MCT')

# Generate the lists of application and experiment SDF file paths
application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-001.sdf' ]
experiment_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / filename for filename in all_mct_sdfs 
                    if filename != 'MIX-COMP-THERM-001-001.sdf' ]

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(application_files, experiment_files, variance=True)

# %%
# Then the matrix of plot objects (that is plotted using :func:`tsunami_ip_utils.viz.viz.matrix_plot`) can be generated using
# :func:`tsunami_ip_utils.viz.plot_utils.generate_plot_objects_array_from_contributions` (note we also label the matrix cells by the
# application and experiment sdf filenames by explicitly passing the ``labels`` dictionary as an argument).

labels = {
    'applications': [ application_file.name for application_file in application_files ],
    'experiments': [ experiment_file.name for experiment_file in experiment_files ],
}

plot_objects_array = generate_plot_objects_array_from_contributions(
    contributions_nuclide,
    integral_index_name='(Δk/k)^2'
)

# %%
# Finally, the matrix plot can be generated and displayed.

fig = matrix_plot(plot_objects_array, plot_type='interactive', labels=labels)
fig.show()

# %%
# Sorting the Matrix
# ------------------
# If you're comparins a single application to multiple experiments, it may be useful to sort the experiments by the computed correlation
# coefficient (pearson, spearman, etc.). Since a plot objects array is first generated, then plotted using the 
# :func:`tsunami_ip_utils.viz.viz.matrix_plot` function, we may arbitrarily reorder the plot objects array before plotting. However,
# we must be careful to reorder the labels dictionary in the same way as the plot objects array. The following code demonstrates how
# to sort the experiments by the pearson correlation coefficient.

import numpy as np

# Sort the plot objects array by the pearson correlation coefficient using numpy.argsort and get the indices representing the sorted
# ordering of plot objects.

sorted_indices = np.argsort([ plot_object.statistics['pearson'] for plot_object in plot_objects_array[:, 0] ])

# Note [:,0] indexes the fist column of the matrix plot, i.e. the column corresponding to the first application.

# %%
# Note that by default :func:`numpy.argsort` returns the indices that would sort the array in ascending order. If you want to sort in
# descending order, you can reverse the order of the indices as follows:

sorted_indices = sorted_indices[::-1]

# Sort the plot objects by the indices
sorted_plot_objects_array = plot_objects_array[sorted_indices]

# Sort the experiment labels by the indices
labels['experiments'] = [ labels['experiments'][index] for index in sorted_indices ]

# now plot the sorted matrix
fig = matrix_plot(sorted_plot_objects_array, plot_type='interactive', labels=labels)
fig.show()