"""
Matrix of Perturbation Correlation Plots
========================================
This example demonstrates how to generate a matrix plot of perturbation plots for a suite of applications and experiments. This method
is detailed in the :ref:`technical manual <sec-consistent-random-variables>`, and requires randomly sampling nuclear data libraries to
generate a scatter plot. It is recommended that the user refers to :ref:`sphx_glr_auto_examples_correlation_plot_perturbation_correlation_plot.py`
for more information on the perturbation correlation plots that compose the matrix plots shown herein before continuing with this
example.
"""

# %%
# Interactive Matrix Plots
# ------------------------
# Interactive matrix plots wrap several Plotly plots in an html grid that's served by a Dash app that allows users to interact
# with each of the individual plots. This is useful for comparing multiple applications and experiments in a single figure.
# Specifically, it is useful for viewing all combinations of a given set of SDF files in a single figure
#
# **Note** that by default the the correlation plots in the off-diagonal cells are 
# :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` plots, which are fully interactive, and dynamically update
# the calculated summary statistics and regression when points are excluded with the legend. However, this interactivity requires
# clientiside python to run, and so is not represented in the documentation. The off-diagonal plot types can be changed, and doing so
# is documented in :ref:`sphx_glr_auto_examples_matrix_plot_changing_plot_defaults.py`.

from paths import EXAMPLES

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'u235-dummy' / f'sphere_model_{i}.sdf' for i in range(1,4) ]
experiment_files = application_files

# %%
# Generate the Points Array
# -------------------------
# First, an array of points must be generated for each application-experiment pair (these are used to construct the perturbation plots).

from tsunami_ip_utils.perturbations import generate_points
from pathlib import Path

base_library = EXAMPLES / 'data' / 'dummy_56_v7.1'
perturbation_factors = Path("~/codes/SCALE-6.3.1/data/perturb/56n.v7.1")
num_perturbations = 100

points = generate_points(application_files, experiment_files, base_library, perturbation_factors, num_perturbations)

# %%
# Generate the Plot Objects
# -------------------------
# Now, the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function needs a numpy array of plot objects to generate the matrix plot.
# We can generate this array using the :func:`tsunami_ip_utils.viz.plot_utils.generate_plot_objects_array_from_perturbations` function.
# This function takes a dictionary of contributions and the name of the integral index, and returns a numpy array of plot objects
# that can be passed to the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function.

from tsunami_ip_utils.viz.plot_utils import generate_plot_objects_array_from_perturbations

plot_objects = generate_plot_objects_array_from_perturbations(points)

# %%
# Create the Matrix Plot
# ----------------------
# Finally, we can create the matrix plot using the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function. This function takes the
# numpy array of plot objects generated in the previous step, and generates a matrix plot of perturbation correlation plots.

from tsunami_ip_utils.viz.viz import matrix_plot

labels = {
    'applications': [ file.name for file in application_files ],
    'experiments': [ file.name for file in experiment_files ]
}
fig = matrix_plot(plot_objects, plot_type='interactive', labels=labels)
fig.show()

# %%
# The plot can be saved to an image by using the :meth:`tsunami_ip_utils.viz.matrix_plot.InteractiveMatrixPlot.to_image` method

fig = matrix_plot(plot_objects, plot_type='interactive', labels=labels)
fig.to_image( EXAMPLES / '_static' / 'perturbation_matrix.png' )

# sphinx_gallery_thumbnail_path = '../../examples/_static/perturbation_matrix.png'