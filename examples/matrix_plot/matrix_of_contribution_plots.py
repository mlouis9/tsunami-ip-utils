"""
Matrix of Contribution Correlation Plots
========================================
The first type of (interactive) similarity plot that we will view in matrix form is the contribution correlation plot discussed
in :ref:`sphx_glr_auto_examples_correlation_plot_contribution_correlation_plot.py`. This shows the contribution correlation plot between all
pairs of applications and experiments in a (symmetric) matrix. First we consider the case of a square matrix where the set
of applications and experiments are the same.
"""

# %%
# Interactive Matrix Plots
# ------------------------
# Interactive matrix plots wrap several Plotly plots in an html grid that's served by a Dash app that allows users to interact
# with each of the individual plots. This is useful for comparing multiple applications and experiments in a single figure.
# Specifically, it is useful for viewing all combinations of a given set of SDF files in a single figure. 
# 
# **Note** that by default the the correlation plots in the off-diagonal cells are 
# :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` plots, which are fully interactive, and dynamically update
# the calculated summary statistics and regression when points are excluded with the legend. However, this interactivity requires
# clientiside python to run, and so is not represented in the documentation. The off-diagonal plot types can be changed, and doing so
# is documented in :ref:`sphx_glr_auto_examples_matrix_plot_changing_plot_defaults.py`.

from tsunami_ip_utils.viz.viz import matrix_plot
from tsunami_ip_utils.viz.plot_utils import generate_plot_objects_array_from_contributions
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-00{i}.sdf' for i in range(1,4) ]
experiment_files = application_files

# %%
# Now we need to get the uncertainty contributions for each experiment and each application using 
# :func:`tsunami_ip_utils.integral_indices.get_uncertainty_contributions` (described in 
# :ref:`sphx_glr_auto_examples_uncertainty_contributions.py`).

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(
    application_files, 
    experiment_files, 
    variance=True
)

# %%
# Nuclide-Wise Contributions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now, the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function needs a numpy array of plot objects to generate the matrix plot.
# We can generate this array using the :func:`tsunami_ip_utils.viz.plot_utils.generate_plot_objects_array_from_contributions` function.
# This function takes a dictionary of contributions and the name of the integral index, and returns a numpy array of plot objects
# that can be passed to the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function.

plot_objects_array_nuclide = generate_plot_objects_array_from_contributions(contributions_nuclide, '%(Δk/k)^2')

# %%
# Now the matrix plot can easily be generated

fig = matrix_plot(plot_objects_array_nuclide, plot_type='interactive')
fig.show()

# %%
# This plot is fully interactive. You can hold shift and scroll using the mousewheel to scroll horizontally, or just use the
# mousewheel alone to scroll vertically. Individual subfigures can also be saved by clicking the camera icon in the top right,
# or the entire plot can be saved by calling the ``to_image`` method on the figure object.

fig = matrix_plot(plot_objects_array_nuclide, plot_type='interactive')
fig.to_image( EXAMPLES / '_static' / 'matrix_plot.png' )

# sphinx_gallery_thumbnail_path = '../../examples/_static/matrix_plot.png'

# %%
# Note if the default, nondescriptive labels aren't preferred, an arbitrary dictionary of labels for the applications and
# experiments may be passed instead. Perhaps the most useful is the SDF file names.

labels = {
    'applications': [ application_file.name for application_file in application_files ],
    'experiments': [ experiment_file.name for experiment_file in experiment_files ],
}
fig = matrix_plot(plot_objects_array_nuclide, plot_type='interactive', labels=labels)
fig.show()

# %%
# Nuclide-Reaction-Wise Contributions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The same can be done for the nuclide-reaction-wise contributions

plot_objects_array_nuclide_reaction = generate_plot_objects_array_from_contributions(contributions_nuclide_reaction, '%(Δk/k)^2')
fig = matrix_plot(plot_objects_array_nuclide_reaction, plot_type='interactive', labels=labels)
fig.show()

# %%
# Static Matrix Plots
# -------------------
# Static matrix plots are used for creating an analogous figure using static matplotlib plots. This implementation would have
# no need for the html of the interactive matrix plot, but as it stands, the interactive plot can contain all static plots as
# well (and a mixture of static and interactive plots), and so this feature has not yet been implemented.