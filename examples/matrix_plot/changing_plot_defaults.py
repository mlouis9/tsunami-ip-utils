"""
Changing Plot Defaults
======================
This example demonstrates how to change the default plot settings for the matrix plot. This includes changing the plot
type on the diagonals, changing the plot type for the off-diagonals.
"""

from tsunami_ip_utils.viz.viz import matrix_plot
from tsunami_ip_utils.viz.plot_utils import generate_plot_objects_array_from_contributions
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-00{i}.sdf' for i in range(1,4) ]
experiment_files = application_files

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(application_files, experiment_files)

labels = {
    'applications': [ application_file.name for application_file in application_files ],
    'experiments': [ experiment_file.name for experiment_file in experiment_files ],
}

# %% 
# Changing the Diagonal Plot Type
# -------------------------------
# On the diagonals, where the application and experiment are the same, the plot is a pie chart of contributions to the
# nuclear-data indcued uncertainty, as shown in :ref:`sphx_glr_auto_examples_visualizing_contributions.py`. However, other
# types of plots, such as static bar plots can be used instead. This is done by specifying the ``diagonal_type`` in 
# the :func:`tsunami_ip_utils.viz.plot_utils.generate_plot_objects_array_from_contributions` function. Note that any plot type
# available in :func:`tsunami_ip_utils.viz.viz.contribution_plot` can be used.

# %%
# Static Pie Chart Diagonals
# ^^^^^^^^^^^^^^^^^^^^^^^^^^
# For a static pie chart, the ``diagonal_type`` should be set to ``'pie'``

plot_objects_array_nuclide = generate_plot_objects_array_from_contributions(
    contributions_nuclide, 
    '(Δk/k)^2', 
    diagonal_type='pie'
)
fig = matrix_plot(plot_objects_array_nuclide, plot_type='interactive', labels=labels)
fig.show()

# %%
# Static Bar Chart Diagonals
# ^^^^^^^^^^^^^^^^^^^^^^^^^^
# For a static pie chart, the ``diagonal_type`` should be set to ``'bar'``

plot_objects_array_nuclide = generate_plot_objects_array_from_contributions(
    contributions_nuclide, 
    '(Δk/k)^2', 
    diagonal_type='bar'
)
fig = matrix_plot(plot_objects_array_nuclide, plot_type='interactive', labels=labels)
fig.show()