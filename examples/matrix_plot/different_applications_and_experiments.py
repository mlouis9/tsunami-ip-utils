"""
Matrices With a Different Set of Applications and Experiments
=============================================================
This example demonstrates the use of the `matrix_plot` function to visualize matrices with a different set of applications and 
experiments. This means that the diagonal of the matrix does not necessarily represent an application/experiment compared with
itself, and so contribution plots aren't shown.
"""

# %%
# The input is exactly as described in :ref:`sphx_glr_auto_examples_matrix_plot_matrix_of_contribution_plots.py`, except that the
# applications and experiments lists include different SDF files.

from tsunami_ip_utils.viz.viz import matrix_plot
from tsunami_ip_utils.viz.plot_utils import generate_plot_objects_array_from_contributions
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-00{i}.sdf' for i in range(1,4) ]
experiment_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-002-00{i}S.sdf' for i in range(1,4) ]

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(application_files, experiment_files, variance=True)

labels = {
    'applications': [ application_file.name for application_file in application_files ],
    'experiments': [ experiment_file.name for experiment_file in experiment_files ],
}

plot_objects_array = generate_plot_objects_array_from_contributions(
    contributions_nuclide,
    integral_index_name='(Δk/k)^2'
)

fig = matrix_plot(plot_objects_array, plot_type='interactive', labels=labels)
fig.show()

# %%
# Now the plot can also be saved to an image

fig = matrix_plot(plot_objects_array, plot_type='interactive')
fig.to_image( EXAMPLES / '_static' / 'different_applications_and_experiments_matrix.png' )

# sphinx_gallery_thumbnail_path = '../../examples/_static/different_applications_and_experiments_matrix.png'