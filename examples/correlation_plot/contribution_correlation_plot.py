"""
Contribution Correlation Plots
==============================
This example demonstrates the use of the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function to generate a contribution 
correlation plot.
"""

# %%
# First, we import the necessary modules, and generate teh uncertainty contributions for the desired application and experiment
# pair (as described in :ref:`sphx_glr_auto_examples_uncertainty_contributions.py`).

from tsunami_ip_utils.viz.viz import correlation_plot
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES
import matplotlib.pyplot as plt

application_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-001-001.sdf' ]
experiment_files = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-004-009.sdf' ]

contributions_nuclide, contributions_nuclide_reaction = get_uncertainty_contributions(
    application_files, 
    experiment_files,
)

# %%
# Nuclide-Wise Contributions
# --------------------------
# Next, we generate the correlation plot using the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function. First we start
# with the nuclide-wise contributions.

# %%
# Matplotlib Scatter Plot
# ^^^^^^^^^^^^^^^^^^^^^^^
# The default plot type of the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function is a matplotlib scatter plot.

application_contributions = contributions_nuclide['application'][0]
experiment_contributions  = contributions_nuclide['experiment'][0]
fig, axs = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='$(Δk/k)^2$'
)

# Note that the generated figure object is special, and has attributes relating to the regression and correlation coefficients
# that are shown in the plot. These are accessible as follows:
print(fig.statistics)
print(fig.regression)

plt.show()

# %%
# From the plot you can clearly see the linear regression, the pearson correlation coefficient, and the spearman rank correlation
# coefficient (more detail given in :ref:`the technical manual <sec-correlation-coefficient>`).

# %%
# Plotly Interactive Scatter Plot
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Matplotlib plots are nice for reports, but it's often useful to be able to interactively explore the data, and see which points
# correspond to which nuclides, and see the x and y ordinates of each point. These are all made possible with Plotly plots. To
# enable Plotly plots, simply set the ``plot_type`` argument to ``'interactive_scatter'``. Note that now, only a single ``fig`` object is
# returned, rather than ``fig`` and ``axs`` objects as in the matplotlib case.

fig = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='(Δk/k)^2', # Note that Plotly plots don't support pylatex in labels, so we just use a plain string
    plot_type='interactive_scatter'
)

# Like before, the generated figure object also has information about the regression and correlation coefficients
print(fig.statistics)
print(fig.regression)

fig.show()

# %%
# Note that ``fig.show()`` starts a local server that serves the web-based Plotly plot (that should pop up in a new browser tab).

# %%
# Plotly Interactive Scatter Plot with Interactive Legend
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# The previous plot allows you to toggle the visibility of specific nuclides, but unfortunately doesn't update the regression
# or summary statistics to reflect that change. To do that, you must pass the ``interactive_legend`` argument as ``True``.

fig = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='(Δk/k)^2', # Note that Plotly plots don't support pylatex in labels, so we just use a plain string
    plot_type='interactive_scatter',
    interactive_legend=True,
)
fig.show()

# %%
# **Note** if you're viewing this example via the online documentation, the regression and summary statistics interactivity
# won't work, because this requires client-side python, and the plots shown above are merely html snapshots. To see the full
# interactivity, you must run the example locally.

# %%
# Nuclide-Reaction-Wise Contributions
# -----------------------------------
# Next, we generate the correlation plot using the :func:`tsunami_ip_utils.viz.viz.correlation_plot` function for the 
# nuclide-reaction-wise contributions. The input is exactly the same, but the contributions are on a nuclide-reaction-wise
# basis.

# %%
# Matplotlib Scatter Plot
# ^^^^^^^^^^^^^^^^^^^^^^^

application_contributions = contributions_nuclide_reaction['application'][0] # Note that the contributions are nuclide_reaction now
experiment_contributions  = contributions_nuclide_reaction['experiment'][0]
fig, axs = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='$(\\Delta k/k)^2$'
)
plt.show()

# %%
# Plotly Interactive Scatter Plot
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fig = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='(Δk/k)^2', # Note that Plotly plots don't support pylatex in labels, so we just use a plain string
    plot_type='interactive_scatter'
)
fig.show()

# %%
# Plotly Interactive Scatter Plot with Interactive Legend
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fig = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='(Δk/k)^2', # Note that Plotly plots don't support pylatex in labels, so we just use a plain string
    plot_type='interactive_scatter',
    interactive_legend=True,
)
fig.show()