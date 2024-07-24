"""
Visualizing Contributions
=========================

This example illustrates the use of the various tools for visualizing contributions to a TSUNAMI-IP integral index (such as
:math:`E`, or even contributions to :math:`\Delta k/k`). This example shows the visualization of :math:`\Delta k/k` contributions,
but the same tools can be used for visualizing contributions to :math:`E`.
"""

# %%
# Isotope-Wise Contributions
# --------------------------
# First, we demonstrate the tools for visualizing `isotope-wise` contributions to the nuclear data induced uncertainty, 
# :math:`\Delta k/k`. First we need to get the contributions (further details are provided in
# :ref:`sphx_glr_auto_examples_uncertainty_contributions.py`).

from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_filenames = [ f"{EXAMPLES}/data/example_sdfs/MCT/MIX-COMP-THERM-001-001.sdf" ]
uncertainty_contributions_nuclide, uncertainty_contributions_nuclide_reaction = \
      get_uncertainty_contributions(application_filenames)

# %%
# The isotope-wise contributions can be plotted using the :func:`tsunami_ip_utils.viz.contribution_plot` function.
# This function takes a contribution list and creates a plot of the specified type. First, we must extract a list of
# contributions for a specific application from the :func:`tsunami_ip_utils.integral_indices.get_uncertainty_contributions`
# output.

contributions = uncertainty_contributions_nuclide['application'][0]
# %%
# Bar Chart
# ^^^^^^^^^
# The contributions can then be plotted as a matplotlib bar chart.

from tsunami_ip_utils.viz.viz import contribution_plot
import matplotlib.pyplot as plt

fig, axs = contribution_plot(contributions, plot_type='bar',integral_index_name='Δk/k')
plt.tight_layout()
plt.show()

# %%
# Pie Chart
# ^^^^^^^^^
# The contributions can also be plotted as a matplotlib pie chart.

fig, axs = contribution_plot(contributions, plot_type='pie',integral_index_name='Δk/k')
plt.tight_layout()
plt.show()

# %%
# Interactive Pie Chart
# ^^^^^^^^^^^^^^^^^^^^^
# The contributions can also be plotted as an interactive plotly pie chart.

fig = contribution_plot(contributions, plot_type='interactive_pie', integral_index_name='Δk/k', interactive_legend=False)
fig.show()

# %%
# Interactive Pie Chart With Interactive Legend
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# The contributions can also be plotted as an interactive plotly pie chart with an interactive legend. This feature has
# yet to be formally implemented in Plotly Express (although an issue is open `here <https://github.com/plotly/plotly.py/issues/4692>`_ ).

fig = contribution_plot(contributions, plot_type='interactive_pie', integral_index_name='Δk/k', interactive_legend=True)
fig.show()

# %%
# Reaction-Wise Contributions
# ---------------------------
# The isotope-reaction-wise contributions can be plotted similarly.

# %%
# Bar Chart
# ^^^^^^^^^

contributions = uncertainty_contributions_nuclide_reaction['application'][0]

fig, axs = contribution_plot(contributions, plot_type='bar',integral_index_name='Δk/k')
plt.tight_layout()
plt.show()

# %%
# Pie Chart
# ^^^^^^^^^

fig, axs = contribution_plot(contributions, plot_type='pie',integral_index_name='Δk/k')
plt.tight_layout()
plt.show()

# %%
# Interactive Pie Chart
# ^^^^^^^^^^^^^^^^^^^^^

fig = contribution_plot(contributions, plot_type='interactive_pie', integral_index_name='Δk/k', interactive_legend=False)
fig.show()

# %%
# Interactive Pie Chart With Interactive Legend
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fig = contribution_plot(contributions, plot_type='interactive_pie', integral_index_name='Δk/k', interactive_legend=True)
fig.show()

# %%
# Test3