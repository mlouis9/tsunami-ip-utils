"""
Saving and Loading Interactive Plots
====================================
This example shows how to save and reopen interactive plots created with ``tsunami_ip_utils``. This is useful for sharing plots
with collaborators or saving plots for later use.
"""

# %%
# Non-Custom Interactive Plots
# ----------------------------
# The interactive plots generated in :ref:`sphx_glr_auto_examples_visualizing_contributions.py` can all be saved using the ``write_html``
# method of the plotly figure object. This method writes the plot to an HTML file, which can be opened in a web browser to view the
# interactive plot. Even for :class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend` plots (which aren't plotly figures, but rather
# a custom wrapper around a plotly figure), the same ``write_html`` method can be used.
#
# Saving Plots
# ------------
# However, for plots that implement interactive python callbacks, such as :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend`
# plots, the interactivity cannot fully be saved to html because the callbacks are implemented in python (though the current static html state of
# the plot can be written to html with the :meth:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend.write_html` method).
# So, to save these plots, a different approach is taken. These plots are serialized as python pickle objects, which can be reloaded
# later by anyone with ``tsunami_ip_utils`` installed. To serialize an interactive scatter plot, the 
# :meth:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend.save_state` method is used. This method just takes a path to
# the output pickle file to save the plot state to.

# %%
# First we need to get some data to create an interactive scatter plot. We'll use the same code as in 
# :ref:`sphx_glr_auto_examples_correlation_plot_contribution_correlation_plot.py`

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

# Now get the application and experiment contributions
application_contributions = contributions_nuclide_reaction['application'][0]
experiment_contributions  = contributions_nuclide_reaction['experiment'][0]

fig = correlation_plot(
    application_contributions, 
    experiment_contributions, 
    integral_index_name='%(Δk/k)^2', # Note that Plotly plots don't support pylatex in labels, so we just use a plain string
    plot_type='interactive_scatter',
    interactive_legend=True,
)
fig.save_state( EXAMPLES / 'scatter_plot.pkl' )

# %%
# Loading Plots
# -------------
# Now, that we've serialized the plot to a ``.pkl`` file, we can reload it using the 
# :func:`tsunami_ip_utils.viz.scatter_plot.load_interactive_scatter_plot` function.

from tsunami_ip_utils.viz.scatter_plot import load_interactive_scatter_plot

fig = load_interactive_scatter_plot( EXAMPLES / 'scatter_plot.pkl' )
fig.show()

# %%
# and as the above shows, the full interactivity of the scatter plot (including the legend updates with the recalculated 
# summary statistics and linear regression) is retained.

# %%
# We can also save the plot to an image

fig = load_interactive_scatter_plot( EXAMPLES / 'scatter_plot.pkl' )
fig.to_image( EXAMPLES / '_static' / 'loaded_scatter_plot.png' )

# sphinx_gallery_thumbnail_path = '../examples/_static/loaded_scatter_plot.png'

# %%
# Other Interactive Plots
# -----------------------
# The only other polts that have to be serialized like above to retain their interactivity are matrix plots (i.e. those created by
# :func:`tsunami_ip_utils.viz.viz.matrix_plot``). These figure objects also have a ``save state`` method (i.e.
# :meth:`tsunami_ip_utils.viz.matrix_plot.InteractiveMatrixPlot.save_state`), which can be used to save the plot state to a pickle file.
# And these plots can be reloaded, likewise, using the :func:`tsunami_ip_utils.viz.matrix_plot.load_interactive_matrix_plot` function.
# An example is shown below for a matrix that is `very` time consuming to generate.

from tsunami_ip_utils.viz.matrix_plot import load_interactive_matrix_plot

fig = load_interactive_matrix_plot( EXAMPLES / 'data' / 'hmf_contribution_comparisons_nuclide_reaction.pkl' )
fig.show()