"""
Simple Plotly Plot Example
==========================

This example demonstrates how to create a simple scatter plot using Plotly in Sphinx Gallery.
"""

# Code source: MyProject
# License: BSD 3 clause 

# sphinx_gallery_thumbnail_path = '_static/boi.png'

import plotly.graph_objects as go

def simple_plotly_scatter():
    """
    Generates a simple scatter plot using Plotly.

    Returns:
        fig: A Plotly figure object.
    """
    # Sample data
    x = [1, 2, 3, 4, 5]
    y = [1, 4, 9, 16, 25]

    # Create a Plotly figure
    fig = go.Figure(
        data=[go.Scatter(x=x, y=y, mode='markers+lines', name='Data')],
        layout=go.Layout(
            title='Simple Scatter Plot',
            xaxis_title='X Axis Label',
            yaxis_title='Y Axis Label'
        )
    )

    # Show plot (this will not be executed in the documentation build process)
    fig.show()

    return fig

###############################################################################
# Let's call the function and store the figure. Sphinx Gallery will use the
# custom scraper to handle this figure object.

fig = simple_plotly_scatter()