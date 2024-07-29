"""This module contains the :class:`InteractiveMatrixPlot` class, which is a wrapper around a Dash app that displays a matrix of
interactive plots. The matrix is constructed from a 2D numpy array of plot objects, where each plot object is an instance of
either :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` or 
:class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend`. The matrix is displayed in a Dash app, where each plot object is
displayed in a separate cell. The user can interact with the plots by clicking on the legend to hide or show traces. The user
can also save the state of the interactive matrix plot to a pickle file and load it back later. The :func:`load_interactive_matrix_plot`
function is a convenience function that loads an interactive matrix plot from a saved state pickle file."""

from __future__ import annotations
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
from .scatter_plot import InteractiveScatterLegend
from .pie_plot import InteractivePieLegend
import webbrowser
import os
import sys
import threading
from .plot_utils import _find_free_port
import pickle
import tsunami_ip_utils
from typing import Union, List
import tsunami_ip_utils.config as config

# Style constants
GRAPH_STYLE = {
    'flex': '1',
    'minWidth': '800px',
    'height': '500px',
    'padding': '10px',
    'borderRight': '1px solid black',
    'borderBottom': '1px solid black',
    'borderTop': '0px',
    'borderLeft': '0px'
}

def _create_app(external_stylesheets):
    return dash.Dash(__name__, external_stylesheets=external_stylesheets)

def _create_column_headers(num_cols: int) -> None:
    """Create column headers for the matrix plot. Each column header is a div element with the text 'Application i' where i is
    the column index. The column headers are styled to be centered and have a border on the right and bottom.
    
    Parameters
    ----------
    num_cols
        Number of columns in the matrix plot."""
    return [html.Div(
        f'Application {i+1}', 
        style={
            'flex': '1', 
            'minWidth': '800px', 
            'textAlign': 'center', 
            'padding': '10px', 
            'borderRight': '1px solid black', 
            'borderBottom': '1px solid black', 
            'display': 'flex', 
            'alignItems': 'center', 
            'justifyContent': 'center'
        }
    ) for i in range(num_cols)]

def _create_row_label(i: int) -> html.Div:
    """Create a row label for the matrix plot. The row label is a div element with the text 'Experiment i' where i is the row
    index. The row label is styled to be centered and have a border on the right and bottom. The text is rotated -90 degrees to
    make it vertical.
    
    Parameters
    ----------
    i
        Row index of the matrix plot.
        
    Returns
    -------
        A div element representing the row label."""
    return html.Div(
        html.Span(
            f'Experiment {i+1}',
            style={
                'display': 'block',
                'overflow': 'visible',
                'transform': 'rotate(-90deg)',
                'transformOrigin': 'center',
                'whiteSpace': 'nowrap',
            }
        ), 
        style={
            'flex': 'none',
            'width': '50px', 
            'textAlign': 'center', 
            'marginRight': '0', 
            'padding': '10px', 
            'borderRight': '1px solid black', 
            'borderBottom': '1px solid black', 
            'display': 'flex', 
            'alignItems': 'center', 
            'justifyContent': 'center'
        }
    )

def _create_plot_element(i: int, j: int, plot_object: Union[InteractiveScatterLegend, InteractivePieLegend, go.Figure]
                        ) -> Union[dcc.Graph, html.Iframe]:
    """Create a plot element based on the plot object. If the plot object is an instance of 
    :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend`, or a ``plotly.graph_objects.Figure``, 
    the plot element is a ``dcc.Graph`` element. If the plot object is an instance of
    :class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend`, the plot element is an ``html.Iframe`` element.
    
    """
    if isinstance(plot_object, InteractiveScatterLegend):
        graph_id = f"interactive-scatter-{i}-{j}"
        return dcc.Graph(id=graph_id, figure=plot_object.fig, style=GRAPH_STYLE)
    elif isinstance(plot_object, InteractivePieLegend):
        with plot_object._app.test_client() as client:
            response = client.get('/')
            html_content = response.data.decode('utf-8')
            return html.Iframe(srcDoc=html_content, style=GRAPH_STYLE)
    else:
        return dcc.Graph(figure=plot_object, style=GRAPH_STYLE)

def _create_update_figure_callback(app, graph_id, app_instance):
    @app.callback(
        Output(graph_id, 'figure'),
        Input(graph_id, 'restyleData'),
        State(graph_id, 'figure')
    )
    def update_figure_on_legend_click(restyleData, current_figure_state):
        if restyleData and 'visible' in restyleData[0]:
            current_fig = go.Figure(current_figure_state)

            # Get the index of the clicked trace
            clicked_trace_index = restyleData[1][0]

            # Get the name of the clicked trace
            clicked_trace_name = current_fig.data[clicked_trace_index].name

            # Update excluded isotopes based on the clicked trace
            if restyleData[0]['visible'][0] == 'legendonly' and clicked_trace_name not in app_instance.excluded_isotopes:
                app_instance.excluded_isotopes.append(clicked_trace_name)
            elif restyleData[0]['visible'][0] == True and clicked_trace_name in app_instance.excluded_isotopes:
                app_instance.excluded_isotopes.remove(clicked_trace_name)

            # Update DataFrame based on excluded isotopes
            updated_df = app_instance.df.copy()
            updated_df = updated_df[~updated_df['Isotope'].isin(app_instance.excluded_isotopes)]

            # Recalculate the regression and summary statistics
            app_instance.add_regression_and_stats(updated_df)

            # Update trace visibility based on excluded isotopes
            for trace in app_instance.fig.data:
                if trace.name in app_instance.excluded_isotopes:
                    trace.visible = 'legendonly'
                else:
                    trace.visible = True

            return app_instance.fig

        return dash.no_update

def _generate_layout(app: dash.Dash, rows: List[html.Div]) -> None:
    """Generate the layout of the Dash app. The layout consists of an H1 element with the title 'Matrix of Plots', followed by
    a div element containing the rows of the matrix plot. The rows are displayed in a flex column with horizontal scrolling.
    The layout also includes a JavaScript script that listens for window resize events and resizes the Plotly plots
    accordingly.
    
    Parameters
    ----------
    app
        Dash app object.
    rows
        List of div elements representing the rows of the matrix plot."""
    app.layout = html.Div([
        html.H1("Matrix of Plots", style={'textAlign': 'center', 'marginLeft': '121px'}),
        html.Div(rows, style={'display': 'flex', 'flexDirection': 'column', 'width': '100%', 'overflowX': 'auto'}),
        html.Script("""
        window.addEventListener('resize', function() {
            const graphs = Array.from(document.querySelectorAll('.js-plotly-plot'));
            graphs.forEach(graph => {
                Plotly.Plots.resize(graph);
            });
        });
        """)
    ])

class InteractiveMatrixPlot:
    """Interactive matrix plot class that displays a matrix of interactive plots in a Dash app. The matrix is constructed from
    a 2D numpy array of plot objects, where each plot object is an instance of either
    :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` or
    :class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend`. The matrix is displayed in a Dash app, where each plot object is
    displayed in a separate cell. The user can interact with the plots by clicking on the legend to hide or show traces. The user
    can also save the state of the interactive matrix plot to a pickle file and load it back later."""
    _app: dash.Dash
    """Dash app object that displays the interactive matrix plot."""
    _plot_objects_array: np.ndarray
    """2D numpy array of plot objects to be displayed in the matrix plot."""
    def __init__(self, app: dash.Dash, plot_objects_array: np.ndarray) -> None:
        """Initialize the InteractiveMatrixPlot object with the Dash app and the 2D numpy array of plot objects.
        
        Parameters
        ----------
        app
            Dash app object that displays the interactive matrix plot.
        plot_objects_array
            2D numpy array of plot objects to be displayed in the matrix plot."""
        self._app = app
        self._plot_objects_array = plot_objects_array
    
    def _open_browser(self, port: int) -> None:
        """Open the browser to display the interactive matrix plot.
        
        Parameters
        ----------
        port
            Port number of the Flask server."""
        print(f"Now running at http://localhost:{port}/")
        webbrowser.open(f"http://localhost:{port}/")
        pass

    def show(self, open_browser: bool=True, silent: bool=False):
        """Start the Flask server and open the browser to display the interactive sunburst chart
        
        Parameters
        ----------
        open_browser
            Whether to open the browser automatically to display the chart.
        silent
            Whether to suppress Flask's startup and runtime messages."""
        # Suppress Flask's startup and runtime messages by redirecting them to dev null
        log = open(os.devnull, 'w')
        sys.stderr = log
        if silent:
            sys.stdout = log

        port = _find_free_port()
        if not config.generating_docs:
            if open_browser:
                threading.Timer(1, self._open_browser(port)).start()
            self._app.run(host='localhost', port=port)
    
    def save_state(self, filename: Union[str, Path]) -> None:
        """Save the state of the interactive matrix plot to a pickle file. The state includes the 2D numpy array of plot objects
        and the types of the plot objects.
        
        Parameters
        ----------
        filename
            Name of the pickle file to save the state to. The file extension should be ``'.pkl'``."""
        # Serialize interactive plots in the plot objects array
        self.plot_types = np.empty_like(self._plot_objects_array, dtype=object)
        for i, row in enumerate(self._plot_objects_array):
            for j, plot_object in enumerate(row):
                if isinstance(plot_object, InteractiveScatterLegend):
                    self._plot_objects_array[i,j] = plot_object.save_state()
                    self.plot_types[i,j] = "InteractiveScatterLegend"
                elif isinstance(plot_object, InteractivePieLegend):
                    self._plot_objects_array[i,j] = plot_object.save_state()
                    self.plot_types[i,j] = "InteractivePieLegend"

        with open(filename, 'wb') as f:
            pickle.dump( ( self._plot_objects_array, self.plot_types ) , f)

    @classmethod
    def load_state(self, filename: Union[str, Path]) -> InteractiveMatrixPlot:
        """Loads an interactive matrix plot from a saved state pickle file.
        
        Parameters
        ----------
        filename
            Name of the pickle file to load the state from.
            
        Returns
        -------
            An reserialized instance of the :class:`InteractiveMatrixPlot` class."""
        with open(filename, 'rb') as f:
            plot_objects_array, plot_types = pickle.load(f)
            # Reserialize interactive scatter legends
            for i, row in enumerate(plot_objects_array):
                for j, plot_object in enumerate(row):
                    if plot_types[i,j] == "InteractiveScatterLegend":
                        plot_objects_array[i,j] = InteractiveScatterLegend.load_state(data_dict=plot_object)
                    elif plot_types[i,j] == "InteractivePieLegend":
                        plot_objects_array[i,j] = InteractivePieLegend.load_state(data_dict=plot_object)

        return _interactive_matrix_plot(plot_objects_array)
    

def load_interactive_matrix_plot(filename):
    """Loads an interactive matrix plot from a saved state pickle file. This function is purely for convenience and is a
    wrapper of the :meth:`InteractiveMatrixPlot.load_state` class method"""
    return InteractiveMatrixPlot.load_state(filename)


def _interactive_matrix_plot(plot_objects_array: np.ndarray) -> InteractiveMatrixPlot:
    """Create an interactive matrix plot from a 2D numpy array of plot objects. This function creates a Dash app that displays
    the matrix plot. The matrix is constructed from the plot objects array, where each plot object is an instance of either
    :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` or :class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend`.
    The matrix is displayed in a Dash app, where each plot object is displayed in a separate cell. The user can interact with
    the plots by clicking on the legend to hide or show traces.

    Parameters
    ----------
    plot_objects_array
        2D numpy array of plot objects to be displayed in the matrix plot.

    Returns
    -------
        An instance of the :class:`InteractiveMatrixPlot` class that wraps the Dash app displaying the matrix plot.
    """
    current_directory = Path(__file__).parent
    external_stylesheets = [str(current_directory / 'css' / 'matrix_plot.css')]
    app = _create_app(external_stylesheets)

    num_rows = plot_objects_array.shape[0]
    num_cols = plot_objects_array.shape[1]

    column_headers = _create_column_headers(num_cols)
    header_row = html.Div([html.Div('', style={'flex': 'none', 'width': '71px', 'borderBottom': '1px solid black'})] + column_headers, style={'display': 'flex'})

    rows = [header_row]
    for i in range(num_rows):
        row = [_create_row_label(i)]
        for j in range(num_cols):
            plot_object = plot_objects_array[i, j]
            plot_element = _create_plot_element(i, j, plot_object) if plot_object else html.Div('Plot not available', style=GRAPH_STYLE)
            row.append(plot_element)
            if isinstance(plot_object, InteractiveScatterLegend):
                _create_update_figure_callback(app, f"interactive-scatter-{i}-{j}", plot_object)
        rows.append(html.Div(row, style={'display': 'flex'}))

    _generate_layout(app, rows)
    return InteractiveMatrixPlot(app, plot_objects_array)