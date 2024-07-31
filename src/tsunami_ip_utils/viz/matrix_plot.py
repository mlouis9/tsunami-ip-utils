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
from typing import Union, List, Dict, Optional
import tsunami_ip_utils.config as config
import time, requests
import signal
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
import multiprocessing
import tempfile
from tsunami_ip_utils.viz.plot_utils import _capture_html_as_image
import base64
import io
import matplotlib
import matplotlib.pyplot as plt
import plotly.io as pio

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
    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
    return app

def _create_column_headers(num_cols: int, labels: Optional[List[str]]=None) -> None:
    """Create column headers for the matrix plot. Each column header is a div element with the text 'Application i' where i is
    the column index. The column headers are styled to be centered and have a border on the right and bottom.
    
    Parameters
    ----------
    num_cols
        Number of columns in the matrix plot.
    labels
        List of labels for the columns of the matrix plot. Default is ``None``."""
    if labels is None:
        div = [html.Div(
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
    else:
        div = [html.Div(
            label, 
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
        ) for label in labels]
    return div

def _create_row_label(i: int, label: Optional[str]=None) -> html.Div:
    """Create a row label for the matrix plot. The row label is a div element with the text 'Experiment i' where i is the row
    index. The row label is styled to be centered and have a border on the right and bottom. The text is rotated -90 degrees to
    make it vertical.
    
    Parameters
    ----------
    i
        Row index of the matrix plot.
    label
        An optional label for the row (to use instead of ``Experiment i``). Default is ``None``.
        
    Returns
    -------
        A div element representing the row label."""
    if label is None:
        div = html.Div(
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
    else:
        div = html.Div(
            html.Span(
                label,
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

    return div

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
    elif isinstance(plot_object, tuple):
        if isinstance(plot_object[0], matplotlib.figure.Figure):
            # Convert matplotlib figure to PNG image
            buf = io.BytesIO()
            plot_object[0].tight_layout()
            plot_object[0].savefig(buf, format='png')
            buf.seek(0)
            plt.close(plot_object[0])  # Close the figure to prevent it from being captured
            encoded_image = base64.b64encode(buf.read()).decode("utf-8")
            return html.Img(src='data:image/png;base64,{}'.format(encoded_image), style=GRAPH_STYLE)
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
            if restyleData[0]['visible'][0] == 'legendonly' and clicked_trace_name not in app_instance._excluded_isotopes:
                app_instance._excluded_isotopes.append(clicked_trace_name)
            elif restyleData[0]['visible'][0] == True and clicked_trace_name in app_instance._excluded_isotopes:
                app_instance._excluded_isotopes.remove(clicked_trace_name)

            # Update DataFrame based on excluded isotopes
            updated_df = app_instance.df.copy()
            updated_df = updated_df[~updated_df['Isotope'].isin(app_instance._excluded_isotopes)]

            # Recalculate the regression and summary statistics
            app_instance._add_regression_and_stats(updated_df)

            # Update trace visibility based on excluded isotopes
            for trace in app_instance.fig.data:
                if trace.name in app_instance._excluded_isotopes:
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
        html.Div(rows, id='matrix-plot',style={'display': 'flex', 'flexDirection': 'column', 'width': '100%', 'overflowX': 'auto'}),
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
    def __init__(self, app: dash.Dash, plot_objects_array: np.ndarray, labels: Optional[Dict[str, List]]=None) -> None:
        """Initialize the InteractiveMatrixPlot object with the Dash app and the 2D numpy array of plot objects.
        
        Parameters
        ----------
        app
            Dash app object that displays the interactive matrix plot.
        plot_objects_array
            2D numpy array of plot objects to be displayed in the matrix plot.
        labels
            Dictionary of lists containing the labels for the rows and columns of the matrix plot. Keys are ``'applications'`` 
            and ``'experiments'``. Default is ``None``.
        """
        self._app = app
        self._plot_objects_array = plot_objects_array
        self._labels = labels
    
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
            self._app.run(host='localhost', port=port, debug=False)
    
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

    def write_html(self, filename: Optional[Union[str, Path]]=None) -> None:
        """Write the HTML content of the interactive matrix plot to a file.

        Parameters
        ----------
        filename
            Name of the HTML file to write the content to. The file extension should be `'.html'`.
        """
        # Create a new Dash app instance
        app = _create_app(external_stylesheets=[])

        num_rows = self._plot_objects_array.shape[0]
        num_cols = self._plot_objects_array.shape[1]

        column_headers = _create_column_headers(num_cols, self._labels['applications']) if self._labels else _create_column_headers(num_cols)
        header_row = html.Div([html.Div('', style={'flex': 'none', 'width': '71px', 'borderBottom': '1px solid black'})] + column_headers, style={'display': 'flex'})

        rows = [header_row]
        for i in range(num_rows):
            row = [_create_row_label(i, self._labels['experiments'][i])] if self._labels else [_create_row_label(i)]
            for j in range(num_cols):
                plot_object = self._plot_objects_array[i, j]
                if isinstance(plot_object, InteractiveScatterLegend):
                    plot_element = html.Iframe(srcDoc=plot_object.write_html(), style=GRAPH_STYLE)
                elif isinstance(plot_object, InteractivePieLegend):
                    plot_element = html.Iframe(srcDoc=plot_object.write_html(), style=GRAPH_STYLE)
                elif isinstance(plot_object, tuple):
                    if isinstance(plot_object[0], matplotlib.figure.Figure):
                        # Convert matplotlib figure to PNG image
                        buf = io.BytesIO()
                        plot_object[0].tight_layout()
                        plot_object[0].savefig(buf, format='png')
                        buf.seek(0)
                        plt.close(plot_object[0])  # Close the figure to prevent it from being captured
                        encoded_image = base64.b64encode(buf.read()).decode("utf-8")
                        plot_element = html.Img(src='data:image/png;base64,{}'.format(encoded_image), style=GRAPH_STYLE)
                else:
                    html_string = pio.to_html(plot_object, full_html=False)
                    plot_element = html.Iframe(srcDoc=html_string, style=GRAPH_STYLE)
                row.append(plot_element)
            rows.append(html.Div(row, style={'display': 'flex'}))

        _generate_layout(app, rows)

        port = _find_free_port()
        def run_app():
            log = open(os.devnull, 'w')
            sys.stdout = log
            sys.stderr = log
            app.run_server('localhost', port, debug=False)

        # Start the app in a separate thread
        process = multiprocessing.Process(target=run_app)
        process.start()

        # Set up a headless Chrome browser using Selenium
        chrome_options = Options()
        chrome_options.add_argument("--headless")  # Run Chrome in headless mode
        driver = webdriver.Chrome(options=chrome_options)

        try:
            # Load the Dash app in the headless browser
            driver.get(f"http://localhost:{port}")

            # Wait for the app to finish loading (you can adjust the sleep time if needed)
            time.sleep(1)

            # Capture the fully rendered HTML content
            html_content = driver.page_source

            # Write the HTML content to the specified file
            if filename is not None:
                with open(filename, "w") as file:
                    file.write(html_content)
            else:
                return html_content

        finally:
            # Quit the headless browser
            driver.quit()

            # Stop the Dash server
            requests.post(f"http://localhost:{port}/shutdown")
            process.terminate()
            process.join()
    
    def to_image(self, filename: Union[str, Path]):
        """Write the initial state of the InteractiveMatrixPlot to an image file. This function saves the matrix plot as an 
        image by using selenium webdriver.
        
        Parameters
        ----------
        filename
            Name of the image file to save the matrix plot to. The file extension should be ``'.png'``, ``.jpg``, etc."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            html_filename = f.name + '.html'
        self.write_html(html_filename)
        _capture_html_as_image(html_filename, filename, matrix=True)
        os.remove(html_filename)


def load_interactive_matrix_plot(filename):
    """Loads an interactive matrix plot from a saved state pickle file. This function is purely for convenience and is a
    wrapper of the :meth:`InteractiveMatrixPlot.load_state` class method"""
    return InteractiveMatrixPlot.load_state(filename)


def _interactive_matrix_plot(plot_objects_array: np.ndarray, labels: Optional[Dict[str, List]]) -> InteractiveMatrixPlot:
    """Create an interactive matrix plot from a 2D numpy array of plot objects. This function creates a Dash app that displays
    the matrix plot. The matrix is constructed from the plot objects array, where each plot object is an instance of either
    :class:`tsunami_ip_utils.viz.scatter_plot.InteractiveScatterLegend` or :class:`tsunami_ip_utils.viz.pie_plot.InteractivePieLegend`.
    The matrix is displayed in a Dash app, where each plot object is displayed in a separate cell. The user can interact with
    the plots by clicking on the legend to hide or show traces.

    Parameters
    ----------
    plot_objects_array
        2D numpy array of plot objects to be displayed in the matrix plot.
    labels
        Dictionary of lists containing the labels for the rows and columns of the matrix plot. Keys are ``'applications'`` and 
        ``'experiments'``. Default is ``None``. Most commonly, this is the sdf filename for the given application/experiment.

    Returns
    -------
        An instance of the :class:`InteractiveMatrixPlot` class that wraps the Dash app displaying the matrix plot.
    """
    current_directory = Path(__file__).parent
    external_stylesheets = [str(current_directory / 'css' / 'matrix_plot.css')]
    app = _create_app(external_stylesheets)
    @app.server.route('/shutdown', methods=['POST'])
    def shutdown():
        """Function to shutdown the server"""
        os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
        return 'Server shutting down...'

    num_rows = plot_objects_array.shape[0]
    num_cols = plot_objects_array.shape[1]

    column_headers = _create_column_headers(num_cols, labels['applications']) if labels else _create_column_headers(num_cols)
    header_row = html.Div([html.Div('', style={'flex': 'none', 'width': '71px', 'borderBottom': '1px solid black'})] + column_headers, style={'display': 'flex'})

    rows = [header_row]
    for i in range(num_rows):
        if labels is not None:
            row = [_create_row_label(i, labels['experiments'][i])]
        else:
            row = [_create_row_label(i)]
        for j in range(num_cols):
            plot_object = plot_objects_array[i, j]
            plot_element = _create_plot_element(i, j, plot_object) if plot_object else html.Div('Plot not available', style=GRAPH_STYLE)
            row.append(plot_element)
            if isinstance(plot_object, InteractiveScatterLegend):
                _create_update_figure_callback(app, f"interactive-scatter-{i}-{j}", plot_object)
        rows.append(html.Div(row, style={'display': 'flex'}))

    _generate_layout(app, rows)
    return InteractiveMatrixPlot(app, plot_objects_array, labels=labels)