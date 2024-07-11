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
from .plot_utils import find_free_port
import pickle

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

def create_app(external_stylesheets):
    return dash.Dash(__name__, external_stylesheets=external_stylesheets)

def create_column_headers(num_cols):
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

def create_row_label(i):
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

def create_plot_element(i, j, plot_object):
    if isinstance(plot_object, InteractiveScatterLegend):
        graph_id = f"interactive-scatter-{i}-{j}"
        return dcc.Graph(id=graph_id, figure=plot_object.fig, style=GRAPH_STYLE)
    elif isinstance(plot_object, InteractivePieLegend):
        with plot_object.app.test_client() as client:
            response = client.get('/')
            html_content = response.data.decode('utf-8')
            return html.Iframe(srcDoc=html_content, style=GRAPH_STYLE)
    else:
        return dcc.Graph(figure=plot_object, style=GRAPH_STYLE)

def create_update_figure_callback(app, graph_id, app_instance):
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

def generate_layout(app, rows):
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
    def __init__(self, app, plot_objects_array):
        self.app = app
        self.plot_objects_array = plot_objects_array
    
    def open_browser(self, port):
        print(f"Now running at http://localhost:{port}/")
        webbrowser.open(f"http://localhost:{port}/")
        pass

    def show(self, open_browser=True, silent=False):
        """Start the Flask server and open the browser to display the interactive sunburst chart
        
        Parameters
        ----------
        - open_browser: bool, whether to open the browser automatically to display the chart
        - silent: bool, whether to suppress Flask's startup and runtime messages"""
        # Suppress Flask's startup and runtime messages by redirecting them to dev null
        log = open(os.devnull, 'w')
        # sys.stdout = log
        sys.stderr = log

        port = find_free_port()
        if open_browser:
            threading.Timer(1, self.open_browser(port)).start()
        self.app.run(host='localhost', port=port)
    
    def save_state(self, filename):
        # Serialize interactive plots in the plot objects array
        self.plot_types = np.empty_like(self.plot_objects_array, dtype=object)
        for i, row in enumerate(self.plot_objects_array):
            for j, plot_object in enumerate(row):
                if isinstance(plot_object, InteractiveScatterLegend):
                    self.plot_objects_array[i,j] = plot_object.save_state()
                    self.plot_types[i,j] = "InteractiveScatterLegend"
                elif isinstance(plot_object, InteractivePieLegend):
                    self.plot_objects_array[i,j] = plot_object.save_state()
                    self.plot_types[i,j] = "InteractivePieLegend"

        with open(filename, 'wb') as f:
            pickle.dump( ( self.plot_objects_array, self.plot_types ) , f)

    @classmethod
    def load_state(self, filename):
        with open(filename, 'rb') as f:
            plot_objects_array, plot_types = pickle.load(f)
            # Reserialize interactive scatter legends
            for i, row in enumerate(plot_objects_array):
                for j, plot_object in enumerate(row):
                    if plot_types[i,j] == "InteractiveScatterLegend":
                        plot_objects_array[i,j] = InteractiveScatterLegend.load_state(data_dict=plot_object)
                    elif plot_types[i,j] == "InteractivePieLegend":
                        plot_objects_array[i,j] = InteractivePieLegend.load_state(data_dict=plot_object)

        return interactive_matrix_plot(plot_objects_array)
    

def load_interactive_matrix_plot(filename):
    """Loads an interactive matrix plot from a saved state pickle file. This function is purely for convenience and is a
    wrapper of the InteractiveScatterLegend.load_state method"""
    return InteractiveMatrixPlot.load_state(filename)


def interactive_matrix_plot(plot_objects_array: np.ndarray):
    current_directory = Path(__file__).parent
    external_stylesheets = [str(current_directory / 'css' / 'matrix_plot.css')]
    app = create_app(external_stylesheets)

    num_rows = plot_objects_array.shape[0]
    num_cols = plot_objects_array.shape[1]

    column_headers = create_column_headers(num_cols)
    header_row = html.Div([html.Div('', style={'flex': 'none', 'width': '71px', 'borderBottom': '1px solid black'})] + column_headers, style={'display': 'flex'})

    rows = [header_row]
    for i in range(num_rows):
        row = [create_row_label(i)]
        for j in range(num_cols):
            plot_object = plot_objects_array[i, j]
            plot_element = create_plot_element(i, j, plot_object) if plot_object else html.Div('Plot not available', style=GRAPH_STYLE)
            row.append(plot_element)
            if isinstance(plot_object, InteractiveScatterLegend):
                create_update_figure_callback(app, f"interactive-scatter-{i}-{j}", plot_object)
        rows.append(html.Div(row, style={'display': 'flex'}))

    generate_layout(app, rows)
    return InteractiveMatrixPlot(app, plot_objects_array)