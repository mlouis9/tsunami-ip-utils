import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
from .base_plotter import Plotter
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import plotly.graph_objects as go
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import pandas as pd
import webbrowser
import os
import signal
import pickle
import threading
from .plot_utils import find_free_port
import sys, os, signal
import threading
import webbrowser
import sys
from plotly.graph_objs import Figure
from pyparsing import *
import re
from uncertainties import ufloat

def replace_spearman_and_pearson(text, new_pearson, new_spearman):
    """Replaces the Spearman and Pearson values in the given text with the new values provided using regex."""
    
    new_pearson = f"{new_pearson:1.6f}"
    new_spearman = f"{new_spearman:1.6f}"

    # Replace Pearson value using regex
    pearson_pattern = r"(Pearson: <b>)[0-9.]+(</b>)"
    text = re.sub(pearson_pattern, r"\g<1>" + new_pearson + r"\2", text)

    # Replace Spearman value using regex
    spearman_pattern = r"(Spearman: <b>)[0-9.]+(</b>)"
    text = re.sub(spearman_pattern, r"\g<1>" + new_spearman + r"\2", text)
    
    return text


def update_percent_difference(text, reference_value):
    # Regex to find the TSUNAMI-IP c_k value
    ck_pattern = r"TSUNAMI-IP c_k: <b>([\d\.]+)\+/-([\d\.]+)</b>"
    # Find the TSUNAMI-IP c_k value in the text
    match = re.search(ck_pattern, text)
    if not match:
        raise ValueError("TSUNAMI-IP c_k value not found in the text")
    
    # Parse the nominal value and the uncertainty
    nominal_value, uncertainty = match.groups()
    tsunami_ck = ufloat(float(nominal_value), float(uncertainty))
    
    # Calculate the percent difference
    if reference_value == 0 or np.isnan(reference_value):
        percent_difference = ufloat(0, 0)
    else:
        percent_difference = (tsunami_ck - reference_value) / reference_value * 100
    
    # Format the new percent difference with uncertainty
    new_percent_diff_text = f"<b>{percent_difference.nominal_value:.2f}+/-{percent_difference.std_dev:.2f}</b>%"
    
    # Regex to replace the Percent Difference in the text
    # Escaping special characters used in the regex pattern
    updated_text = re.sub(r"Percent Difference: <b>[\d\.\+\-/]+</b>+%", f"Percent Difference: {new_percent_diff_text}", text)
    
    return updated_text, percent_difference

class EnhancedPlotlyFigure(Figure):
    """This class wraps a plotly express figure object (intended for a scatter plot) and adds additional metadata for the
    summary statistics and linear regression data. This class is intended to be used with the ``InteractiveScatterPlotter``
    class."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Directly set the attributes using object's __setattr__ to bypass Plotly's checks
        object.__setattr__(self, 'statistics', None)
        object.__setattr__(self, 'regression', None)

    def __setattr__(self, name, value):
        if name in ['statistics', 'regression']:
            # Handle custom attributes internally
            object.__setattr__(self, name, value)
        else:
            # Use the super class's __setattr__ for all other attributes
            super().__setattr__(name, value)

class ScatterPlot(Plotter):
    """This class exists to add some additional functionality for calculating regressions and summary statistics that's
    common to all types of scatter plots, interactive or otherwise"""
    def get_summary_statistics(self, x, y):
        """Calculates the Pearson correlation coefficient, Spearman rank correlation coefficient, and linear regression
        parameters for the given x and y datasets. The linear regression parameters are the slope and intercept of the
        regression line. The Pearson and Spearman coefficients are also stored in the class instance as 'pearson' and
        'spearman' respectively. The slope and intercept are stored as 'slope' and 'intercept' respectively. The linear
        regression is stored as 'regression'"""
        self.regression = stats.linregress(x, y)
        self.pearson    = stats.pearsonr(x, y).statistic
        self.spearman   = stats.spearmanr(x, y).statistic
        self.slope      = self.regression.slope
        self.intercept  = self.regression.intercept

        # If the figure has been plotted (and is an enhanced plot which supports adding this metadata), add the regression 
        # and correlation statistics to the figure
        is_metadata_plot = isinstance(getattr(self, 'fig', None), EnhancedPlotlyFigure)
        if hasattr(self, 'fig') and is_metadata_plot:
            self.fig.statistics = {
                'pearson': self.pearson,
                'spearman': self.spearman
            }
            self.fig.regression = {
                'slope': self.slope,
                'intercept': self.intercept,
            }

        # Now create teh summary statistics text for figure annotation
        pearson_text = f"Pearson: <b>{self.pearson:1.6f}</b>"
        spearman_text = f"Spearman: <b>{self.spearman:1.6f}</b>"
        self.summary_stats_text = f"{pearson_text} {spearman_text}"


class ScatterPlotter(ScatterPlot):
    def __init__(self, integral_index_name, nested, plot_redundant=False, **kwargs):
        self.nested = nested
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contribution_pairs, isotopes, reactions):
        self.fig, self.axs = plt.subplots()

        # Extract the x and y values from the contribution pairs
        application_points        = [ contribution[0].n for contribution in contribution_pairs ]
        application_uncertainties = [ contribution[0].s for contribution in contribution_pairs ]
        experiment_points         = [ contribution[1].n for contribution in contribution_pairs ]
        experiment_uncertainties  = [ contribution[1].s for contribution in contribution_pairs ]

        self.fig = plt.errorbar(application_points, experiment_points, xerr=application_uncertainties, \
                               yerr=experiment_uncertainties, fmt='.', capsize=5)
        
        # Linear regression
        self.get_summary_statistics(application_points, experiment_points)

        # Plot the regression line
        x = np.linspace(min(application_points), max(application_points), 100)
        y = self.slope * x + self.intercept
        self.axs.plot(x, y, 'r', label='Linear fit')

        self.axs.text(0.05, 0.95, self.summary_stats_text, transform=self.axs.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

        self.style()

    def get_plot(self):
        return self.fig, self.axs
        
    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
    
    def style(self):
        title_text = f'Contributions to {self.index_name}'
        self.axs.set_title(title_text)
        self.axs.set_ylabel(f"Experiment {self.index_name} Contribution")
        self.axs.set_xlabel(f"Application {self.index_name} Contribution")
        self.axs.grid()


def load_interactive_scatter_plot(filename):
    """Loads an interactive scatter plot from a saved state file. This function is purely for convenience and is a
    wrapper of the InteractiveScatterLegend.load_state method"""
    return InteractiveScatterLegend.load_state(filename)


class InteractiveScatterPlotter(ScatterPlot):
    def __init__(self, integral_index_name, nested, plot_redundant=False, **kwargs):
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = False
        self.nested = nested
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contribution_pairs, isotopes, reactions):
        self.fig = make_subplots()

        # Extract isotope and reaction pairs from the given list of isotopes and reactions
        df = self._create_scatter_data(contribution_pairs, isotopes, reactions)

        hover_data_dict = {
            'Isotope': True  # Always include Isotope
        }

        if 'Reaction' in df.columns:
            hover_data_dict['Reaction'] = True  # Include Reaction only if it exists

        # Create scatter plot with error bars using Plotly Express
        self.fig = px.scatter(
            df, 
            x=f'Application {self.index_name} Contribution', 
            y=f'Experiment {self.index_name} Contribution',
            error_x='Application Uncertainty', 
            error_y='Experiment Uncertainty',
            color='Isotope',
            labels={
                "color": "Isotope"
            },
            title=f'Contributions to {self.index_name}',
            hover_data=hover_data_dict
        )

        # Wrap the plotly express figure in a MetadataPlotly object
        self.fig = EnhancedPlotlyFigure(self.fig.to_dict())

        self.add_regression_and_stats(df)

        # Now style the plot
        self.style()

        if self.interactive_legend:
            self.fig = InteractiveScatterLegend(self, df)

    def add_regression_and_stats(self, df):
        # Calculate the linear regression and correlation statistics
        self.get_summary_statistics(df[f'Application {self.index_name} Contribution'], \
                                    df[f'Experiment {self.index_name} Contribution'])

        # Prepare data for the regression line
        x_reg = np.linspace(df[f'Application {self.index_name} Contribution'].min(), 
                            df[f'Application {self.index_name} Contribution'].max(), 100)
        y_reg = self.slope * x_reg + self.intercept

        # Convert self.fig.data to a list for mutability
        current_traces = list(self.fig.data)

        # Remove existing regression line if it exists
        traces_to_keep = [trace for trace in current_traces if not trace.name.startswith('Regression Line')]

        # Set the modified list of traces back to the figure
        self.fig.data = tuple(traces_to_keep)

        # Add new linear regression to the plot
        self.fig.add_trace(go.Scatter(x=x_reg, y=y_reg, mode='lines', 
                                    name=f'Regression Line y={self.slope:1.4E}x + {self.intercept:1.4E}'))

        # Remove existing annotation if it exists
        if hasattr(self.fig, 'layout') and hasattr(self.fig.layout, 'annotations'):
            if len(self.fig.layout.annotations) != 0:
                # Get the text of the old annotation
                old_annotation_text = \
                    [ann for ann in self.fig.layout.annotations if ann.text.startswith('Pearson')][0].text
                annotation_text = replace_spearman_and_pearson(old_annotation_text, self.pearson, self.spearman)
                bordercolor = '#444'

                self.fig.layout.annotations = [ann for ann in self.fig.layout.annotations if not ann.text.startswith('Pearson')]
                if "TSUNAMI-IP" in annotation_text:
                    try:
                        annotation_text, percent_difference = update_percent_difference(annotation_text, self.pearson)
                    except Exception as e:
                        print(f"Error updating percent difference: {e}")
                    if abs(percent_difference.nominal_value) > 5:
                        bordercolor = 'red'
            else:
                bordercolor = '#444'
                annotation_text = self.summary_stats_text
        else:
            bordercolor = '#444'
            annotation_text = self.summary_stats_text

        # Add correlation statistics to the plot
        self.fig.add_annotation(
            x=0.05, xref="paper", 
            y=0.95, yref="paper",
            text=annotation_text,
            bordercolor=bordercolor,
            showarrow=False, 
            font=dict(size=12),
            align='left',
            bgcolor="white", 
            opacity=0.8
        )

    def _create_scatter_data(self, contribution_pairs, isotopes, reactions):
        data = {
            f'Application {self.index_name} Contribution': [cp[0].n for cp in contribution_pairs],
            f'Experiment {self.index_name} Contribution': [cp[1].n for cp in contribution_pairs],
            'Application Uncertainty': [cp[0].s for cp in contribution_pairs],
            'Experiment Uncertainty': [cp[1].s for cp in contribution_pairs],
            'Isotope': [],
        }

        # Add nuclides and reactions (if they exist) to the data dictionary
        if reactions == []:
            for isotope in isotopes:
                data['Isotope'].append(isotope)
        else:
            data['Reaction'] = []
            for isotope in isotopes:
                for reaction in reactions:
                    data['Isotope'].append(isotope)
                    data['Reaction'].append(reaction)

        # Now filter out (0,0) points, which don't contribute to either the application or the experiment, these are
        # usually chi, nubar, or fission reactions for nonfissile isotopes that are added for consistency with the set
        # of reactions only
        data = { key: [val for val, app, exp in zip(data[key], data[f'Application {self.index_name} Contribution'], \
                    data[f'Experiment {self.index_name} Contribution']) if app != 0 or exp != 0] for key in data }

        return pd.DataFrame(data)

    def add_to_subplot(self, fig, position):
        for trace in self.fig.data:
            fig.add_trace(trace, row=position[0], col=position[1])

        # Transfer annotations
        if hasattr(self.fig, 'layout') and hasattr(self.fig.layout, 'annotations'):
            for ann in self.fig.layout.annotations:
                # Adjust annotation references to new subplot
                new_ann = ann.update(xref=f'x{position[1]}', yref=f'y{position[1]}')
                fig.add_annotation(new_ann, row=position[0], col=position[1])
        return fig

    def get_plot(self):
        return self.fig
    
    def style(self):
        title_text = f'Contributions to {self.index_name}'
        self.fig.update_layout(title_text=title_text, title_x=0.5)  # 'title_x=0.5' centers the title


class InteractivePerturbationScatterPlotter(ScatterPlot):
    def __init__(self, **kwargs):
        pass

    def create_plot(self, points):
        self.fig = make_subplots()

        # Extract isotope and reaction pairs from the given list of isotopes and reactions
        df = pd.DataFrame({
            'Application': [point[0].n for point in points],
            'Experiment': [point[1].n for point in points],
            'Application Uncertainty': [point[0].s for point in points],
            'Experiment Uncertainty': [point[1].s for point in points]
        })

        # Create scatter plot with error bars using Plotly Express
        self.fig = px.scatter(
            df, 
            x=f'Application', 
            y=f'Experiment',
            error_x='Application Uncertainty', 
            error_y='Experiment Uncertainty',
            title=f'Correlation Plot',
        )

        # Wrap the plotly express figure in a MetadataPlotly object
        self.fig = EnhancedPlotlyFigure(self.fig.to_dict())

        self.add_regression_and_stats(df)

        # Now style the plot
        self.style()

    def add_regression_and_stats(self, df):
        # Calculate the linear regression and correlation statistics
        self.get_summary_statistics(df[f'Application'], df[f'Experiment'])

        # Prepare data for the regression line
        x_reg = np.linspace(df[f'Application'].min(), df[f'Application'].max(), 100)
        y_reg = self.slope * x_reg + self.intercept

        # Convert self.fig.data to a list for mutability
        current_traces = list(self.fig.data)

        # Remove existing regression line if it exists
        traces_to_keep = [trace for trace in current_traces if not trace.name.startswith('Regression Line')]

        # Set the modified list of traces back to the figure
        self.fig.data = tuple(traces_to_keep)

        # Add new linear regression to the plot
        self.fig.add_trace(go.Scatter(x=x_reg, y=y_reg, mode='lines', 
                                    name=f'Regression Line y={self.slope:1.4E}x + {self.intercept:1.4E}'))

        # Add correlation statistics to the plot
        self.fig.add_annotation(
            x=0.05, xref="paper", 
            y=0.95, yref="paper",
            text=self.summary_stats_text,
            showarrow=False,
            align='left',
            font=dict(size=12),
            bgcolor="white", 
            opacity=0.8
        )

    def add_to_subplot(self, fig, position):
        for trace in self.fig.data:
            fig.add_trace(trace, row=position[0], col=position[1])
        return fig

    def get_plot(self):
        return self.fig
    
    def style(self):
        pass


class InteractiveScatterLegend(InteractiveScatterPlotter):
    def __init__(self, interactive_scatter_plot, df):
        self.interactive_scatter_plot = interactive_scatter_plot
        self.fig = interactive_scatter_plot.fig
        self.index_name = interactive_scatter_plot.index_name
        self.df = df
        self.excluded_isotopes = []  # Keep track of excluded isotopes
        self.app = dash.Dash(__name__)
        self.app.layout = html.Div([
            dcc.Graph(id='interactive-scatter', figure=self.fig, style={'height': '100vh'})
        ], style={'margin': 0})
        self.setup_callbacks()

    def setup_callbacks(self):
        @self.app.callback(
            Output('interactive-scatter', 'figure'),
            Input('interactive-scatter', 'restyleData'),
            State('interactive-scatter', 'figure')
        )
        def update_figure_on_legend_click(restyleData, current_figure_state):
            if restyleData and 'visible' in restyleData[0]:
                current_fig = go.Figure(current_figure_state)

                # Get the index of the clicked trace
                clicked_trace_index = restyleData[1][0]

                # Get the name of the clicked trace
                clicked_trace_name = current_fig.data[clicked_trace_index].name

                # Update excluded isotopes based on the clicked trace
                if restyleData[0]['visible'][0] == 'legendonly' and clicked_trace_name not in self.excluded_isotopes:
                    self.excluded_isotopes.append(clicked_trace_name)
                elif restyleData[0]['visible'][0] == True and clicked_trace_name in self.excluded_isotopes:
                    self.excluded_isotopes.remove(clicked_trace_name)

                # Update DataFrame based on excluded isotopes
                updated_df = self.df.copy()
                updated_df = updated_df[~updated_df['Isotope'].isin(self.excluded_isotopes)]

                # Recalculate the regression and summary statistics
                self.add_regression_and_stats(updated_df)

                # Update trace visibility based on excluded isotopes
                for trace in self.fig.data:
                    if trace.name in self.excluded_isotopes:
                        trace.visible = 'legendonly'
                    else:
                        trace.visible = True

                return self.fig

            return dash.no_update

        @self.app.server.route('/shutdown', methods=['POST'])
        def shutdown():
            os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
            return 'Server shutting down...'

    def update_figure_on_legend_click(self, restyleData, current_figure_state):
        if restyleData and 'visible' in restyleData[0]:
            current_fig = go.Figure(current_figure_state)

            # Get the index of the clicked trace
            clicked_trace_index = restyleData[1][0]

            # Get the name of the clicked trace
            clicked_trace_name = current_fig.data[clicked_trace_index].name

            # Update excluded isotopes based on the clicked trace
            if restyleData[0]['visible'][0] == 'legendonly' and clicked_trace_name not in self.excluded_isotopes:
                self.excluded_isotopes.append(clicked_trace_name)
            elif restyleData[0]['visible'][0] == True and clicked_trace_name in self.excluded_isotopes:
                self.excluded_isotopes.remove(clicked_trace_name)

            # Update DataFrame based on excluded isotopes
            updated_df = self.df.copy()
            updated_df = updated_df[~updated_df['Isotope'].isin(self.excluded_isotopes)]

            # Create a new _InteractiveScatterPlotter instance with the updated data
            updated_plotter = InteractiveScatterPlotter(self.index_name, self.interactive_scatter_plot.nested)
            updated_plotter.create_plot(self.interactive_scatter_plot._create_scatter_data(updated_df), updated_df['Isotope'].unique(), updated_df['Reaction'].unique() if 'Reaction' in updated_df.columns else [])

            # Update the current figure with the new traces and layout
            current_fig.data = updated_plotter.fig.data
            current_fig.layout = updated_plotter.fig.layout

            # Update trace visibility based on excluded isotopes
            for trace in current_fig.data:
                if trace.name in self.excluded_isotopes:
                    trace.visible = 'legendonly'
                else:
                    trace.visible = True

            return current_fig

        return dash.no_update

    def show(self):
        port = find_free_port()
        # Function to open the browser
        def open_browser():
            if not os.environ.get("WERKZEUG_RUN_MAIN"):
                print(f"Now running at http://localhost:{port}/")
                webbrowser.open(f"http://localhost:{port}/")

        # Silence the Flask development server logging
        log = open(os.devnull, 'w')
        # sys.stdout = log
        sys.stderr = log

        # Disable Flask development server warning
        os.environ['FLASK_ENV'] = 'development'

        # JavaScript code to detect when the tab or window is closed
        self.app.index_string = '''
        <!DOCTYPE html>
        <html>
            <head>
                {%metas%}
                <title>{%title%}</title>
                {%favicon%}
                {%css%}
            </head>
            <body style="margin: 0;">
                {%app_entry%}
                <footer>
                    {%config%}
                    {%scripts%}
                    <script type="text/javascript">
                        window.addEventListener("beforeunload", function (e) {
                            var xhr = new XMLHttpRequest();
                            xhr.open("POST", "/shutdown", false);
                            xhr.send();
                        });
                    </script>
                    {%renderer%}
                </footer>
            </body>
        </html>
        '''

        # Timer to open the browser shortly after the server starts
        threading.Timer(1, open_browser).start()

        self.app.run_server(debug=False, host='localhost', port=port)

    def save_state(self, filename=None):
        state = {
            'fig': self.fig.to_dict(),
            'df': self.df.to_dict(),
            'excluded_isotopes': self.excluded_isotopes,
            'index_name': self.index_name,
            'nested': self.interactive_scatter_plot.nested
        }
        if filename is None:
            return state
        else:
            with open(filename, 'wb') as f:
                pickle.dump(state, f)

    @classmethod
    def load_state(cls, filename=None, data_dict=None):
        if filename is None and data_dict is None:
            raise ValueError("Either a filename or a data dictionary must be provided")
        if filename is not None:
            with open(filename, 'rb') as f:
                state = pickle.load(f)
        else:
            state = data_dict

        # Recreate the _InteractiveScatterPlotter instance from the saved state
        fig = go.Figure(state['fig'])
        index_name = state['index_name']
        nested = state['nested']
        interactive_scatter_plot = InteractiveScatterPlotter(index_name, nested)
        interactive_scatter_plot.fig = fig

        # Recreate the InteractiveScatterLegend instance from the saved state
        instance = cls(interactive_scatter_plot, pd.DataFrame.from_dict(state['df']))
        instance.excluded_isotopes = state['excluded_isotopes']

        # Update trace visibility based on excluded isotopes
        for trace in instance.fig.data:
            if trace.name in instance.excluded_isotopes:
                trace.visible = 'legendonly'
            else:
                trace.visible = True

        return instance

    def write_html(self, filename):
        # Utilize Plotly's write_html to save the current state of the figure
        self.fig.write_html(filename)
