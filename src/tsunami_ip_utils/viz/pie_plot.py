"""Tools for creating pie charts of contributions to integral indices."""

from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from ._base_plotter import _Plotter
from plotly.subplots import make_subplots
import os, sys, signal
import threading
import webbrowser
from flask import Flask, render_template_string
import uuid
from .plot_utils import _find_free_port
import pickle
from typing import Tuple, Dict, Union, Optional
from uncertainties import ufloat
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import plotly
from pathlib import Path
from tsunami_ip_utils import config
import tempfile
from tsunami_ip_utils.viz.plot_utils import _capture_html_as_image

plt.rcParams['hatch.linewidth'] = 0.6

class _PiePlotter(_Plotter):
    """A class for creating static pie charts of contributions to integral indices."""
    _index_name: str
    """The name of the integral index whose contributions are being plotted (e.g. ``'E'`` or ``'c_k'``)."""
    _plot_redundant: bool
    """Whether to plot redundant/irrelevant reactions in the pie chart."""
    def __init__(self, integral_index_name, plot_redudant=False, **kwargs):
        """Initializes a pie chart of the contributions to the given integral index.
        
        Parameters
        ----------
        integral_index_name
            The name of the integral index whose contributions are to be plotted.
        plot_redundant
            Wether to include redundant/irrelevant reactions in the plot. NOTE: this only applies to nested plots, and
            only affects the plot title; it is expected that the provided data is consistent with the flag.
            
        Notes
        -----
        * Redundant reactions are defined as those which are derived from other reactions, e.g. 'total' and 'capture' reactions
          in SCALE.
        * Irrelevant reactions are defined as those which are not directly cross sections (but rather other nuclear data 
          parameters), e.g. 'chi' in SCALE.
        * A flag for including/excluding redundant/irrelevant reactions was provided since, if the user is expecting the
          the contributions to add up nicely, then the redundant reactions should be excluded, and if only cross sections are
          being considered, then the irrelevant reactions should be excluded.
        """
        self._index_name = integral_index_name
        self._plot_redundant = plot_redudant
    
    def _create_plot(self, contributions: Union[Dict[str, ufloat], Dict[str, Dict[str, ufloat]]], nested: bool):
        """Creates a pie chart of the given contributions to the integral index.
        
        Parameters
        ----------
        contributions
            * If ``nested`` is ``False``, then this should be a dictionary of the form ``{nuclide: contribution}``, where 
              contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index.
            * If ``nested`` is ``True``, then this should be a dictionary of the form ``{nuclide: {reaction: contribution}}``,
              where contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index 
              through the given reaction.
        nested
            Wether the contributions are on a reaction-wise basis or not."""
        
        self._nested = nested
        self._fig, self._axs = plt.subplots(figsize=(10, 10))
        self._textprops = dict(va="center", rotation_mode = 'anchor', fontsize=12)
        if nested:
            self._nested_pie_chart(contributions)
        else:
            self._pie_chart(contributions)

        self._style()

    def _add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.ax, sharey=self.ax)

    def _get_plot(self) -> Tuple[Figure, Axes]:
        return self._fig, self._axs

    def _nested_pie_chart(self, contributions: Dict[str, Dict[str, ufloat]]):
        """Create a pie chart of the contributions to the integral index on a nuclide-reaction-wise basis.

        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: {reaction: contribution}}``, where contribution is a ``ufloat`` object
            representing the contribution of the nuclide to the integral index through the given reaction."""
        
        def blend_colors(color1, color2, alpha):
            return np.array( [ alpha * c1 + (1 - alpha) * c2 for c1, c2 in zip(color1, color2 ) ] )
        
        # Create a nested ring chart
        nuclide_colors = plt.get_cmap('rainbow')(np.linspace(0, 1, len(contributions.keys())))
        nuclide_colors = [ blend_colors(color, [1, 1, 1, 1], 0.8) for color in nuclide_colors ]

        nuclide_totals = { nuclide: sum(contribution.n for contribution in contributions[nuclide].values()) \
                        for nuclide in contributions }
        
        # Sort nuclides by absolute magnitude of their total contributions
        nuclide_labels, nuclide_totals_sorted = zip(*sorted(nuclide_totals.items(), key=lambda x: abs(x[1]), reverse=True))

        # Now, deal with negative values

        nuclides_with_opposite_sign_contributions = []
        for nuclide, contribution in contributions.items():
            contribution_values = [contribution[reaction].n for reaction in contribution]
            if not (all(v >= 0 for v in contribution_values) or all(v <= 0 for v in contribution_values)):
                nuclides_with_opposite_sign_contributions.append(nuclide)
            
        # For nuclides with opposite sign contributions, we distinguish the positive and negative contributions
        # by coloring some of the inner ring a lighter color to indicate the negative contributions in the outer ring
        wedge_widths = list(nuclide_totals_sorted)
        inner_wedge_hatches = [None] * len(wedge_widths)

        if len(nuclides_with_opposite_sign_contributions) > 0:
            for nuclide in nuclides_with_opposite_sign_contributions:
                # First, determine the fraction of the contributions that are opposite (in sign) to the total
                total_sign = np.sign(nuclide_totals[nuclide])
                
                # Now, we want to plot the "lost" wedge width in white, i.e. the width lost from cancellations between the
                # positive and negative contributions. This will be colored a lighter color. The absolute sum of the
                # contributions represents the wedge width if there were no cancellations, so the total wedge width
                # minus the absolute sum of the contributions is "lost" wedge width.

                absolute_sum_of_contributions = sum(np.abs(contribution.n) for contribution in contributions[nuclide].values())
                
                # NOTE the sign function is needed to handle the case when the nuclide total is negative
                lost_wedge_width = absolute_sum_of_contributions - total_sign * nuclide_totals[nuclide]

                # Now, insert the lost wedge width into the wedge widths list right after the nuclide
                nuclide_index = list(nuclide_labels).index(nuclide)
                wedge_widths.insert(nuclide_index + 1, lost_wedge_width)
                nuclide_labels = list(nuclide_labels)
                nuclide_labels.insert(nuclide_index + 1, '')
                
                # The color of the lost wedge width will be a blend of the nuclide color and white
                white_color = np.array([1, 1, 1, 1])
                opacity = 1.0
                blended_color = blend_colors(white_color, nuclide_colors[nuclide_index], opacity)
                nuclide_colors = np.insert(nuclide_colors, nuclide_index + 1, blended_color, axis=0)
                
                # Add hatches to the negative total sum wedge
                if nuclide_totals[nuclide] < 0:
                    inner_wedge_hatches[nuclide_index] = '//'

        # Now make everything positive for the pie chart
        wedge_widths = np.abs(wedge_widths)

        # Plot the inner ring for nuclide totals
        inner_ring, _ = self._axs.pie(wedge_widths, radius=0.7, labels=nuclide_labels, \
                                colors=nuclide_colors, labeldistance=0.4, textprops={'fontsize': 8}, \
                                    rotatelabels=True, wedgeprops=dict(width=0.4, edgecolor='w'))

        # Add hatches to the negative total sum wedges
        for wedge, hatch in zip(inner_ring, inner_wedge_hatches):
            if hatch:
                wedge.set_hatch(hatch)

        # Plot the outer ring for reaction-specific contributions
        outer_labels = []
        outer_colors = []
        outer_sizes = []
        outer_hatches = []
        nuclide_colors = plt.get_cmap('rainbow')(np.linspace(0, 1, len(contributions.keys())))
        nuclide_colors = [ blend_colors(color, [1, 1, 1, 1], 0.8) for color in nuclide_colors ]
        for i, nuclide in enumerate( [ nuclide_label for nuclide_label in nuclide_labels if nuclide_label != '' ] ):
            blended_color = blend_colors(nuclide_colors[i], [1, 1, 1, 1], 0.8)
            
            # Sort reactions by absolute magnitude of their contributions
            reactions = contributions[nuclide]
            sorted_reactions = sorted(reactions.items(), key=lambda x: abs(x[1].n), reverse=True)
            
            for j, (reaction, contribution) in enumerate(sorted_reactions):
                outer_labels.append(reaction)
                
                outer_colors.append(blended_color)
                outer_sizes.append(np.abs(contribution.n))
                
                if contribution.n < 0:
                    outer_hatches.append('//')
                else:
                    outer_hatches.append(None)

        outer_ring, _ = self._axs.pie(outer_sizes, radius=1, labels=outer_labels, labeldistance=0.7, colors=outer_colors, \
                textprops={'fontsize': 8}, startangle=inner_ring[0].theta1, counterclock=True, \
                    rotatelabels=True, wedgeprops=dict(width=0.4, edgecolor='w'))

        # Add hatches to the negative contribution wedges
        for wedge, hatch in zip(outer_ring, outer_hatches):
            if hatch:
                wedge.set_hatch(hatch)
        
    def _pie_chart(self, contributions: Dict[str, ufloat]):
        """Create a pie chart of the contributions to the integral index on a nuclide-wise basis.
        
        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: contribution}``, where contribution is a ``ufloat`` object representing the
            contribution of the nuclide to the integral index."""
        labels, values = zip(*sorted(contributions.items(), key=lambda x: abs(x[1].n), reverse=True))
        values = [values.n for values in values]
        
        # Determining hatching patterns: empty string for positive, cross-hatch for negative
        hatches = ['//' if contributions[key].n < 0 else '' for key in labels]

        # Creating the pie chart  
        wedges, _ = self._axs.pie(np.abs(values), labels=labels, startangle=90, rotatelabels =True, textprops=self._textprops)

        # Applying hatching patterns to the wedges
        for wedge, hatch in zip(wedges, hatches):
            wedge.set_hatch(hatch)

    def _style(self):
        if self._plot_redundant and self._nested:
            title_text = f'Contributions to {self._index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self._index_name}'
        self._axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self._axs.set_title(title_text, pad=80) # Pad set to avoid overlap with labels

class _InteractivePiePlotter(_Plotter):
    """A class for creating interactive pie charts of contributions to integral indices."""
    _index_name: str
    """The name of the integral index whose contributions are being plotted (e.g. ``'E'`` or ``'c_k'``)."""
    _plot_redundant: bool
    """Whether to plot redundant/irrelevant reactions in the pie chart."""
    def __init__(self, integral_index_name: str, plot_redundant: bool=False, **kwargs: dict):
        """Initializes a pie chart of the contributions to the given integral index.
        
        Parameters
        ----------
        integral_index_name
            The name of the integral index whose contributions are to be plotted.
        plot_redundant
            Wether to include redundant/irrelevant reactions in the plot. NOTE: this only applies to nested plots, and
            only affects the plot title; it is expected that the provided data is consistent with the flag.
        kwargs
            Additional keyword arguments to control the behavior of the interactive legend.

            - interactive_legend (bool)
                Wether to include an interactive legend in the plot. Default is ``True``.
        Notes
        -----
        * Redundant reactions are defined as those which are derived from other reactions, e.g. 'total' and 'capture' reactions
          in SCALE.
        * Irrelevant reactions are defined as those which are not directly cross sections (but rather other nuclear data 
          parameters), e.g. 'chi' in SCALE.
        * A flag for including/excluding redundant/irrelevant reactions was provided since, if the user is expecting the
          the contributions to add up nicely, then the redundant reactions should be excluded, and if only cross sections are
          being considered, then the irrelevant reactions should be excluded.
        """
        # Check if the user wants an interactive legend
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = True
        
        self._index_name = integral_index_name
        self._plot_redundant = plot_redundant

    def _create_plot(self, contributions: Union[Dict[str, ufloat], Dict[str, Dict[str, ufloat]]], nested: bool):
        """Creates an interactive pie chart of the given contributions to the integral index.
        
        Parameters
        ----------
        contributions
            * If ``nested`` is ``False``, then this should be a dictionary of the form ``{nuclide: contribution}``, where 
              contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index.
            * If ``nested`` is ``True``, then this should be a dictionary of the form ``{nuclide: {reaction: contribution}}``,
              where contribution is a ``ufloat`` object representing the contribution of the nuclide to the integral index 
              through the given reaction.
        nested
            Wether the contributions are on a reaction-wise basis or not."""
        self.fig = make_subplots()

        # Prepare data for the sunburst chart
        self.nested = nested
        if nested:
            df = self._create_nested_sunburst_data(contributions)
        else:
            df = self._create_sunburst_data(contributions)
        
        # Create a sunburst chart
        self.fig = px.sunburst(
            data_frame=df,
            names='labels',
            parents='parents',
            ids='ids',
            values='normalized_values',
            custom_data=['values', 'uncertainties']
        )

        # Update hovertemplate with correct syntax
        self.fig.update_traces(
            hovertemplate=(
                "<b>%{label}</b><br>"
                "Value: %{customdata[0]:1.4E} +/- %{customdata[1]:1.4E}"  # Corrected format specifiers
                "<extra></extra>"  # This hides the trace info
            )
        )

        # Now style the plot
        self._style()

        self.fig.update_layout(
            autosize=True,
            width=None,  # Removes fixed width
            height=None,  # Removes fixed height
            margin=dict(l=5, r=5, t=30, b=5)
        )

        if self.interactive_legend:
            self.fig = InteractivePieLegend(self.fig, df)


    
    def _add_to_subplot(self, fig, position):
        if self.interactive_legend:
            raise ValueError("Interactive legends are not supported when adding to a subplot")
        else:
            for trace in self.fig.data:
                fig.add_trace(trace, row=position[0], col=position[1])
            return fig

    def _get_plot(self) -> Union[plotly.graph_objs.Figure, InteractivePieLegend]:
        return self.fig

    def _create_sunburst_data(self, contributions: Dict[str, ufloat]) -> pd.DataFrame:
        """Create a pandas dataframe for a (not nested) sunburst chart of contributions.
        
        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: contribution}``, where contribution is a ``ufloat`` object representing the
            contribution of the nuclide to the integral index.
        
        Returns
        -------
            The dataframe for creating the sunburst chart."""
        data = {
            'labels': [], 
            'ids': [], 
            'parents': [], 
            'values': [], 
            'uncertainties': [],
            'normalized_values': [],
            'nuclide': []
        }

        abs_sum_of_nuclide_totals = sum( abs(contribution.n) for contribution in contributions.values())

        for nuclide, nuclide_total in contributions.items():
            # Caclulate the nuclide total, and the positive and negative contributions
            norm_nuclide_total = abs(nuclide_total) / abs_sum_of_nuclide_totals

            # Add the nuclide as a parent
            data['labels'].append(nuclide)
            data['ids'].append(nuclide)
            data['parents'].append('')
            data['values'].append(nuclide_total.n)
            data['uncertainties'].append(nuclide_total.s)
            data['normalized_values'].append(norm_nuclide_total.n)
            data['nuclide'].append(nuclide)

        return pd.DataFrame(data)

    def _create_nested_sunburst_data(self, contributions: Dict[str, Dict[str, ufloat]]) -> pd.DataFrame:
        """Create a pandas dataframe for a nested sunburst chart of contributions.
        
        Parameters
        ----------
        contributions
            A dictionary of the form ``{nuclide: {reaction: contribution}}``, where contribution is a ``ufloat`` object
            representing the contribution of the nuclide to the integral index through the given reaction.
        
        Returns
        -------
            The dataframe for creating the sunburst chart."""
        data = {
            'labels': [], 
            'ids': [], 
            'parents': [], 
            'values': [], 
            'uncertainties': [],
            'normalized_values': [],
            'nuclide': []
        }

        abs_sum_of_nuclide_totals = sum(sum(abs(contribution.n) for contribution in reactions.values()) \
                                    for reactions in contributions.values())

        for nuclide, reactions in contributions.items():
            # Caclulate the nuclide total, and the positive and negative contributions
            nuclide_total = sum(contribution for contribution in reactions.values())
            if abs_sum_of_nuclide_totals != 0:
                norm_nuclide_total = abs(nuclide_total) / abs_sum_of_nuclide_totals
            else:
                norm_nuclide_total = 0

            positive_contributions = { reaction: contribution for reaction, contribution in reactions.items() \
                                      if contribution.n >= 0 }
            negative_contributions = { reaction: contribution for reaction, contribution in reactions.items() \
                                      if contribution.n < 0 }
            positive_total = sum(contribution for contribution in positive_contributions.values())
            negative_total = sum(contribution for contribution in negative_contributions.values())

            # Add the nuclide as a parent
            data['labels'].append(nuclide)
            data['ids'].append(nuclide)
            data['parents'].append('')
            data['values'].append(nuclide_total.n)
            data['uncertainties'].append(nuclide_total.s)
            data['normalized_values'].append(norm_nuclide_total.n)
            data['nuclide'].append(nuclide)
    
            # --------------------------------------------------------
            # Add the positive and negative contributions as children
            # --------------------------------------------------------

            # Normalize the contributions by the absolute value of the nuclide total 
            absolute_sum = positive_total + abs(negative_total)
            if absolute_sum != 0:
                normalization_factor = abs(norm_nuclide_total) / absolute_sum
            else:
                normalization_factor = 0

            # Positive contributions
            if positive_total != 0:
                norm_positive_total = positive_total * normalization_factor
                data['labels'].append('Positive')
                data['ids'].append(f"{nuclide}-Positive")
                data['parents'].append(nuclide)
                data['values'].append(positive_total.n)
                data['uncertainties'].append(positive_total.s)
                data['normalized_values'].append( norm_positive_total.n )
                data['nuclide'].append(nuclide)
            else:
                norm_positive_total = 0

            # Negative contributions
            if negative_total != 0:
                norm_negative_total = abs(negative_total) * normalization_factor
                data['labels'].append('Negative')
                data['ids'].append(f"{nuclide}-Negative")
                data['parents'].append(nuclide)
                data['values'].append(negative_total.n)
                data['uncertainties'].append(negative_total.s)
                data['normalized_values'].append( norm_negative_total.n )
                data['nuclide'].append(nuclide)
            else:
                norm_negative_total = 0

            # -------------------------------
            # Add the reaction contributions
            # -------------------------------
            # NOTE: Plotly express apparently has issues dealing with small numbers, so unless the contribution is
            # multiplied by a sufficiently large scale factor, the data won't be displayed correctly
            scale_factor = 10000
            for reaction, contribution in positive_contributions.items():
                # Now normalize contributions so they sum to the "normalized_positive_total
                if positive_total != 0:
                    normalization_factor = norm_positive_total / positive_total
                else:
                    normalization_factor = 0
                norm_reaction_contribution = contribution.n * normalization_factor
                
                if contribution.n != 0:
                    data['labels'].append(reaction)
                    data['ids'].append(f"{nuclide}-{reaction}")
                    data['parents'].append(f"{nuclide}-Positive")
                    data['values'].append(contribution.n)
                    data['uncertainties'].append(contribution.s)
                    data['normalized_values'].append(scale_factor*norm_reaction_contribution.n)
                    data['nuclide'].append(nuclide)

            for reaction, contribution in negative_contributions.items():
                # Now normalize contributions so they sum to the "normalized_negative_total"
                normalization_factor = norm_negative_total / abs(negative_total)
                norm_reaction_contribution = abs(contribution.n) * normalization_factor

                if contribution.n != 0:
                    data['labels'].append(reaction)
                    data['ids'].append(f"{nuclide}-{reaction}")
                    data['parents'].append(f"{nuclide}-Negative")
                    data['values'].append(contribution.n)
                    data['uncertainties'].append(contribution.s)
                    data['normalized_values'].append(scale_factor*norm_reaction_contribution.n)
                    data['nuclide'].append(nuclide)


        return pd.DataFrame(data)

    def _style(self):
        if self._plot_redundant and self.nested:
            title_text = f'Contributions to {self._index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self._index_name}'
        self.fig.update_layout(title_text=title_text, title_x=0.5)  # 'title_x=0.5' centers the title


class InteractivePieLegend:
    """A class for creating an interactive legend for a sunburst chart."""
    fig: plotly.graph_objs.Figure
    """The sunburst chart for which the interactive legend is being created."""
    df: pd.DataFrame
    """The dataframe used to create the sunburst chart."""
    _app: Flask
    """The Flask webapp that will display the interactive legend."""
    def __init__(self, fig: plotly.graph_objs.Figure, df: pd.DataFrame):
        """Create a flask webapp that will display an interactive legend for the sunburst chart.
        
        Parameters
        ----------
        fig
            The sunburst chart for which the interactive legend is being created.
        df
            The dataframe used to create the sunburst chart."""
        self.fig = fig
        self.df = df
        self._app = Flask(__name__)

        @self._app.route('/shutdown', methods=['POST'])
        def shutdown():
            """Function to shutdown the server"""
            os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
            return 'Server shutting down...'

        @self._app.route('/')
        def show_sunburst() -> str:
            """Function to display the sunburst chart with an interactive legend
            
            Returns
            -------
                The HTML content for the sunburst chart with the interactive legend."""
            # Extract root nodes (nodes without parents)
            root_nodes = self.df[self.df['parents'] == '']

            # Generate a unique ID for the container
            container_id = f"container-{uuid.uuid4()}"

            # Generate legend HTML with a title
            legend_html = f'<div id="{container_id}-legend" style="border: 2px solid black; padding: 10px;"><h3 style="margin-top: 0; text-align: center;">Legend</h3>\n'
            for _, row in root_nodes.iterrows():
                legend_html += f'<div class="{container_id}-legend-item" style="cursor: pointer; margin-bottom: 5px;" data-target="{row["ids"]}">{row["ids"]}: {row["values"]:1.4E}</div>\n'
            legend_html += '</div>\n'

            # JavaScript for interactivity and shutdown
            script_html = f"""
            <script>
            window.addEventListener('beforeunload', (event) => {{
                navigator.sendBeacon('/shutdown');
            }});
            document.addEventListener('DOMContentLoaded', function () {{
                const legendItems = document.querySelectorAll('.{container_id}-legend-item');
                legendItems.forEach(item => {{
                    item.addEventListener('mouseenter', function() {{
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {{
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {{
                                path.style.opacity = 0.5; // Highlight by changing opacity
                            }}
                        }});
                    }});
                    item.addEventListener('mouseleave', function() {{
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {{
                            path.style.opacity = 1; // Reset opacity
                        }});
                    }});
                    item.addEventListener('click', function() {{
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {{
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {{
                                path.dispatchEvent(new MouseEvent('click', {{ 'view': window, 'bubbles': true, 'cancelable': true }}));
                            }}
                        }});
                    }});
                }});
                // Force Redraw/Reflow
                setTimeout(() => {{
                    window.dispatchEvent(new Event('resize'));
                }}, 100); // Delay may be adjusted based on actual rendering time
            }});
            </script>
            """

            # Save the chart with interactivity and layout adjustments
            fig_html = self.fig.to_html(full_html=False, include_plotlyjs='cdn')
            full_html = f"""
            <!DOCTYPE html>
            <html>
            <head>
            <title>Interactive Sunburst Chart</title>
            <style>
                #{container_id} {{
                    display: flex;
                    flex-direction: row; /* Align children horizontally */
                    height: 100%;
                    width: 100%; /* Ensure the container takes full width */
                    margin: 0;
                    font-family: Arial, sans-serif;
                }}
                #{container_id} > div {{
                    display: flex;
                    justify-content: space-between; /* Space out the chart and legend */
                    align-items: flex-start; /* Align items at the start of the cross axis */
                    width: 100%;
                    overflow: hidden; /* Hide overflow to prevent breaking layout */
                }}
                #{container_id}-chart {{
                    flex: 1 1 70%; /* Allow chart to grow and shrink but base at 70% width */
                    padding: 10px;
                }}
                #{container_id}-legend {{
                    flex: 0 1 30%; /* Start with 30% width but allow shrinking */
                    padding: 5px;
                    max-height: calc(100vh - 20px); /* Limit height to viewport height minus some margin */
                    overflow: auto; /* Scroll internally if content overflows */
                }}
            </style>
            </head>
            <body>
            <div id="{container_id}">
                <div>
                    <div id="{container_id}-chart">{fig_html}</div>
                    <div id="{container_id}-legend">{legend_html}</div>
                </div>
            </div>
            {script_html}
            </body>
            </html>
            """

            return render_template_string(full_html)
        
    def _open_browser(self, port: int):
        """Open the browser to display the interactive sunburst chart

        Parameters
        ----------
        port
            The port at which the Flask server is running.
        """
        print(f"Now running at http://localhost:{port}/")
        webbrowser.open(f"http://localhost:{port}/")

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
        if silent:
            sys.stdout = log
        sys.stderr = log

        port = _find_free_port()
        if not config.generating_docs:
            if open_browser:
                threading.Timer(1, self._open_browser(port)).start()

            
            self._app.run(host='localhost', port=port)

    def serve(self):
        """Start the Flask server to display the interactive sunburst chart without a browser tab."""
        port = _find_free_port()
        log = open(os.devnull, 'w')
        # sys.stdout = log
        # sys.stderr = log

        # Run the Flask application in a separate thread
        thread = threading.Thread(target=lambda: self._app.run(host='localhost', port=port))
        print(f"Now running at http://localhost:{port}/")
        thread.daemon = True  # This ensures thread exits when main program exits
        thread.start()

    def write_html(self, filename: Optional[Union[str, Path]]=None) -> Union[None, str]:
        """Write the HTML content of the interactive sunburst chart to a file.
        
        Parameters
        ----------
        filename
            The path of the file to which the HTML content will be written. If not provided, the content will be returned.
            
        Returns
        -------
            If ``filename`` is not provided, the HTML content of the interactive sunburst chart."""
        with self._app.test_client() as client:
            response = client.get('/')
            html_content = response.data.decode('utf-8')
            if filename is None:
                return html_content
            else:
                with open(filename, 'w') as f:
                    f.write(html_content)

    def save_state(self, filename: Optional[Union[str, Path]]=None
                   ) -> Union[None, Dict[str, Union[plotly.graph_objs.Figure, pd.DataFrame]]]:
        """Save the figure (with interactive legend) as a pickle that can be deserialized for plotting later with full
        interactivity.
        
        Parameters
        ----------
        filename
            Path to the file to store the pickle of the plot object (should have a ``.pkl`` extension).

        Returns
        -------
            If ``filename`` is not provided, a dictionary containing the figure and dataframe used to create the plot. This
            data can be used to researialize the plot later.
            
        Notes
        -----
        This data can be loaded later using the :meth:`load_state` method."""
        state = {
            'fig': self.fig,
            'df': self.df,
        }
        if filename is None:
            return state
        else:
            with open(filename, 'wb') as f:
                pickle.dump(state, f)

    @classmethod
    def load_state(cls, filename: Optional[Union[str, Path]]=None, data_dict: Optional[dict]=None) -> InteractivePieLegend:
        """Researilize an :class:`InteractivePieLegend` plot from a pickle file or a data dictionary (produced by 
        :meth:`save_state`).
        
        Parameters
        ----------
        filename
            Path to the file containing the pickle of the plot object.
        data_dict
            A dictionary containing the figure and dataframe used to create the plot.
        
        Returns
        -------
            A new instance of the :class:`InteractivePieLegend` class with the same state as the original instance whose data
            was provided.
        """
        if filename is None and data_dict is None:
            raise ValueError("Either a filename or a data dictionary must be provided")
        if filename is not None:
            with open(filename, 'rb') as f:
                state = pickle.load(f)
        else:
            state = data_dict
        fig = state['fig']
        df = state['df']
        instance = cls(fig, df)
        return instance
    
    def to_image(self, filename: Union[str, Path]) -> None:
        """Write the initial state of the interactive plot to an image file. This function saves the plot as an 
        image by using selenium webdriver.
        
        Parameters
        ----------
        filename
            Name of the image file to save the matrix plot to. The file extension should be ``'.png'``, ``.jpg``, etc."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            html_filename = f.name + '.html'
        self.write_html(html_filename)
        _capture_html_as_image(html_filename, filename)
        os.remove(html_filename)