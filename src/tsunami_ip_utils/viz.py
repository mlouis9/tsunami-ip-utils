from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import pickle

from plotly.subplots import make_subplots
import plotly.express as px
from plotly import graph_objects as go
import pandas as pd
import scipy.stats as stats
from uncertainties import umath, unumpy, ufloat
from tsunami_ip_utils.utils import isotope_reaction_list_to_nested_dict
from tsunami_ip_utils.integral_indices import add_missing_reactions_and_nuclides
from tsunami_ip_utils.utils import filter_redundant_reactions

# Imports for the interactive legend
import webbrowser
import os, sys, signal
from flask import Flask, render_template_string
import threading
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import socket
import uuid
from pathlib import Path

plt.rcParams['hatch.linewidth'] = 0.5


PLOTTING_PORT = 8050

class Plotter(ABC):
    @abstractmethod
    def create_plot(self, data, nested):
        pass

    @abstractmethod
    def add_to_subplot(self, fig, position):
        pass

    @abstractmethod
    def get_plot(self):
        pass

    @abstractmethod
    def style(self):
        pass


class BarPlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested):
        self.nested = nested
        self.fig, self.axs = plt.subplots()
        if nested:
            self.nested_barchart(contributions)
        else:
            self.barchart(contributions)

        self.style()

    def get_plot(self):
        return self.fig, self.axs
        
    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
        
    def barchart(self, contributions):
        self.axs.bar(contributions.keys(), [contribution.n for contribution in contributions.values()],
            yerr=[contribution.s for contribution in contributions.values()], capsize=5, error_kw={'elinewidth': 0.5})

    def nested_barchart(self, contributions):
        # Colors for each reaction type
        num_reactions = len(next(iter(contributions.values())))
        cmap = plt.get_cmap('Set1')
        colors = cmap(np.linspace(0, 1, num_reactions))

        # Variables to hold the bar positions and labels
        indices = range(len(contributions))
        labels = list(contributions.keys())

        # Bottom offset for each stack
        bottoms_pos = [0] * len(contributions)
        bottoms_neg = [0] * len(contributions)

        color_index = 0
        for reaction in next(iter(contributions.values())).keys():
            values = [contributions[nuclide][reaction].n for nuclide in contributions]
            errs = [contributions[nuclide][reaction].s for nuclide in contributions]
            # Stacking positive values
            pos_values = [max(0, v) for v in values]
            neg_values = [min(0, v) for v in values]
            self.axs.bar(indices, pos_values, label=reaction, bottom=bottoms_pos, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            self.axs.bar(indices, neg_values, bottom=bottoms_neg, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            # Update the bottom positions
            bottoms_pos = [bottoms_pos[i] + pos_values[i] for i in range(len(bottoms_pos))]
            bottoms_neg = [bottoms_neg[i] + neg_values[i] for i in range(len(bottoms_neg))]
            color_index += 1

        # Adding 'effective' box with dashed border
        total_values = [sum(contributions[label][r].n for r in contributions[label]) for label in labels]
        for idx, val in zip(indices, total_values):
            self.axs.bar(idx, abs(val), bottom=0 if val > 0 else val, color='none', edgecolor='black', hatch='///', linewidth=0.5)

        self.axs.set_xticks(indices)
        self.axs.set_xticklabels(labels)
        self.axs.legend()

    def style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.set_ylabel(f"Contribution to {self.index_name}")
        self.axs.set_xlabel("Isotope")
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)


class PiePlotter(Plotter):
    def __init__(self, integral_index_name, plot_redudant=False, **kwargs):
        self.index_name = integral_index_name
        self.plot_redundant = plot_redudant
    
    def create_plot(self, contributions, nested):
        self.nested = nested
        self.fig, self.axs = plt.subplots()
        if nested:
            self.nested_pie_chart(contributions)
        else:
            self.pie_chart(contributions)

        self.style()

    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.ax, sharey=self.ax)

    def get_plot(self):
        return self.fig, self.axs

    def nested_pie_chart(self, contributions):
        # Create a nested ring chart
        num_reactions = len(next(iter(contributions.values())))
        nuclide_colors = plt.get_cmap('rainbow')(np.linspace(0, 1, len(contributions.keys())))
        nuclide_totals = { nuclide: sum(contribution.n for contribution in contributions[nuclide].values()) \
                        for nuclide in contributions }
        nuclide_labels = list(nuclide_totals.keys())

        # Now, deal with negative values

        nuclides_with_opposite_sign_contributions = []
        for nuclide, contribution in contributions.items():
            contribution_values = [contribution[reaction].n for reaction in contribution]
            if not (all(v >= 0 for v in contribution_values) or all(v <= 0 for v in contribution_values)):
                nuclides_with_opposite_sign_contributions.append(nuclide)
            
        # For nuclides with opposite sign contributions, we distinguish the positive and negative contributions
        # by coloring some of the inner ring a lighter color to indicate the negative contributions in the outer ring
        wedge_widths = list(nuclide_totals.values())
        inner_wedge_hatches = [None] * len(wedge_widths)

        def blend_colors(color1, color2, alpha):
            return np.array( [ alpha * c1 + (1 - alpha) * c2 for c1, c2 in zip(color1, color2 ) ] )

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
                nuclide_index = list(nuclide_totals.keys()).index(nuclide)
                wedge_widths.insert(nuclide_index + 1, lost_wedge_width)
                nuclide_labels.insert(nuclide_index + 1, '')
                
                # The color of the lost wedge width will be a blend of the nuclide color and white
                white_color = np.array([1, 1, 1, 1])
                opacity = 0.8
                blended_color = blend_colors(white_color, nuclide_colors[nuclide_index], opacity)
                nuclide_colors = np.insert(nuclide_colors, nuclide_index + 1, blended_color, axis=0)
                
                # Add hatches to the negative total sum wedge
                if nuclide_totals[nuclide] < 0:
                    inner_wedge_hatches[nuclide_index] = '//'

        # Now make everything positive for the pie chart
        wedge_widths = np.abs(wedge_widths)

        # Plot the inner ring for nuclide totals
        inner_ring, _ = self.axs.pie(wedge_widths, radius=0.7, labels=nuclide_labels, \
                                colors=nuclide_colors, labeldistance=0.6, textprops={'fontsize': 8}, \
                                    wedgeprops=dict(width=0.3, edgecolor='w'))

        # Add hatches to the negative total sum wedges
        for wedge, hatch in zip(inner_ring, inner_wedge_hatches):
            if hatch:
                wedge.set_hatch(hatch)

        # Get colors for reactions from the "rainbow" colormap
        reaction_colors = plt.get_cmap('Set1')(np.linspace(0, 1, num_reactions))

        # Plot the outer ring for reaction-specific contributions
        outer_labels = []
        outer_colors = []
        outer_sizes = []
        outer_hatches = []
        for i, (nuclide, reactions) in enumerate(contributions.items()):
            for j, (reaction, contribution) in enumerate(list(reactions.items())):
                outer_labels.append(f"{nuclide} - {reaction}")
                
                outer_colors.append(reaction_colors[j])
                outer_sizes.append(np.abs(contribution.n))
                
                if contribution.n < 0:
                    outer_hatches.append('//')
                else:
                    outer_hatches.append(None)

        outer_ring, _ = self.axs.pie(outer_sizes, radius=1, labels=outer_labels, labeldistance=0.9, colors=outer_colors, \
                textprops={'fontsize': 6}, startangle=inner_ring[0].theta1, counterclock=True, \
                    wedgeprops=dict(width=0.3, edgecolor='w'))

        # Add hatches to the negative contribution wedges
        for wedge, hatch in zip(outer_ring, outer_hatches):
            if hatch:
                wedge.set_hatch(hatch)
        
    def pie_chart(self, contributions):
        labels = list(contributions.keys())
        values = [abs(contributions[key].n) for key in labels]

        # Determining hatching patterns: empty string for positive, cross-hatch for negative
        hatches = ['//' if contributions[key].n < 0 else '' for key in labels]

        # Creating the pie chart
        wedges, _ = self.axs.pie(values, labels=labels, startangle=90)

        # Applying hatching patterns to the wedges
        for wedge, hatch in zip(wedges, hatches):
            wedge.set_hatch(hatch)

    def style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)


class InteractivePieLegend:
    def __init__(self, fig, df):
        """Return a flask webapp that will display an interactive legend for the sunburst chart"""
        self.fig = fig
        self.df = df
        self.port = find_free_port()
        self.app = Flask(__name__)

        @self.app.route('/shutdown', methods=['POST'])
        def shutdown():
            os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
            return 'Server shutting down...'

        @self.app.route('/')
        def show_sunburst():
            # Extract root nodes (nodes without parents)
            root_nodes = self.df[self.df['parents'] == '']

            # Generate a unique ID for the container
            container_id = f"container-{uuid.uuid4()}"

            # Generate legend HTML with a title
            legend_html = f'<div id="{container_id}-legend" style="margin-left: 20px; border: 2px solid black; padding: 10px;"><h3 style="margin-top: 0; text-align: center;">Legend</h3>\n'
            for _, row in root_nodes.iterrows():
                legend_html += f'<div class="{container_id}-legend-item" style="cursor: pointer; margin-bottom: 5px;" data-target="{row["ids"]}">{row["ids"]}: {row["values"]:1.4E}</div>\n'
            legend_height = len(root_nodes) * 30  # Calculate a dynamic height based on number of items
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
                    margin-left: 20px;
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
        
    def open_browser(self):
        print(f"Now running at http://localhost:{PLOTTING_PORT}/")
        webbrowser.open(f"http://localhost:{PLOTTING_PORT}/")
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

        if open_browser:
            threading.Timer(1, self.open_browser).start()
        self.app.run(host='localhost', port=self.port)

    def serve(self):
        """Start the Flask server to display the interactive sunburst chart"""
        log = open(os.devnull, 'w')
        # sys.stdout = log
        # sys.stderr = log

        # Run the Flask application in a separate thread
        thread = threading.Thread(target=lambda: self.app.run(host='localhost', port=self.port))
        print(f"Now running at http://localhost:{self.port}/")
        thread.daemon = True  # This ensures thread exits when main program exits
        thread.start()

    def write_html(self, filename=None):
        with self.app.test_client() as client:
            response = client.get('/')
            html_content = response.data.decode('utf-8')
            if filename is None:
                return html_content
            else:
                with open(filename, 'w') as f:
                    f.write(html_content)


class InteractivePiePlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        # Check if the user wants an interactive legend
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = True
        
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested=True):
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
        self.style()

        self.fig.update_layout(
            autosize=True,
            width=None,  # Removes fixed width
            height=None,  # Removes fixed height
            margin=dict(l=5, r=5, t=30, b=5)
        )

        if self.interactive_legend:
            self.fig = InteractivePieLegend(self.fig, df)


    
    def add_to_subplot(self, fig, position):
        if self.interactive_legend:
            raise ValueError("Interactive legends are not supported when adding to a subplot")
        else:
            for trace in self.fig.data:
                fig.add_trace(trace, row=position[0], col=position[1])
            return fig

    def get_plot(self):
        return self.fig

    def _create_sunburst_data(self, contributions):
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

    def _create_nested_sunburst_data(self, contributions):
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

    def style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.fig.update_layout(title_text=title_text, title_x=0.5)  # 'title_x=0.5' centers the title


def determine_plot_type(contributions, plot_redundant_reactions):
    """Determines whether the contributions are nuclide-wise or nuclide-reaction-wise and whether to plot redundant
    reactions or not
    
    Parameters
    ----------
    - contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair
    - plot_redundant_reactions: bool, whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions
        
    Returns
    -------
    - contributions: dict, contributions to the similarity parameter keyed by isotope then by reaction type"""
    if 'reaction_type' in contributions[0]: # Nuclide-reaction-wise contributions
        nested_plot = True # Nested plot by nuclide then by reaction type

        # Create a dictionary of contributions keyed by isotope then by reaction type
        contributions = isotope_reaction_list_to_nested_dict(contributions, 'contribution')

        # If viewing nuclide-reaction wise contributions, it's important (at least for the visualizations in this function)
        # that if viewing the true contributions to the nuclide total, that redundant interactions (e.g. capture and fission
        # + (n, g)) and irrelevant interactions (e.g. chi and nubar) are not plotted.

        if not plot_redundant_reactions:
            # Remove redundant interactions
            contributions = filter_redundant_reactions(contributions)
    else: # Nuclide-wise contributions
        nested_plot = False
        contributions = { contribution['isotope']: contribution['contribution'] for contribution in contributions }

    return contributions, nested_plot

def contribution_plot(contributions, plot_type='bar', integral_index_name='E', plot_redundant_reactions=True, **kwargs):
    """Plots the contributions to an arbitrary similarity parameter for a single experiment application pair
    
    Parameters
    ----------
    - contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair
    - plot_type: str, type of plot to create. Default is 'bar' which creates a bar plot. Other option is 'pie' which creates
        a pie chart
    - integral_index_name: str, name of the integral index being plotted. Default is 'E'
    - plot_redundant_reactions: bool, whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is True
    - kwargs: additional keyword arguments to pass to the plotting function"""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    contributions, nested_plot = determine_plot_type(contributions, plot_redundant_reactions)

    plotters = {
        'bar': BarPlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'pie': PiePlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'interactive_pie': InteractivePiePlotter(integral_index_name, plot_redundant_reactions, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError("Unsupported plot type")

    # Create the plot and style it
    plotter.create_plot(contributions, nested_plot)

    return plotter.get_plot()


def manual_pearson(x, y):
    """Calculates the Pearson correlation coefficient between two sets of values x and y"""
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_diff = x - x_mean
    y_diff = y - y_mean
    numerator = np.sum(x_diff * y_diff)
    denominator = umath.sqrt(np.sum(x_diff**2) * np.sum(y_diff**2))
    return numerator / denominator


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
        self.pearson = stats.pearsonr(x, y).statistic
        self.spearman = stats.spearmanr(x, y).statistic
        self.slope = self.regression.slope
        self.intercept = self.regression.intercept

        # Now create teh summary statistics text for figure annotation
        pearson_text = f"Pearson: {self.pearson:1.6f}"
        spearman_text = f"Spearman: {self.spearman:1.6f}"
        self.summary_stats_text = f"{pearson_text}\n{spearman_text}"


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

        # Add correlation statistics to the plot
        self.fig.add_annotation(
            x=0.05, xref="paper", 
            y=0.95, yref="paper",
            text=self.summary_stats_text,
            showarrow=False, 
            font=dict(size=12),
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

            # Create a new InteractiveScatterPlotter instance with the updated data
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
        # Function to open the browser
        def open_browser():
            if not os.environ.get("WERKZEUG_RUN_MAIN"):
                print(f"Now running at http://localhost:{PLOTTING_PORT}/")
                webbrowser.open(f"http://localhost:{PLOTTING_PORT}/")

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

        self.app.run_server(debug=False, host='localhost', port=PLOTTING_PORT)

    def save_state(self, filename):
        state = {
            'fig': self.fig.to_dict(),
            'df': self.df.to_dict(),
            'excluded_isotopes': self.excluded_isotopes,
            'index_name': self.index_name,
            'nested': self.interactive_scatter_plot.nested
        }
        with open(filename, 'wb') as f:
            pickle.dump(state, f)

    @classmethod
    def load_state(cls, filename):
        with open(filename, 'rb') as f:
            state = pickle.load(f)

        # Recreate the InteractiveScatterPlotter instance from the saved state
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

def load_interactive_scatter_plot(filename):
    """Loads an interactive scatter plot from a saved state file. This function is purely for convenience and is a
    wrapper around the InteractiveScatterLegend.load_state method"""
    return InteractiveScatterLegend.load_state(filename)

def correlation_plot(application_contributions, experiment_contributions, plot_type='scatter', integral_index_name='E', \
                     plot_redundant_reactions=True, **kwargs):
    """Creates a correlation plot for a given application-experiment pair for which the contributions to the similarity
    parameter are given.
    
    Parameters
    ----------
    - application_contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter
        for the application
    - experiment_contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter
        for the experiment
    - plot_type: str, type of plot to create. Default is 'scatter' which creates a matplotlib scatter plot. Other options
        are interactive_scatter, which creates a Plotly scatter plot.
    - integral_index_name: str, name of the integral index being plotted. Default is 'E'
    - plot_redundant_reactions: bool, whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is True
        
    Returns
    -------
    - fig: matplotlib.figure.Figure, the figure containing the correlation plot"""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    application_contributions, application_nested = determine_plot_type(application_contributions, plot_redundant_reactions)
    experiment_contributions, experiment_nested = determine_plot_type(experiment_contributions, plot_redundant_reactions)

    if application_nested != experiment_nested:
        raise ValueError("Application and experiment contributions must have the same nested structure")
    else:
        nested = application_nested # They are be the same, so arbitrarily choose one

    # Get the list of isotopes for which contributions are available
    isotopes = list(set(application_contributions.keys()).union(experiment_contributions.keys()))

    all_reactions = add_missing_reactions_and_nuclides(application_contributions, experiment_contributions, isotopes, mode='contribution')

    # Now convert the contributions for the application and experiment into a list of x, y pairs for plotting
    contribution_pairs = []
    if nested:
        for isotope in isotopes:
            for reaction in all_reactions:
                contribution_pairs.append((application_contributions[isotope][reaction], \
                                           experiment_contributions[isotope][reaction]))
    else:
        for isotope in isotopes:
            contribution_pairs.append((application_contributions[isotope], experiment_contributions[isotope]))

    plotters = {
        'scatter': ScatterPlotter(integral_index_name, plot_redundant_reactions, nested, **kwargs),
        'interactive_scatter': InteractiveScatterPlotter(integral_index_name, plot_redundant_reactions, nested, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError("Unsupported plot type")

    # Create the plot and style it
    plotter.create_plot(contribution_pairs, isotopes, all_reactions)

    return plotter.fig

def perturbation_plot(points):
    """Plots the perturbation points for a given application-experiment pair for which the perturbation points have already
    been calculated
    
    Parameters
    ----------
    - points: list of tuple of ufloat, list of tuples containing the perturbation points for the application-experiment pair"""
    
    # Extracting data
    plotter = InteractivePerturbationScatterPlotter()
    plotter.create_plot(points)

    return plotter.get_plot(), plotter.pearson

def find_free_port():
    """Finds a free port on localhost for running a Flask server."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))  # Let the OS pick an available port
    port = s.getsockname()[1]
    s.close()
    return port

def matrix_plot(plot_type: str, plot_objects_array: np.ndarray):
    """Creates a Dash app to display a matrix of plots from a numpy object array of figure objects."""
    current_directory = Path(__file__).parent
    external_stylesheets = [ str( current_directory / 'css' / 'matrix_plot.css' )]  # Modify this path to your actual CSS file
    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    num_rows = plot_objects_array.shape[0]
    num_cols = plot_objects_array.shape[1]

    # Create column headers
    column_headers = [html.Div(f'Application {i+1}', style={'flex': '1', 'minWidth': '700px', 'textAlign': 'center'}) for i in range(num_cols)]
    header_row = html.Div([html.Div('', style={'flex': 'none', 'width': '100px'})] + column_headers, 
                          style={'display': 'flex', 'marginBottom': '10px'})

    # Create rows of plots
    rows = [header_row]
    for i in range(num_rows):
        # Create a row of plots
        row = [html.Div(f'Experiment {i+1}', style={'flex': 'none', 'width': '100px', 'textAlign': 'center', 'marginRight': '10px'})]

        for j in range(num_cols):
            plot_object = plot_objects_array[i, j]
            if plot_object is not None:
                # Style parameters for plotly graphs and iframes
                graph_style = {'flex': '1', 'minWidth': '700px', 'height': '400px'}

                if isinstance(plot_object, InteractiveScatterLegend):
                    # Handle the special case of interactive legend plots
                    interactive_legend_app = plot_object
                    graph_id = f"interactive-scatter-{i}-{j}"
                    
                    # Add the plot to the row using dcc.Graph
                    row.append(dcc.Graph(id=graph_id, figure=interactive_legend_app.fig, style=graph_style))
                    
                    # Create a closure to capture the correct interactive_legend_app instance
                    def create_update_figure_callback(app_instance):
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
                    
                    # Create the callback using the closure
                    create_update_figure_callback(interactive_legend_app)
                elif isinstance(plot_object, InteractivePieLegend):
                    # Embed the Flask app using iframe with srcdoc
                    with plot_object.app.test_client() as client:
                        response = client.get('/')
                        html_content = response.data.decode('utf-8')
                        iframe_html = html.Iframe(srcDoc=html_content, style=graph_style)
                    row.append(iframe_html)
                else:
                    row.append(dcc.Graph(figure=plot_object, style=graph_style))
            else:
                row.append(html.Div('Plot not available', style=graph_style))

        # Append the row to the rows list
        rows.append(html.Div(row, style={'display': 'flex', 'marginBottom': '10px'}))

    # Generate the layout for the app
    app.layout = html.Div([
        html.H1("Matrix of Plots", style={'textAlign': 'center'}),
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

    return app