from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import numpy as np

from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
import scipy.stats as stats
from uncertainties import umath, unumpy, ufloat

# Imports for the interactive legend
import webbrowser
import os, sys, signal
from flask import Flask, render_template_string
import threading

from tsunami_ip_utils.utils import isotope_reaction_list_to_nested_dict

plt.rcParams['hatch.linewidth'] = 0.5


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


class InteractiveLegend:
    def __init__(self, fig, df):
        """Return a flask webapp that will display an interactive legend for the sunburst chart"""
        self.fig = fig

        self.app = Flask(__name__)

        @self.app.route('/shutdown', methods=['POST'])
        def shutdown():
            os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
            return 'Server shutting down...'

        @self.app.route('/')
        def show_sunburst():
            # Extract root nodes (nodes without parents)
            root_nodes = df[df['parents'] == '']

            # Generate legend HTML with a title
            legend_html = '<div id="legend" style="margin-left: 20px; border: 2px solid black; padding: 10px;"><h3 style="margin-top: 0; text-align: center;">Legend</h3>\n'
            for _, row in root_nodes.iterrows():
                legend_html += f'    <div id="legend-item" class="legend-item" style="cursor: pointer; margin-bottom: 5px;" data-target="{row["ids"]}">{row["ids"]}: {row["values"]:1.4E}</div>\n'
            legend_html += '</div>\n'

            # JavaScript for interactivity and shutdown
            script_html = """
            <script>
            window.addEventListener('beforeunload', (event) => {
                navigator.sendBeacon('/shutdown');
            });
            document.addEventListener('DOMContentLoaded', function () {
                const legendItems = document.querySelectorAll('.legend-item');
                legendItems.forEach(item => {
                    item.addEventListener('mouseenter', function() {
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {
                                path.style.opacity = 0.5; // Highlight by changing opacity
                            }
                        });
                    });
                    item.addEventListener('mouseleave', function() {
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            path.style.opacity = 1; // Reset opacity
                        });
                    });
                    item.addEventListener('click', function() {
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {
                                path.dispatchEvent(new MouseEvent('click', { 'view': window, 'bubbles': true, 'cancelable': true }));
                            }
                        });
                    });
                });
            });
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
                body, html {{
                    height: 100%;
                    margin: 0;
                    font-family: Arial, sans-serif;
                }}
                #container {{
                    display: flex;
                    justify-content: center;
                    align-items: center;
                    width: 100%;
                }}
                div:not(#legend, 
                .modebar-container) {{
                    height: 100%;
                    min-width: 0;    /* Prevent flex item from overflowing its container */
                }}
                #legend {{
                    flex: 0 1 auto;  /* Do not grow, allow shrink */
                    margin-left: 20px;
                    padding: 10px;
                    border: 2px solid black;
                    max-height: 100%; /* Make sure it does not overflow vertically */
                    overflow: auto;   /* Enable scrolling if content is too large */
                }}
            </style>
            </head>
            <body>
            <div id="container">
                <div id="chart">{fig_html}</div>
                {legend_html}
            </div>
            {script_html}
            </body>
            </html>
            """

            return render_template_string(full_html)
        
    def open_browser(self):
        print("Now running at http://localhost:5000/")
        webbrowser.open("http://localhost:5000/")
        pass

    def show(self):
        # Suppress Flask's startup and runtime messages by redirecting them to dev null
        log = open(os.devnull, 'w')
        # sys.stdout = log
        sys.stderr = log

        threading.Timer(1, self.open_browser).start()
        self.app.run()

    def write_html(self, filename):
        with self.app.test_client() as client:
            response = client.get('/')
            html_content = response.data.decode('utf-8')
            with open(filename, 'w') as f:
                f.write(html_content)


class InteractivePiePlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.fig = make_subplots()

        # Check if the user wants an interactive legend
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = True
        
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested=True):
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
            self.fig = InteractiveLegend(self.fig, df)


    
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
            redundant_interactions = ['chi', 'capture', 'nubar', 'total']
            contributions = { isotope: { reaction: contributions[isotope][reaction] for reaction in contributions[isotope] \
                                if reaction not in redundant_interactions } for isotope in contributions }
    else: # Nuclide-wise contributions
        nested_plot = False
        contributions = { contribution['isotope']: contribution['contribution'] for contribution in contributions }

    return contributions, nested_plot

def contribution_plot(contributions, plot_type='bar', integral_index_name='E', plot_redundant_reactions=False, **kwargs):
    """Plots the contributions to an arbitrary similarity parameter for a single experiment application pair
    
    Parameters
    ----------
    - contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair
    - plot_type: str, type of plot to create. Default is 'bar' which creates a bar plot. Other option is 'pie' which creates
        a pie chart
    - integral_index_name: str, name of the integral index being plotted. Default is 'E'
    - plot_redundant_reactions: bool, whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is False
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

class ScatterPlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contribution_pairs, nested):
        self.nested = nested
        self.fig, self.axs = plt.subplots()

        application_points        = [ contribution[0].n for contribution in contribution_pairs ]
        application_uncertainties = [ contribution[0].s for contribution in contribution_pairs ]
        experiment_points         = [ contribution[1].n for contribution in contribution_pairs ]
        experiment_uncertainties  = [ contribution[1].s for contribution in contribution_pairs ]

        application               = unumpy.uarray(application_points, application_uncertainties)
        experiment                = unumpy.uarray(experiment_points, experiment_uncertainties)
        application               = np.array([ufloat(v.n, v.s) for v in application])
        experiment                = np.array([ufloat(v.n, v.s) for v in experiment])

        self.fig = plt.errorbar(application_points, experiment_points, xerr=application_uncertainties, \
                               yerr=experiment_uncertainties, fmt='.', capsize=5)
        
        # Linear regression
        test = manual_pearson(application, experiment)
        print(f"Pearson correlation: {test:1.3f}")
        regression = stats.linregress(application_points, experiment_points)
        pearson = regression.rvalue
        slope = regression.slope
        intercept = regression.intercept

        # Plot the regression line
        x = np.linspace(min(application_points), max(application_points), 100)
        y = slope * x + intercept
        self.axs.plot(x, y, 'r', label='Linear fit')
        self.axs.text(0.05, 0.95, f"Pearson correlation: {pearson:1.3f}", transform=self.axs.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

        self.style()

    def scatter_plot():
        pass

    def get_plot(self):
        return self.fig, self.axs
        
    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
    
    def style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.set_title(title_text)
        self.axs.set_ylabel(f"Experiment {self.index_name} Contribution")
        self.axs.set_xlabel(f"Application {self.index_name} Contribution")
        self.axs.grid()

def correlation_plot(application_contributions, experiment_contributions, plot_type='scatter', integral_index_name='E', \
                     plot_redundant_reactions=False, **kwargs):
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
        nuclide-reaction-wise contributions. Default is False
        
    Returns
    -------
    - fig: matplotlib.figure.Figure, the figure containing the correlation plot"""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    application_contributions, application_nested = determine_plot_type(application_contributions, False)
    experiment_contributions, experiment_nested = determine_plot_type(experiment_contributions, False)

    if application_nested != experiment_nested:
        raise ValueError("Application and experiment contributions must have the same nested structure")
    else:
        nested = application_nested # They are be the same, so arbitrarily choose one

    # Get the list of isotopes for which contributions are available
    isotopes = list(application_contributions.keys())

    # Now convert the contributions for the application and experiment into a list of x, y pairs for plotting
    contribution_pairs = []
    if nested:
        reactions = list(application_contributions[isotopes[0]].keys())
        for isotope in isotopes:
            for reaction in reactions:
                contribution_pairs.append((application_contributions[isotope][reaction], \
                                           experiment_contributions[isotope][reaction]))
    else:
        for isotope in isotopes:
            contribution_pairs.append((application_contributions[isotope], experiment_contributions[isotope]))

    plotters = {
        'scatter': ScatterPlotter(integral_index_name, plot_redundant_reactions, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError("Unsupported plot type")

    # Create the plot and style it
    plotter.create_plot(contribution_pairs, nested)

    return plotter.get_plot()