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
from .plot_utils import find_free_port
import pickle


class PiePlotter(_Plotter):
    def __init__(self, integral_index_name, plot_redudant=False, **kwargs):
        self.index_name = integral_index_name
        self.plot_redundant = plot_redudant
    
    def _create_plot(self, contributions, nested):
        self.nested = nested
        self.fig, self.axs = plt.subplots()
        if nested:
            self.nested_pie_chart(contributions)
        else:
            self.pie_chart(contributions)

        self._style()

    def _add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.ax, sharey=self.ax)

    def _get_plot(self):
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

    def _style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)


class InteractivePiePlotter(_Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        # Check if the user wants an interactive legend
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = True
        
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def _create_plot(self, contributions, nested=True):
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

    def _get_plot(self):
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

    def _style(self):
        if self.plot_redundant and self.nested:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.fig.update_layout(title_text=title_text, title_x=0.5)  # 'title_x=0.5' centers the title


class InteractivePieLegend:
    def __init__(self, fig, df):
        """Return a flask webapp that will display an interactive legend for the sunburst chart"""
        self.fig = fig
        self.df = df
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

    def serve(self):
        """Start the Flask server to display the interactive sunburst chart"""
        port = find_free_port()
        log = open(os.devnull, 'w')
        # sys.stdout = log
        # sys.stderr = log

        # Run the Flask application in a separate thread
        thread = threading.Thread(target=lambda: self.app.run(host='localhost', port=port))
        print(f"Now running at http://localhost:{port}/")
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

    def save_state(self, filename=None):
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
    def load_state(cls, filename=None, data_dict=None):
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