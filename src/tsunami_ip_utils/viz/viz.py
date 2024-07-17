from .bar_plot import BarPlotter
from .pie_plot import PiePlotter, InteractivePiePlotter
from .scatter_plot import ScatterPlotter, InteractiveScatterPlotter, InteractivePerturbationScatterPlotter
from .matrix_plot import interactive_matrix_plot
from .plot_utils import determine_plot_type
from tsunami_ip_utils.integral_indices import _add_missing_reactions_and_nuclides
import numpy as np
from uncertainties import ufloat, unumpy
from typing import List, Dict, Tuple, Union
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from tsunami_ip_utils._utils import _parse_ufloats
from tsunami_ip_utils import config
from matplotlib.figure import Figure

def contribution_plot(contributions: List[Dict], plot_type: str='bar', integral_index_name: str='E', 
                      plot_redundant_reactions: bool=True, **kwargs: dict):
    """Plots the contributions to an arbitrary similarity parameter for a single experiment application pair
    
    Parameters
    ----------
    contributions
        list of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair
    plot_type
        type of plot to create. Default is 'bar' which creates a bar plot. Other option is 'pie' which creates
        a pie chart
    integral_index_name
        name of the integral index being plotted. Default is `'E'`
    plot_redundant_reactions
        whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is True
    kwargs
        additional keyword arguments to pass to the plotting function"""

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

    all_reactions = _add_missing_reactions_and_nuclides(application_contributions, experiment_contributions, isotopes, mode='contribution')

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
        raise ValueError(f"Unsupported plot type: {plot_type}")

    # Create the plot and style it
    plotter.create_plot(contribution_pairs, isotopes, all_reactions)

    return plotter.fig


def perturbation_plot(points: List[Tuple[ufloat, ufloat]]):
    """Plots the perturbation points for a given application-experiment pair for which the perturbation points have already
    been calculated
    
    Parameters
    ----------
    - points: list of tuple of ufloat, list of tuples containing the perturbation points for the application-experiment pair"""
    
    # Extracting data
    plotter = InteractivePerturbationScatterPlotter()
    plotter.create_plot(points)

    return plotter.get_plot()


def matrix_plot(plot_objects_array: np.ndarray, plot_type: str):
    """Creates a Dash app to display a matrix of plots from a numpy object array of figure objects."""
    if plot_type == 'interactive':
        return interactive_matrix_plot(plot_objects_array)
    elif plot_type == 'static':
        raise NotImplementedError("Static matrix plots are not yet supported")


class BlockingFigureWrapper:
    def __init__(self, figure):
        if not isinstance(figure, Figure):
            raise ValueError("Expected a matplotlib.figure.Figure instance")
        self._figure = figure

    def show(self):
        # Show the figure
        self._figure.show()

        input("Press Enter to continue...")

    def __getattr__(self, name):
        return getattr(self._figure, name)


def generate_heatmap_from_comparison(comparison_excel_path: Union[str, Path], base_fontsize: int=6
                                     ) -> Dict[str, Tuple[plt.Figure, plt.Axes]]:
    """Generates a heatmap from a comparison excel file"""
    
    comparison = pd.read_excel(comparison_excel_path, header=[0,1], index_col=0)
    
    # Get unique secondary headers from level 1 of the MultiIndex
    unique_headers = comparison.columns.get_level_values(1).unique()

    plot_dict = {}
    # Assuming `config` is a predefined dictionary that includes 'COMPARISON_HEATMAP_LABELS'
    heatmap_labels = config['COMPARISON_HEATMAP_LABELS']

    # Extract a separate 2D numpy array for each unique header in level 1
    for index, header in enumerate(unique_headers):
        # Filter the DataFrame for the current header
        filtered_df = comparison.xs(header, level=1, axis=1)
        # Convert the filtered DataFrame to a numpy array
        array = _parse_ufloats(filtered_df.to_numpy())
        
        # Now plot the array
        data = unumpy.nominal_values(array)
        plt.figure(index)
        cax = plt.matshow(data, cmap='viridis')
        cbar = plt.colorbar(cax)
        cbar.set_label(heatmap_labels[header])

        # Annotate cells with values
        font_size = min(base_fontsize, 200 / max(data.shape))
        for (i, j), val in np.ndenumerate(data):
            plt.text(j, i, f"{val:.2f}", ha='center', va='center', color='white', fontsize=font_size)

        # Set axis labels and ticks
        plt.xlabel('Application')
        plt.ylabel('Experiment')
        plt.xticks(range(data.shape[1]))
        plt.yticks(range(data.shape[0]))
        plt.tick_params(axis="x", bottom=True, labelbottom=True, top=False, labeltop=False)

        plot_dict[header] = (BlockingFigureWrapper(plt.gcf()), plt.gca())

    return plot_dict