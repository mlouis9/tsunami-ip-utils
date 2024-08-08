from ._bar_plot import _BarPlotter
from .pie_plot import _PiePlotter, _InteractivePiePlotter
from .scatter_plot import _ScatterPlotter, _InteractiveScatterPlotter, _InteractivePerturbationScatterPlotter
from .matrix_plot import _interactive_matrix_plot
from .plot_utils import _determine_plot_type
from tsunami_ip_utils.integral_indices import _add_missing_reactions_and_nuclides
import numpy as np
from uncertainties import ufloat, unumpy
from typing import List, Dict, Tuple, Union, Optional
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from tsunami_ip_utils.utils import _parse_ufloats
from tsunami_ip_utils import config
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from tsunami_ip_utils.viz.pie_plot import InteractivePieLegend
import plotly
from tsunami_ip_utils.viz.scatter_plot import InteractiveScatterLegend
from tsunami_ip_utils.viz.scatter_plot import EnhancedPlotlyFigure
from tsunami_ip_utils.viz.scatter_plot import _PerturbationScatterPlotter
from tsunami_ip_utils.viz.matrix_plot import InteractiveMatrixPlot
import copy

def contribution_plot(contributions: List[Dict], plot_type: str='bar', integral_index_name: str='E', 
                      plot_redundant_reactions: bool=True, **kwargs: dict
                      ) -> Union[Tuple[Figure, Axes], plotly.graph_objs.Figure, InteractivePieLegend]:
    """Plots the contributions to an arbitrary similarity parameter for a single experiment application pair
    
    Parameters
    ----------
    contributions
        List of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair.
    plot_type
        Type of plot to create. Default is ``'bar'`` which creates a matplotlib bar plot. Other options are:
        
        - ``'pie'``
            Creates a ``matplotlib`` pie chart.
        - ``'interactive_pie'``
            Creates an interactive Plotly sunburst chart.
    integral_index_name
        Name of the integral index being plotted. Default is ``'E'``.
    plot_redundant_reactions
        Whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is ``True``.
    kwargs
        Additional keyword arguments to pass to the plotting function.
        
    Returns
    -------
        * If ``plot_type`` is ``'bar'``, returns a tuple containing the matplotlib figure and axes objects. These
          can be used to further style the plot or save it to a file. Can be plotted by calling ``plt.show()``.
        * If ``plot_type`` is ``'pie'``, returns a Plotly figure object which can be displayed by calling ``fig.show()``
        * If ``plot_type`` is ``'interactive_pie'``, returns a :class:`.InteractivePieLegend` object which can be displayed by
          calling ``fig.show()``"""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    contributions, nested_plot = _determine_plot_type(contributions, plot_redundant_reactions)

    plotters = {
        'bar': _BarPlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'pie': _PiePlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'interactive_pie': _InteractivePiePlotter(integral_index_name, plot_redundant_reactions, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError("Unsupported plot type")

    # Create the plot and style it
    plotter._create_plot(contributions, nested_plot)

    return plotter._get_plot()


def correlation_plot(application_contributions: List[dict], experiment_contributions: List[dict], plot_type: str='scatter', 
                     integral_index_name: str='E', plot_redundant_reactions: bool=True, **kwargs: dict
                     ) -> Union[Tuple[Figure, Axes], EnhancedPlotlyFigure, InteractiveScatterLegend]:
    """Creates a correlation plot for a given application-experiment pair for which the contributions to the similarity
    parameter are given.
    
    Parameters
    ----------
    application_contributions
        List of dictionaries containing the contributions to the similarity parameter
        for the application.
    experiment_contributions
        List of dictionaries containing the contributions to the similarity parameter
        for the experiment.
    plot_type
        Type of plot to create. Default is ``'scatter'`` which creates a matplotlib scatter plot. Other options
        are ``'interactive_scatter'``, which creates a Plotly scatter plot.
    integral_index_name
        Name of the integral index being plotted. Default is ``'E'``
    plot_redundant_reactions
        Whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is ``True``.
    kwargs
        Additional keyword arguments to pass to the plotting function.

        - ``interactive_legend``: bool
            Whether to include an interactive legend in the plot. Default is ``False``. Only applicable if
            ``plot_type`` is ``'interactive_scatter'``.
        
    Returns
    -------
        * If ``plot_type`` is ``'scatter'``, returns a tuple containing the matplotlib figure and axes objects. These
          can be used to further style the plot or save it to a file. Can be plotted by calling ``plt.show()``.
        * If ``plot_type`` is ``'interactive_scatter'`` and ``interactive_legend`` is ``False``, returns a Plotly figure 
          object which can be displayed by calling ``fig.show()``
        * If ``plot_type`` is ``'interactive_scatter'`` and ``interactive_legend`` is ``True``, returns a 
          :class:`.InteractiveScatterLegend` object which can be displayed by calling ``fig.show()``
          
    Notes
    -----
    Technically, if ``plot_type`` is ``interactive_scatter`` and ``interactive_legend`` is ``False``, the function will
    return a :class:`.EnhancedPlotlyFigure` object, which is a wrapper of Plotly's ``Figure`` class, which has additional
    attributes related to the statistics of the scatter plot."""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    application_contributions, application_nested = _determine_plot_type(application_contributions, plot_redundant_reactions)
    experiment_contributions, experiment_nested = _determine_plot_type(experiment_contributions, plot_redundant_reactions)

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
        'scatter': _ScatterPlotter(integral_index_name, plot_redundant_reactions, nested, **kwargs),
        'interactive_scatter': _InteractiveScatterPlotter(integral_index_name, plot_redundant_reactions, nested, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError(f"Unsupported plot type: {plot_type}")

    # Create the plot and style it
    plotter._create_plot(copy.deepcopy(contribution_pairs), isotopes, all_reactions)

    return plotter._get_plot()


def perturbation_plot(points: List[Tuple[ufloat, ufloat]], plot_type: str='interactive_scatter') -> EnhancedPlotlyFigure:
    """Plots the perturbation points for a given application-experiment pair for which the perturbation points have already
    been calculated.
    
    Parameters
    ----------
    points: 
        List of tuples containing the perturbation points for the application-experiment pair.
    plot_type
        Type of plot to create. Default is ``'scatter'`` which creates a matplotlib scatter plot. Other options
        are ``'interactive_scatter'``, which creates a Plotly scatter plot.
        
    Returns
    -------
        A (enhanced) Plotly figure object which can be displayed by calling ``fig.show()``."""
    
    # Extracting data
    plotters = {
        'interactive_scatter': _InteractivePerturbationScatterPlotter(),
        'scatter': _PerturbationScatterPlotter(),
    }

    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError(f"Unsupported plot type: {plot_type}")

    plotter._create_plot(points)

    return plotter._get_plot()


def matrix_plot(plot_objects_array: np.ndarray, plot_type: str, labels: Optional[Dict[str, List]]=None
                ) -> InteractiveMatrixPlot:
    """Creates a Dash app to display a matrix of plots from a numpy object array of figure objects.
    
    Parameters
    ----------
    plot_objects_array
        A 2D numpy array of figure objects.
    plot_type
        Type of plot to create. Default is ``'interactive'`` which creates an interactive matrix plot (served by a Dash app). 
        Other options are ``'static'``, which creates a static matplotlib matrix plot.
    labels
        Dictionary of lists containing the labels for the rows and columns of the matrix plot. Keys are ``'expeirments'`` and 
        ``'applications'``. Default is ``None``. Most commonly, this is the sdf filename for the given application/experiment.
        
    Returns
    -------
        * If ``plot_type`` is ``'interactive'``, returns a :class:`.InteractiveMatrixPlot` object which can be run by 
          calling ``app.show()``.
        * If ``plot_type`` is ``'static'``, raises a ``NotImplementedError``."""
    if plot_type == 'interactive':
        return _interactive_matrix_plot(plot_objects_array, labels)
    elif plot_type == 'static':
        raise NotImplementedError("Static matrix plots are not yet supported")


class BlockingFigureWrapper:
    """This class wraps a ``matplotlib.figure.Figure`` instance and blocks the execution of the program until the figure
    is closed. This is useful for displaying figures in a Jupyter notebook or other interactive environments where the
    program would otherwise continue executing and close the figure immediately, i.e. ``figure.show()`` is non-blocking. 
    This is nice for creating multiple plots programmatically via a plotting function, then returning the relevant figure 
    and axs objects, after which the user can make changes and plot or save the figures individually. This functionality
    is not allowed with ``plt.show()``, because it shows all plots created in a given session (and so will show multiple plots
    instead of just a selected one)."""
    def __init__(self, figure: Figure):
        if not isinstance(figure, Figure):
            raise ValueError("Expected a matplotlib.figure.Figure instance")
        self._figure = figure

    def show(self):
        """Display the figure and block the program until the user presses Enter."""
        # Show the figure
        self._figure.show()

        input("Press Enter to continue...")

    def __getattr__(self, name):
        """Delegate all other attribute access to the wrapped matplotlib figure instance."""
        return getattr(self._figure, name)


def generate_heatmap_from_comparison(comparison_excel_path: Union[str, Path], base_fontsize: int=6
                                     ) -> Dict[str, Tuple[plt.Figure, plt.Axes]]:
    """Generates a heatmap from a comparison excel file generated by :func:`tsunami_ip_utils.comparisons.correlation_comparison`.
    
    Parameters
    ----------
    comparison_excel_path
        Path to the comparison excel file.
    base_fontsize
        Base font size for the heatmap (there is some heuristic for scaling fontsize with matrix size). Default is 6.
    
    Returns
    -------
        A dictionary containing the matplotlib figure and axes objects for each unique header in the comparison excel file."""
    
    comparison = pd.read_excel(comparison_excel_path, header=[0,1], index_col=0)
    
    # Get unique secondary headers from level 1 of the MultiIndex
    unique_headers = comparison.columns.get_level_values(1).unique()

    plot_dict = {}
    heatmap_labels = config.COMPARISON_HEATMAP_LABELS

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