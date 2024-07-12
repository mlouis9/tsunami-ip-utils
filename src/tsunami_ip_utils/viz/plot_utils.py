import socket
from tsunami_ip_utils.utils import filter_redundant_reactions, isotope_reaction_list_to_nested_dict
import numpy as np
from typing import Dict, List
from uncertainties import unumpy
from pathlib import Path

def find_free_port():
    """Finds a free port on localhost for running a Flask server."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))  # Let the OS pick an available port
    port = s.getsockname()[1]
    s.close()
    return port

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

def generate_plot_objects_from_array_contributions(contributions: Dict[ str, List[ unumpy.uarray ] ], integral_index_name: str, 
                                                   **kwargs: dict) -> np.ndarray:
    """Generate a matrix of plot objects (for creating a matrix plot) for the given contributions to an arbitrary integral index.
    This is valid for plots of %Δk/k, E contributions, c_k contributions, etc..
    
    Parameters
    ----------
    contributions
        Dictionary of a list of contributions to the integral index keyed by application or experiment.
    integral_index_name
        Name of the integral index being plotted.
    kwargs
        Additional keyword arguments. The following are supported:
        * diagonal_type : str
            Type of plot to create on the diagonal. Default is ``'interactive_pie'`` which creates an interactive
            pie chart.
        * interactive_contribution_legend : bool
            Whether to make the legend interactive for the contribution plots. Default is ``True``.
        * interactive_correlation_legend : bool
            Whether to make the legend interactive for the correlation plots. Default is ``True``.
            
    Returns
    -------
        2D numpy array of plot objects to be plotted with the matrix plot function."""
    from tsunami_ip_utils.viz import contribution_plot, correlation_plot # Import here to avoid circular import
    
    # Get options for legend interactivity and the diagonal plot type if supplied
    diagonal_type = kwargs.get('diagonal_type', 'interactive_pie')
    interactive_correlation_legend = kwargs.get('interactive_correlation_legend', True)
    interactive_contribution_legend = kwargs.get('interactive_contribution_legend', True)
    
    num_applications = len(contributions['application'])
    num_experiments = len(contributions['experiment'])

    # Construct plot matrix
    plot_objects_array = np.empty( ( num_applications, num_experiments ), dtype=object )

    for application_index in range(num_applications):
        for experiment_index in range(num_experiments):
            if experiment_index == application_index:
                # On the diagonal, make a contribution plot, as a correlation plot is not useful when comparing the same
                # application and experiment
                plot_objects_array[application_index, experiment_index] = \
                contribution_plot(
                    contributions['application'][application_index],
                    plot_type=diagonal_type,
                    integral_index_name=integral_index_name,
                    interactive_legend=interactive_contribution_legend,     
                )
            else:
                plot_objects_array[application_index, experiment_index] = \
                correlation_plot(
                    contributions['application'][application_index], 
                    contributions['experiment'][experiment_index], 
                    plot_type=diagonal_type,
                    integral_index_name=integral_index_name, 
                    plot_redundant_reactions=True, 
                    interactive_legend=interactive_correlation_legend
                )

    return plot_objects_array

def generate_plot_objects_array_from_perturbations(points_array: np.ndarray) -> np.ndarray:
    """Generate a matrix of plot objects (for creating a matrix plot) for the given contributions to an arbitrary integral index.
    This is valid for plots of %Δk/k, E contributions, c_k contributions, etc..
    
    Parameters
    ----------
    points_array
        Array of points generated from the perturbation test.
            
    Returns
    -------
        2D numpy array of plot objects to be plotted with the matrix plot function."""
    from tsunami_ip_utils.viz import perturbation_plot # Import here to avoid circular import
    
    # Construct plot matrix
    num_applications, num_experiments, _, _ = np.shape(points_array)
    plot_objects_array = np.empty( ( num_applications, num_experiments ), dtype=object)

    for i, row in enumerate(points_array):
        for j, _ in enumerate(row):
            fig = perturbation_plot(points_array[i, j])
            plot_objects_array[i, j] = fig

    return plot_objects_array
    
