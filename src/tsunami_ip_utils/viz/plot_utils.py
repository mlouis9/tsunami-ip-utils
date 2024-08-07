"""Various utility functions for creating interactive plots, including functions for generating plot objects from
contributions and perturbations. These functions are used in the interactive plotting functions in the ``viz`` module."""

import socket
from tsunami_ip_utils.utils import _filter_redundant_reactions, _isotope_reaction_list_to_nested_dict
import numpy as np
from typing import Dict, List, Tuple
from uncertainties import unumpy
from pathlib import Path
from uncertainties import ufloat
import tsunami_ip_utils
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import time
from typing import Union

def _find_free_port() -> int:
    """Finds a free port on localhost for running a Flask/Dash server. This is done by creating a socket and binding it
    to an available port. The socket is then closed and the port number is returned.
    
    Returns
    -------
        A free port number."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(("", 0))  # Let the OS pick an available port
    port = s.getsockname()[1]
    s.close()
    return port

def _determine_plot_type(contributions: List[ dict ], plot_redundant_reactions: bool) -> Tuple[dict, bool]:
    """Determines whether the contributions are nuclide-wise or nuclide-reaction-wise and whether to plot redundant
    reactions or not. Then converts the contributions (list of dictionaries) to a dictionary of contributions keyed by
    isotope then by reaction type if necessary.
    
    Parameters
    ----------
    contributions
        List of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair.
    plot_redundant_reactions
        Whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions.
        
    Returns
    -------
        - contributions
            Contributions to the similarity parameter keyed by isotope then by reaction type (if necessary).
        - nested_plot
            Whether the plot is nested by nuclide then by reaction type."""
    if 'reaction_type' in contributions[0]: # Nuclide-reaction-wise contributions
        nested_plot = True # Nested plot by nuclide then by reaction type

        # Create a dictionary of contributions keyed by isotope then by reaction type
        contributions = _isotope_reaction_list_to_nested_dict(contributions, 'contribution')

        # If viewing nuclide-reaction wise contributions, it's important (at least for the visualizations in this function)
        # that if viewing the true contributions to the nuclide total, that redundant interactions (e.g. capture and fission
        # + (n, g)) and irrelevant interactions (e.g. chi and nubar) are not plotted.

        if not plot_redundant_reactions:
            # Remove redundant interactions
            contributions = _filter_redundant_reactions(contributions)
    else: # Nuclide-wise contributions
        nested_plot = False
        contributions = { contribution['isotope']: contribution['contribution'] for contribution in contributions }

    return contributions, nested_plot

def generate_plot_objects_array_from_contributions(contributions: Dict[ str, List[ unumpy.uarray ] ], 
                                                   integral_index_name: str, **kwargs: dict) -> np.ndarray:
    """Generate a matrix of plot objects (for creating a matrix plot) for the given contributions to an arbitrary integral index.
    This is valid for plots of :math:`\\Delta k/k` contributions, :math:`E` contributions, :math:`c_k` contributions, etc..
    
    Parameters
    ----------
    contributions
        Dictionary of a list of contributions to the integral index keyed by application or experiment.
    integral_index_name
        Name of the integral index being plotted.
    kwargs
        Additional keyword arguments. The following are supported:
        
        - diagonal_type (str)
            Type of plot to create on the diagonal. Default is ``'interactive_pie'`` which creates an interactive
            pie chart.
        - off_diagonal_type (str)
            Type of plot to create off the diagonal. Default is ``'interactive_scatter'`` which creates an interactive
            scatter plot.
        - interactive_contribution_legend (bool)
            Whether to make the legend interactive for the contribution plots. Default is ``True``.
        - interactive_correlation_legend (bool)
            Whether to make the legend interactive for the correlation plots. Default is ``True``.
            
    Returns
    -------
        2D numpy array of plot objects to be plotted with the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function."""
    from tsunami_ip_utils.viz import contribution_plot, correlation_plot # Import here to avoid circular import
    
    # Get options for legend interactivity and the diagonal plot type if supplied
    diagonal_type = kwargs.get('diagonal_type', 'interactive_pie')
    off_diagonal_type = kwargs.get('off_diagonal_type', 'interactive_scatter')
    interactive_correlation_legend = kwargs.get('interactive_correlation_legend', True)
    interactive_contribution_legend = kwargs.get('interactive_contribution_legend', True)
    
    num_applications = len(contributions['application'])
    num_experiments = len(contributions['experiment'])
    application_files = contributions['filenames']['application']
    experiment_files = contributions['filenames']['experiment']

    # Construct plot matrix
    plot_objects_array = np.empty( ( num_experiments, num_applications ), dtype=object )

    for application_index, application_file in enumerate(application_files):
        for experiment_index, experiment_file in enumerate(experiment_files):
            if experiment_file == application_file:
                # On the diagonal, make a contribution plot, as a correlation plot is not useful when comparing the same
                # application and experiment
                plot_objects_array[experiment_index, application_index] = \
                contribution_plot(
                    contributions['application'][application_index],
                    plot_type=diagonal_type,
                    integral_index_name=integral_index_name,
                    interactive_legend=interactive_contribution_legend,     
                )
            else:
                plot_objects_array[experiment_index, application_index] =  \
                correlation_plot(
                    contributions['application'][application_index], 
                    contributions['experiment'][experiment_index], 
                    plot_type=off_diagonal_type,
                    integral_index_name=integral_index_name, 
                    plot_redundant_reactions=True, 
                    interactive_legend=interactive_correlation_legend
                )

    return plot_objects_array

def generate_plot_objects_array_from_perturbations(points_array: np.ndarray) -> np.ndarray:
    """Generate a matrix of plot objects (for creating a matrix plot) from a numpy array of perturbation points. This is
    used for a matrix of perturbation plots only.
    
    Parameters
    ----------
    points_array
        2D numpy array of points generated from the perturbation test. Shape ``( num_applications, num_experiments )``.
            
    Returns
    -------
        2D numpy array of plot objects to be plotted with the :func:`tsunami_ip_utils.viz.viz.matrix_plot` function."""
    from tsunami_ip_utils.viz import perturbation_plot # Import here to avoid circular import
    
    # Construct plot matrix
    num_applications, num_experiments, _, _ = np.shape(points_array)
    plot_objects_array = np.empty( ( num_experiments, num_applications ), dtype=object)

    for application_index, row in enumerate(points_array):
        for experiment_index, _ in enumerate(row):
            fig = perturbation_plot(points_array[application_index, experiment_index])
            plot_objects_array[experiment_index, application_index] = fig

    return plot_objects_array

def _capture_html_as_image(html_file: Path, output_image: Path, matrix: bool=False) -> None:
    """Captures an HTML file as an image using Selenium WebDriver. This is useful for saving interactive plots as static
    images.

    Parameters
    ----------
    html_file
        Path to the HTML file to capture.
    output_image
        Path to save the output image
    matrix
        Whether or not the html file is a matrix plot. If ``True``, the function searches for a div element with a specific
        id to determine the necessary size of the webbrowser.
    """
    # Setup WebDriver
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    options.add_argument("--headless")
    driver = webdriver.Chrome(service=service, options=options)

    # Load HTML file
    html_file = Path("file://" + str(html_file))
    driver.get(str(html_file))

    def get_scroll_dimension(axis):
        return driver.execute_script(f"return document.body.parentNode.scroll{axis}")
    
    def get_element_scroll_dimension(element_id, property):
        script = f"return document.getElementById('{element_id}').scroll{property};"
        return driver.execute_script(script)

    # get the page scroll dimensions
    if matrix:
        width = get_element_scroll_dimension("matrix-plot", "Width")
    else:
        width = get_scroll_dimension("Width")
    height = get_scroll_dimension("Height")

    # set the browser window size
    driver.set_window_size(width, height)
    time.sleep(1)

    # Take screenshot
    driver.save_screenshot(str(output_image))
    driver.quit()