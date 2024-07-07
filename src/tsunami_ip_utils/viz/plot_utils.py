import socket
from tsunami_ip_utils.utils import filter_redundant_reactions, isotope_reaction_list_to_nested_dict

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