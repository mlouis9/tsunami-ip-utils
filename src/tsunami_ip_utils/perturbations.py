from tsunami_ip_utils.readers import RegionIntegratedSdfReader
from tsunami_ip_utils.utils import filter_redundant_reactions
from pathlib import Path
from tsunami_ip_utils.xs import read_multigroup_xs

"""This module is used for generating cross section perturbations and combining them with the sensitivity profiles for a given application
experiment pair to generate a similarity scatter plot"""

def generate_points(application_filename: str, experiment_filename: str, multigroup_library: Path, n_perturbations: int):
    # Read the sdfs for the application and experiment
    application = RegionIntegratedSdfReader(application_filename).convert_to_dict('numbers')
    experiment  = RegionIntegratedSdfReader(experiment_filename).convert_to_dict('numbers')

    # Filter out redundant reactions, which will introduce bias into the similarity scatter plot
    # Absorption, or "capture" as it's referred to in SCALE, total, np and nα are excluded because no data available
    redundant_reactions = ['101', '1'] 
    application = filter_redundant_reactions(application.sdf_data, redundant_reactions=redundant_reactions)
    experiment  = filter_redundant_reactions(experiment.sdf_data,  redundant_reactions=redundant_reactions)

    # Create a nuclide reaction dict for the application and experiment
    application_nuclide_reactions = {nuclide: list(reactions.keys()) for nuclide, reactions in application.items()}
    experiment_nuclide_reactions  = {nuclide: list(reactions.keys()) for nuclide, reactions in experiment.items()}

    # Take the union of the nuclide reactions for the application and experiment
    all_nuclide_reactions = application_nuclide_reactions.copy()
    all_nuclide_reactions.update(experiment_nuclide_reactions)

    # Get the base multigroup cross sections for each nuclide reaction
    base_xs = read_multigroup_xs(multigroup_library, all_nuclide_reactions, method='small')
    print(base_xs)
