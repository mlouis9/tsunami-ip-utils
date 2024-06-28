from tsunami_ip_utils.readers import RegionIntegratedSdfReader
from tsunami_ip_utils.utils import filter_redundant_reactions
from pathlib import Path
from tsunami_ip_utils.xs import read_multigroup_xs
import pickle
import os
import tempfile
from string import Template
import subprocess

"""This module is used for generating cross section perturbations and combining them with the sensitivity profiles for a given application
experiment pair to generate a similarity scatter plot"""

def generate_and_read_perturbed_library(base_library: Path, perturbation_factors: Path, all_nuclide_reactions: dict):
    # Read the SCALE input template
    current_dir = Path(__file__).parent
    template_filename = current_dir / 'input_files' / 'generate_perturbed_library.inp'
    with open(template_filename, 'r') as f:
        template = Template(f.read())

    # Open a temporary file to store the output file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name

    # Substitute the input file template values
    input_file = template.substitute(
        base_library=str(base_library),
        perturbation_factors=str(perturbation_factors),
        output=output_file_path
    )

    # Open a temporary file to store the input file, then one to store the output file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as input_temp_file:
        input_temp_file.write(input_file)
        input_temp_file_path = input_temp_file.name

    # Run the executable
    command = ['scalerte', input_temp_file_path]

    proc = subprocess.Popen(command)
    proc.wait()

    # Now delete the input file
    os.remove(input_temp_file_path)

    # Read the perturbed library
    perturbed_xs = read_multigroup_xs(Path(output_file_path), all_nuclide_reactions)

    # Now delete the output file
    os.remove(output_file_path)

    return perturbed_xs

def generate_points(application_filename: str, experiment_filename: str, base_library: Path, perturbed_library: Path, \
                     num_perturbations: int):
    # Read the sdfs for the application and experiment
    application = RegionIntegratedSdfReader(application_filename).convert_to_dict('numbers')
    experiment  = RegionIntegratedSdfReader(experiment_filename).convert_to_dict('numbers')

    # Filter out redundant reactions, which will introduce bias into the similarity scatter plot
    # Absorption, or "capture" as it's referred to in SCALE, total, (n,p) and (n,α) are excluded because no data available
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
    base_xs = read_multigroup_xs(base_library, all_nuclide_reactions)

    # Make a directory to store the cached perturbed multigroup libraries if it doesn't already exist
    current_dir = Path(__file__).parent
    library_name = base_library.name
    cache_directory_name = str(current_dir / f'cached_{library_name}_perturbations')
    if not os.path.isdir(cache_directory_name):
        os.mkdir(cache_directory_name)

    # --------------------------------
    # Main loop for generating points
    # --------------------------------
    points = []
    perturbed_xs_sample_names = [f'Sample{i}' for i in range(1, num_perturbations+1)] 
    for i in range(num_perturbations):
        # Make path to the cross section perturbation factors for this random sample
        perturbed_xs_sample_name = perturbed_xs_sample_names[i]
        perturbed_library_sample = perturbed_library / perturbed_xs_sample_name
        
        # Cache the perturbed cross section libraries if not already cached
        perturbed_xs_cache = str(current_dir / f'cached_{library_name}_perturbations' / f'perturbed_xs_{i}.pkl')
        if not os.path.exists(perturbed_xs_cache):
            perturbed_xs = generate_and_read_perturbed_library(base_library, perturbed_library_sample, all_nuclide_reactions)
            with open(perturbed_xs_cache, 'wb') as f:
                pickle.dump(perturbed_xs, f)
        else:
            with open(perturbed_xs_cache, 'rb') as f:
                perturbed_xs = pickle.load(f)

        