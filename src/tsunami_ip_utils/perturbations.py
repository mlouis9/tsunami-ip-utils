from tsunami_ip_utils.readers import RegionIntegratedSdfReader
from tsunami_ip_utils.utils import filter_redundant_reactions
from pathlib import Path
from tsunami_ip_utils.xs import read_multigroup_xs
import pickle
import os
import tempfile
from string import Template
import subprocess
import time
from tsunami_ip_utils.utils import filter_by_nuclie_reaction_dict

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

    # Make a directory to store all cached cross section libraries if it doesn't already exist
    cache_dir = Path.home() / ".tsunami_ip_utils_cache"
    if not cache_dir.exists():
        os.mkdir(cache_dir)

    # Make a directory to store the cached perturbed multigroup libraries if it doesn't already exist
    library_name = base_library.name

    # Get the base multigroup cross sections for each nuclide reaction and the list of all available nuclide reactions for
    # caching the sampled perturbed cross sections

    if not (cache_dir / f'cached_{library_name}.pkl').exists():
        base_xs, available_nuclide_reactions = read_multigroup_xs(base_library, all_nuclide_reactions, \
                                                                return_available_nuclide_reactions=True)
        # Now read all cross sections and cache the base cross section library
        base_xs = read_multigroup_xs(base_library, available_nuclide_reactions)
        with open(cache_dir / f'cached_{library_name}.pkl', 'wb') as f:
            pickle.dump(base_xs, f)
    else:
        with open(cache_dir / f'cached_{library_name}.pkl', 'rb') as f:
            base_xs = pickle.load(f)

        # Get available nuclide reactions
        available_nuclide_reactions = { nuclide: list( reactions.keys() ) for nuclide, reactions in base_xs.items() }

        # Filter out the desired nuclide reactions
        base_xs = filter_by_nuclie_reaction_dict(base_xs, all_nuclide_reactions)


    perturbed_cache = cache_dir / f'cached_{library_name}_perturbations'
    if not perturbed_cache.exists():
        os.mkdir(perturbed_cache)

    # --------------------------------
    # Main loop for generating points
    # --------------------------------
    points = []
    for i in range(1, num_perturbations + 1):
        # Make path to the cross section perturbation factors for this random sample
        perturbed_xs_sample_name = f"Sample{i}"
        perturbed_library_sample = perturbed_library / perturbed_xs_sample_name
        
        # Cache the perturbed cross section libraries if not already cached
        perturbed_xs_cache = perturbed_cache / f'perturbed_xs_{i}.pkl'
        if not perturbed_xs_cache.exists():
            perturbed_xs = generate_and_read_perturbed_library(base_library, perturbed_library_sample, available_nuclide_reactions)
            with open(perturbed_xs_cache, 'wb') as f:
                pickle.dump(perturbed_xs, f)
        else:
            with open(perturbed_xs_cache, 'rb') as f:
                perturbed_xs = pickle.load(f)

        # Now filter out the desired nuclide reactions
        perturbed_xs = filter_by_nuclie_reaction_dict(perturbed_xs, all_nuclide_reactions)

def cache_all_libraries(base_library: Path, perturbed_library: Path, reset_cache=False):
    """Caches the base and perturbed cross section libraries for a given base library and perturbed library paths
    
    Parameters
    ----------
    - base_library: Path, path to the base cross section library
    - perturbed_library: Path, path to the cross section perturbation factors (used to generate the perturbed libraries)
    - reset_cache: bool, whether to reset the cache or not (default is False)"""
    # Read the base library, use an arbitrary nuclide reaction dict just to get the available reactions
    all_nuclide_reactions = { '92235': ['18'] } # u-235 fission

    # Make a directory to store all cached cross section libraries if it doesn't already exist
    cache_dir = Path.home() / ".tsunami_ip_utils_cache"
    if not cache_dir.exists():
        os.mkdir(cache_dir)

    if reset_cache:
        # Remove all cached cross section libraries
        for f in os.listdir(cache_dir):
            if f.endswith('_perturbations'): # A perturbed library directory
                for p in os.listdir(cache_dir / f):
                    os.remove(cache_dir / f / p)
            else:
                os.remove(cache_dir / f)

    # Cache base library if not already cached (requires reading twice)
    library_name = base_library.name
    base_library_cache = cache_dir / f'cached_{library_name}.pkl'
    print("Reading base library... ")
    if not base_library_cache.exists():
        _, available_nuclide_reactions = read_multigroup_xs(base_library, all_nuclide_reactions, \
                                                                return_available_nuclide_reactions=True)

        start = time.time()
        print("Caching base library... ", end='')
        base_xs = read_multigroup_xs(base_library, available_nuclide_reactions)
        with open( base_library_cache, 'wb') as f:
            pickle.dump(base_xs, f)
        end = time.time()
        print(f"Done in {end - start} seconds")
    else:
        with open(cache_dir / f'cached_{library_name}.pkl', 'rb') as f:
            base_xs = pickle.load(f)

        # Get available nuclide reactions
        available_nuclide_reactions = { nuclide: list( reactions.keys() ) for nuclide, reactions in base_xs.items() }
        
    
    # Make a directory to store the cached perturbed multigroup libraries if it doesn't already exist
    perturbed_cache = cache_dir / f'cached_{library_name}_perturbations'
    if not perturbed_cache.exists():
        os.mkdir( perturbed_cache )

    # ------------------------------------------
    # Main loop for caching perturbed libraries
    # ------------------------------------------
    NUM_SAMPLES = 1000 # Number of xs perturbation samples available in SCLAE
    for i in range(1, NUM_SAMPLES + 1):
        # Make path to the cross section perturbation factors for this random sample
        perturbed_xs_sample_name = f'Sample{i}'
        perturbed_library_sample = perturbed_library / perturbed_xs_sample_name

        # Cache the perturbed cross section libraries if not already cached
        perturbed_xs_cache = perturbed_cache / f'perturbed_xs_{i}.pkl'
        if not perturbed_xs_cache.exists():
            start = time.time()
            print(f"Caching perturbed library {i}... ", end='')
            perturbed_xs = generate_and_read_perturbed_library(base_library, perturbed_library_sample, available_nuclide_reactions)
            with open(perturbed_xs_cache, 'wb') as f:
                pickle.dump(perturbed_xs, f)
            
            end = time.time()
            print(f"Done in {end - start} seconds")
        else:
            print(f"Perturbed library {i} already cached...")