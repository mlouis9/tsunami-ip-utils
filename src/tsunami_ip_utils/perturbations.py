"""This module is used for generating cross section perturbations and combining them with the sensitivity profiles for a given application
experiment pair to generate a similarity scatter plot"""

from tsunami_ip_utils.readers import RegionIntegratedSdfReader, read_region_integrated_h5_sdf
from tsunami_ip_utils._utils import _filter_redundant_reactions
from pathlib import Path
from tsunami_ip_utils.xs import read_multigroup_xs
import pickle
import os
import tempfile
from string import Template
import subprocess
import time
from tsunami_ip_utils._utils import _filter_by_nuclie_reaction_dict
import multiprocessing
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
from tqdm.contrib.concurrent import process_map
from typing import List, Tuple, Dict, Union
from . import config
from numpy.typing import ArrayLike

# Number of xs perturbation samples available in SCALE
NUM_SAMPLES = config['NUM_SAMPLES']

def _generate_and_read_perturbed_library(base_library: Path, perturbation_factors: Path, sample_number: int, 
                                         all_nuclide_reactions: dict) -> dict:
    """Generates and reads perturbed multigroup cross section libraries.
    
    Parameters
    ----------
    base_library
        Path to the base cross section library.
    perturbation_factors
        Path to the perturbation factors directory (corresponding to the base library).
    sample_number
        The sample number to use for generating the perturbed library. Must be from 1 - ``NUM_SAMPLES``, where 
        ``NUM_SAMPLES`` = :globalparam:`NUM_SAMPLES` is the number of perturbation factor samples provided in 
        the user's current version of SCALE. (``0`` :math:`\\leq` ``sample_number`` :math:`\\leq` :globalparam:`NUM_SAMPLES`)
    all_nuclide_reactions
        A dictionary containing the nuclide reactions that are read from the perturbed library.
    
    Returns
    -------
        A dictionary containing the perturbed cross section libraries for each nuclide reaction."""
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
        perturbation_factors=str( perturbation_factors / f'Sample{sample_number}'),
        sample_number=sample_number,
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

def generate_points(application_path: Union[Path, List[Path]], experiment_path: Union[Path, List[Path]], base_library: Path, 
                    perturbation_factors: Path, num_perturbations: int
                    ) -> Union[ List[ Tuple[ float, float ] ], 
                               np.ndarray[ List[ Tuple[ float, float ] ] ]]:
    """Generates points for a similarity scatter plot using the nuclear data sampling method.

    Parameters
    ----------
    application_path
        Path(s) to the application sensitivity profile.
    experiment_path
        Path(s) to the experiment sensitivity profile.
    base_library
        Path to the base cross section library.
    perturbation_factors
        Path to the perturbation factors directory.
    num_perturbations
        Number of perturbation points to generate.
    
    Returns
    -------
        A list of points for the similarity scatter plot.

    Notes
    -----    
    * This function will automatically cache the base cross section library and the perturbed cross section libraries in the 
      user's home directory under the ``.tsunami_ip_utils_cache`` directory if not already cached. Caching is recommended if
      perturbation points are to be generated multiple times, because the I/O overhead of dumping and reading the base and
      perturbed cross section libraries can be significant.

    * This function can also generate a matrix of points for a given set of experiment and applications for making a matrix plot
      done by passing a list of paths for the application and experiment sensitivity profiles.
        
    Theory
    ======
    The nuclear data sampling method
    involves randomly sampling cross section libraries using the AMPX tool ``clarolplus`` to calculate a perturbed cross section
    library :math:`\\Delta \\boldsymbol{\\sigma}_n = \\overline{\\boldsymbol{\\sigma}} - \\boldsymbol{\\sigma}_n`, where
    :math:`\\boldsymbol{\\sigma}_n` is the :math:`n` th randomly sampled cross section library, and 
    :math:`\\overline{\\boldsymbol{\\sigma}}` is the base cross section library, consisting of the mean values of all of the 
    cross sections (e.g. the SCALE 252-group ENDF-V7.1 library). This perturbed cross section library (a vector consisting 
    of the nuclide-reaction-group-wise perturbations to the cross sections) is dotted with the sensitivity vector for the 
    application: :math:`x_n = \\boldsymbol{S}_A \\cdot \\Delta \\boldsymbol{\\sigma}`, and the experiment
    :math:`y_n = \\boldsymbol{S}_A \\cdot \\Delta \\boldsymbol{\\sigma}`, and the resulting points :math:`(x_n, y_n)` are
    plotted on a scatter plot whose Pearson correlation coefficient is meant to correspond to the :math:`c_k` value computed
    by TSUNAMI-IP.
    """
    
    if isinstance(application_path, list) and isinstance(experiment_path, list):
        points_array = np.empty( ( len(application_path), len(experiment_path), num_perturbations, 2), dtype=object )
        for i, application in enumerate(application_path):
            for j, experiment in enumerate(experiment_path):
                if i < j: # Skip upper right triangle of the matrix (since it's symmetric)
                    continue

                points_array[i, j] = generate_points(
                    application, 
                    experiment, 
                    base_library, 
                    perturbation_factors, 
                    num_perturbations
                )
        # Now popupate the upper right triangle of the matrix with the transpose of the lower left triangle
        for i, application in enumerate(application_path):
            for j, experiment in enumerate(experiment_path):
                if i < j:
                    points_array[i, j] = points_array[j, i]

        return points_array
    elif isinstance(application_path, list) or isinstance(experiment_path, list):
        raise ValueError("Both application and experiment paths must be lists or neither.")

    # Read the sdfs for the application and experiment
    if application_path.suffix == '.sdf':
        application = RegionIntegratedSdfReader(application_path).convert_to_dict('numbers').sdf_data
        application = { nuclide_name: { reaction_name: reaction['sensitivities'] } \
                       for nuclide_name, nuclide in application.items() \
                       for reaction_name, reaction in nuclide.items() }
    elif application_path.suffix == '.h5':
        application = read_region_integrated_h5_sdf(application_path)
    else:
        raise ValueError("The application path must be either an sdf or h5 file.")
    
    if experiment_path.suffix == '.sdf':
        experiment  = RegionIntegratedSdfReader(experiment_path).convert_to_dict('numbers').sdf_data
        experiment = { nuclide_name: { reaction_name: reaction['sensitivities'] } \
                      for nuclide_name, nuclide in experiment.items() \
                      for reaction_name, reaction in nuclide.items() }
    elif experiment_path.suffix == '.h5':
        experiment = read_region_integrated_h5_sdf(experiment_path)
    else:
        raise ValueError("The experiment path must be either an sdf or h5 file.")
    

    # Filter out redundant reactions, which will introduce bias into the similarity scatter plot
    # Absorption, or "capture" as it's referred to in SCALE and total
    # redundant_reactions = ['101', '1'] 
    # application = filter_redundant_reactions(application, redundant_reactions=redundant_reactions)
    # experiment  = filter_redundant_reactions(experiment,  redundant_reactions=redundant_reactions)

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
        base_xs = _filter_by_nuclie_reaction_dict(base_xs, all_nuclide_reactions)

    perturbed_cache = cache_dir / f'cached_{library_name}_perturbations'
    if not perturbed_cache.exists():
        os.mkdir(perturbed_cache)

    # --------------------------------
    # Main loop for generating points
    # --------------------------------
    points = []
    for i in tqdm(range(1, num_perturbations + 1), desc="Generating perturbation points"):
        # Cache the perturbed cross section libraries if not already cached
        perturbed_xs_cache = perturbed_cache / f'perturbed_xs_{i}.pkl'
        if not perturbed_xs_cache.exists():
            perturbed_xs = _generate_and_read_perturbed_library(base_library, perturbation_factors, i, \
                                                               available_nuclide_reactions)
            with open(perturbed_xs_cache, 'wb') as f:
                pickle.dump(perturbed_xs, f)
        else:
            with open(perturbed_xs_cache, 'rb') as f:
                perturbed_xs = pickle.load(f)

        # Now filter out the desired nuclide reactions
        perturbed_xs = _filter_by_nuclie_reaction_dict(perturbed_xs, all_nuclide_reactions)

        # ----------------------------------------------
        # Compute S ⋅ Δσ for application and experiment
        # ----------------------------------------------
        running_total_application = 0
        running_total_experiment = 0
        for isotope, reactions in all_nuclide_reactions.items():
            for reaction in reactions:
                # First compute the cross section delta, then multiply by the sensitivity profile
                delta_xs = perturbed_xs[isotope][reaction] - base_xs[isotope][reaction]
                if isotope in application and reaction in application[isotope]:
                    running_total_application += np.dot(application[isotope][reaction], delta_xs)
                if isotope in experiment and reaction in experiment[isotope]:
                    running_total_experiment  += np.dot(experiment[isotope][reaction], delta_xs)

        points.append((running_total_application, running_total_experiment))

    return points
        

def _cache_perturbed_library(args: Tuple[int, Path, Path, int, Dict[str, List[str]], Path]) -> float:
    """Caches a single perturbed cross section library.

    Parameters
    ----------
    args
        A tuple containing all necessary components to perform the caching of a 
        perturbed library.
        
        - i (int):
            The sample number to use for generating the perturbed library.
        - base_library (Path):
            Path to the base cross section library.
        - perturbation_factors (Path):
            Path to the cross section perturbation factors (used to generate the perturbed libraries).
        - sample_number (int):
            The sample number to use for generating the perturbed library. Must be from 1 - ``NUM_SAMPLES``, where ``NUM_SAMPLES``
            is the number of perturbation factor samples provided in the user's current version of SCALE. 
            (0 :math:`\\leq` ``sample_number`` :math:`\\leq` ``NUM_SAMPLES``)
        - available_nuclide_reactions (Dict[str, List[str]]):
            A dictionary containing the nuclide reactions that are read from the perturbed library.

    Returns
    -------
        The time taken to cache the perturbed library.
    """    
    i, base_library, perturbation_factors, sample_number, available_nuclide_reactions, perturbed_cache = args
    perturbed_xs_cache = perturbed_cache / f'perturbed_xs_{i}.pkl'
    if not perturbed_xs_cache.exists():
        start = time.time()
        perturbed_xs = _generate_and_read_perturbed_library(base_library, perturbation_factors, sample_number, available_nuclide_reactions)
        with open(perturbed_xs_cache, 'wb') as f:
            pickle.dump(perturbed_xs, f)
        
        end = time.time()
        return end - start
    else:
        return 0

def cache_all_libraries(base_library: Path, perturbation_factors: Path, reset_cache: bool=False,
                        num_cores: int=multiprocessing.cpu_count()//2) -> None:
    """Caches the base and perturbed cross section libraries for a given base library and perturbed library paths.

    Parameters
    ----------
    base_library
        Path to the base cross section library.
    perturbation_factors
        Path to the cross section perturbation factors (used to generate the perturbed libraries).
    reset_cache
        Whether to reset the cache or not (default is ``False``).
    num_cores
        The number of cores to use for caching the perturbed libraries in parallel (default is half the number of cores 
        available).
        
    Returns
    -------
        This function does not return a value and has no return type.
        
    Notes
    -----
    * This function will cache the base cross section library and the perturbed cross section libraries in the user's home
      directory under the ``.tsunami_ip_utils_cache`` directory. If the user wishes to reset the cache, they can do so by
      setting the ``reset_cache`` parameter to ``True`` in the :func:`cache_all_libraries` function.
    * The caches can be `very` large, so make sure that sufficient space is available. For example, caching SCALE's 252-group
      ENDF-v7.1 library and all of the perturbed libraries currently available in SCALE (1000 samples) requires 48 GB of space,
      and for ENDF-v8.0, it requires 76 GB of space.
    * The time taken to cache the libraries can be significant (~5 hours on 6 cores, but this is hardware dependent), but when 
      caching the libraries a progress bar will be displayed with a time estimate.
    * Note, if using ``num_cores`` greater than half the number of cores available on your system, you may experience excessive
      memory usage, so proceed with caution.
    """
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
        os.mkdir(perturbed_cache)

    # ------------------------------------------
    # Main loop for caching perturbed libraries
    # ------------------------------------------

    # Create a pool of worker processes
    pool = Pool(processes=num_cores)

    # Create a list of arguments for each perturbed library
    args_list = [(i, base_library, perturbation_factors, i, available_nuclide_reactions, perturbed_cache)
                 for i in range(1, NUM_SAMPLES + 1)]

    # Use the pool to cache the perturbed libraries in parallel with a progress bar
    process_map(_cache_perturbed_library, args_list, max_workers=num_cores // 2, chunksize=1, desc='Caching perturbed libraries')

    # Close the pool
    pool.close()
    pool.join()