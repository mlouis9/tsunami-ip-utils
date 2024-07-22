import re
import numpy as np
from uncertainties import ufloat
from pathlib import Path
from typing import Callable
import functools
from typing import List, Union
import tempfile
from string import Template
import subprocess
from tsunami_ip_utils.readers import read_integral_indices
import os

def _isotope_reaction_list_to_nested_dict(isotope_reaction_list, field_of_interest):
    """Converts a list of dictionaries containing isotope-reaction pairs (and some other key that represents a value of
    interest, e.g. an sdf profile or a contribution) to a nested dictionary
    
    Parameters
    ----------
    - isotope_reaction_list: list of dict, list of dictionaries containing isotope-reaction pairs and some other key
    - field_of_interest: str, the key in the dictionary that represents the value of interest

    Returns
    -------
    - nested_dict: dict, nested dictionary containing the isotope-reaction pairs and the value of interest"""

    isotope_reaction_dict = {}

    def get_atomic_number(isotope):
        matches = re.findall(r'\d+', isotope)
        return int(matches[0] if matches else -1) # Return -1 if no atomic number is found, these isotopes will be sorted last
                                                  # Should only be applicable to carbon in ENDF 7.1
    
    # Sort isotopes by atomic number so plots will have similar colors across different calls
    all_isotopes = list(set([isotope_reaction['isotope'] for isotope_reaction in isotope_reaction_list]))
    all_isotopes.sort(key=get_atomic_number)
    isotope_reaction_dict = { isotope: {} for isotope in all_isotopes }

    for isotope_reaction in isotope_reaction_list:
        isotope = isotope_reaction['isotope']
        reaction = isotope_reaction['reaction_type']
        value = isotope_reaction[field_of_interest]

        isotope_reaction_dict[isotope][reaction] = value

    return isotope_reaction_dict

def _filter_redundant_reactions(data_dict, redundant_reactions=['chi', 'capture', 'nubar', 'total']):
    """Filters out redundant reactions from a nested isotope-reaction dictionary
    
    Parameters
    ----------
    - data_dict: dict, nested dictionary containing isotope-reaction pairs
    - redundant_reactions: list of str, list of reactions to filter out"""
    return { isotope: { reaction: data_dict[isotope][reaction] for reaction in data_dict[isotope] \
                        if reaction not in redundant_reactions } for isotope in data_dict }

def _filter_by_nuclie_reaction_dict(data_dict, nuclide_reactions):
    """Filters out isotopes that are not in the nuclide_reactions dictionary
    
    Parameters
    ----------
    - data_dict: dict, nested dictionary containing isotope-reaction pairs
    - nuclide_reactions: dict, dictionary containing isotopes and their reactions"""
    return {nuclide: {reaction: xs for reaction, xs in reactions.items() if reaction in nuclide_reactions[nuclide]} \
                        for nuclide, reactions in data_dict.items() if nuclide in nuclide_reactions.keys()}


def _parse_ufloats(array_of_strings):
    """Parses a 2D array of strings into a 2D array of ufloats, assuming zero uncertainty if '+/-' is not found. """
    def to_ufloat(s):
        if isinstance(s, float):
            return ufloat(s, 0)
        
        parts = s.split('+/-')
        if len(parts) == 2:
            value, error = parts
        else:
            value = s
            error = 0
        return ufloat(float(value), float(error))
    
    return np.vectorize(to_ufloat)(array_of_strings)

def _convert_paths(func: Callable) -> Callable:
    """
    Decorator to ensure that any list argument passed to the decorated function,
    which contains strings, has those strings converted to pathlib.Path objects.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        new_args = []
        for arg in args:
            if isinstance(arg, list):
                # Convert all string items in the list to Path objects
                new_arg = [Path(item) if isinstance(item, str) else item for item in arg]
                new_args.append(new_arg)
            else:
                new_args.append(arg)
        
        new_kwargs = {k: [Path(item) if isinstance(item, str) else item for item in v] if isinstance(v, list) else v
                      for k, v in kwargs.items()}
        
        return func(*new_args, **new_kwargs)
    
    return wrapper

@_convert_paths
def _run_and_read_TSUNAMI_IP(application_filenames: Union[List[str], List[Path]], 
                            experiment_filenames: Union[List[str], List[Path]], coverx_library: str):
    """Runs TSUNAMI-IP and reads the output file using the :func:`tsunami_ip_utils.readers.read_integral_indices` function.
    
    Parameters
    ----------
    application_filenames
        List of paths to the application SDF files.
    experiment_filenames
        List of paths to the experiment SDF files.
    coverx
        The coverx library to use for TSUNAMI-IP.
        
    Returns
    -------
        Integral matrices for each integral index type. The shape of the matrices are ``(num_applications, num_experiments)``. 
        Keys are ``'c_k'``, ``'E_total'``, ``'E_fission'``, ``'E_capture'``, and ``'E_scatter'``."""
    
    current_dir = Path(__file__).parent
    template_filename = current_dir / 'input_files' / 'tsunami_ip_base_input.inp'
    with open(template_filename, 'r') as f:
        template = Template(f.read())

    # Convert filenames to strings
    application_filenames = [str(filename.absolute()) for filename in application_filenames]
    experiment_filenames  = [str(filename.absolute()) for filename in experiment_filenames]
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        # Create the input file
        input_file = template.substitute(
            application_filenames='\n'.join(application_filenames),
            experiment_filenames='\n'.join(experiment_filenames),
            coverx_library=coverx_library
        )
        f.write(input_file)
        input_filename = f.name

    # Run the TSUNAMI-IP calculation
    process = subprocess.Popen( ['scalerte', input_filename], cwd=str( Path( input_filename ).parent ) )
    process.wait()

    # Read the output file
    output_filename = f"{input_filename}.out"
    tsunami_ip_output = read_integral_indices(output_filename)

    # Clean up the temporary files
    os.remove(input_filename)
    os.remove(output_filename)

    return tsunami_ip_output

