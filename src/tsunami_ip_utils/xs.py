"""This module contains the functions necessary for processing the binary SCALE multigroup cross section libraries
(using the AMPX tools: charmin and tabasco) into python friendly dictionaries of numpy arrays."""

from pyparsing import *
import numpy as np
from pathlib import Path
from string import Template
import tempfile
import subprocess, os
import multiprocessing
from functools import partial
import re
from typing import Union, Tuple, Dict, Any, Callable, List

def _parse_nuclide_reaction(filename: Union[str, Path], energy_boundaries: bool=False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Reads a multigroup cross section ``.txt`` file produced by the extractor function and returns the energy-dependent 
    cross sections as a numpy array.
    
    Parameters
    ----------
    filename
        The filename of the cross section file
    energy_boundaries
        If True, the energies at which the cross sections are defined are returned as well
    
    Returns
    -------
        * If ``energy_boundaries`` is ``False``, a numpy array of cross sections
        * If ``energy_boundaries`` is ``True``, a numpy array of cross sections and a numpy array of energy boundaries"""
    xs = {}
    with open(filename, 'r') as f:
        data = f.read()

    # ---------------------------
    # Define grammar for parsing
    # ---------------------------

    xs_data_line = Suppress(pyparsing_common.sci_real) + pyparsing_common.sci_real + Suppress(LineEnd())

    # Note that the output is formatted such that the same cross section value is printed for both energy boundaries of the 
    # group to avoid duplicating the cross section data, skip every other data line
    xs_parser = OneOrMore(xs_data_line + Suppress(xs_data_line))

    xs = np.array(xs_parser.parseString(data).asList())

    if energy_boundaries:
        # Define a parser that reads the energy boundaries of the groups
        energy_data_line = pyparsing_common.sci_real + Suppress(pyparsing_common.sci_real + LineEnd())
        energy_parser = OneOrMore(energy_data_line)
        energy_boundaries = np.unique(energy_parser.parseString(data).asList())
        return xs, energy_boundaries
    else:
        return xs

def _parse_reactions_from_nuclide(filename: Union[str, Path], **kwargs: dict) -> Dict[str, np.ndarray]:
    """Reads a set of reactions (given by the list of reaction mt's) from a dump of all reactions for a single nuclide from
    a SCALE library. Note this function requires that the dump included headers
    
    Parameters
    ----------
    filename
        The filename of the dump file.
    **kwargs
        Additional keyword arguments:

        - reaction_mts (List[str])
            The list of reaction MTs to read (**required** kwarg).
        - energy_boundaries (bool)
            If ``True``, the energies at which the cross sections are defined are returned as well (optional kwarg).
    
    Returns
    -------
        A dictionary containing the cross sections for each reaction MT"""

    if 'reaction_mts' not in kwargs:
        raise ValueError("Missing required keyword argument: reaction_mts")
    
    reaction_mts = kwargs['reaction_mts']
    energy_boundaries = kwargs.get('energy_boundaries', False)

    if energy_boundaries:
        raise NotImplementedError("Energy boundaries are not yet supported for this function")

    with open(filename, 'r') as f:
        data = f.read()

    # ===========================
    # Define grammar for parsing
    # ===========================

    zaid = Word(nums, max=7)
    reaction_mt = Word(nums, max=4)
    fido_field = Word(nums + '$')
    fido_subfield = Word(nums + '#')

    # -------------------------------
    # Define the header line by line
    # -------------------------------
    subfield_end = Literal('t') + LineEnd()
    other_subfield_end = Literal('e') + Literal('t') + LineEnd()

    # Define a field bundle
    bundle_line1 = Suppress(fido_field) + Suppress(zaid) + Suppress(Word(nums)) + reaction_mt
    bundle_line2 = Suppress(OneOrMore(Word(nums)))
    bundle_line3 = Suppress(fido_subfield + Word(alphanums) + OneOrMore(pyparsing_common.fnumber) + other_subfield_end)
    field_bundle = bundle_line1 + bundle_line2 + bundle_line3
    
    misc_field = fido_field + Word(nums) + Word(nums)

    header = Suppress(field_bundle) + \
             field_bundle +\
             Suppress(Optional(field_bundle)) + \
             Suppress(misc_field) + \
             Suppress(fido_subfield)

    # -------------------------------------------
    # Define the cross section data line by line
    # -------------------------------------------
    xs_data_line = Suppress(pyparsing_common.sci_real) + pyparsing_common.sci_real + Suppress(LineEnd())

    # -------------------------------------------
    # Now define the total parser for a reaction
    # -------------------------------------------
    reaction_parser = header + Group(OneOrMore(xs_data_line + Suppress(xs_data_line))) + Suppress(subfield_end)

    #--------------------------------
    # Parse the data and postprocess
    #--------------------------------
    parsed_data = reaction_parser.searchString(data)
    parsed_data = { match[0]: np.array(match[1]) for match in parsed_data }
    all_mts = parsed_data.keys()
    parsed_data = { mt: data for mt, data in parsed_data.items() if mt in reaction_mts }
    if parsed_data.keys() != set(reaction_mts):
        raise ValueError(f"Not all reaction MTs were found in the data. Missing MTs: {set(reaction_mts) - set(parsed_data.keys())}. "
                         f"This nuclide has the available MTs: {list(all_mts)}")
    return parsed_data

def _parse_from_total_library(filename: Union[str, Path], **kwargs: dict
                              ) -> Union[Dict[str, np.ndarray], Tuple[Dict[str, np.ndarray], Dict[str, list]]]:
    """Parse selected nuclide-reactions from an entire cross section library dump.
    
    Parameters
    ----------
    filename
        Path to the file containing the library dump.
    **kwargs
        Additional keyword arguments:
        
        - nuclide_reaction_dict (dict)
            A dictionary mapping nuclides to a list of reaction MTs to read (**required** kwarg).
        - return_available_nuclide_reactions (bool)
            If ``True``, the available nuclide reactions (i.e. all of the nuclide reactions for which data exists in the
            given multigroup librari) are returned as well (optional kwarg).
        - energy_boundaries (bool)
            If ``True``, the energies at which the cross sections are defined are returned as well (optional kwarg).
            
    Returns
    -------
        * If ``return_available_nuclide_reactions`` is ``False``, a doubly nested dictionary containing the cross sections for each
          nuclide-reaction pair. Keyed first by nuclide, then by reaction.
        * If ``return_available_nuclide_reactions`` is ``True``, a tuple containing the above dictionary and a dictionary containing
          all available nuclide reactions."""
    if 'nuclide_reaction_dict' not in kwargs:
        raise ValueError("Missing required keyword argument: nuclide_reaction_dict")
    
    if 'return_available_nuclide_reactions' in kwargs:
        return_available_nuclide_reactions = kwargs['return_available_nuclide_reactions']
    else:
        return_available_nuclide_reactions = False

    nuclide_reaction_dict = kwargs['nuclide_reaction_dict']
    energy_boundaries = kwargs.get('energy_boundaries', False)

    if energy_boundaries:
        raise NotImplementedError("Energy boundaries are not yet supported for this function")

    with open(filename, 'r') as f:
        data = f.read()

    # ===========================
    # Define regex patterns
    # ===========================

    header_pattern = re.compile(r'2\$\$\s+(\d+)\s+\d+\s+(\d+).*?(?=2##)', re.DOTALL)
    xs_data_pattern = re.compile(r'4##\s*((?:\s*[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?\s+[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?\s*\n)+)', re.MULTILINE)

    # ===========================
    # Parse the data
    # ===========================

    parsed_data_dict = {}
    all_nuclide_reactions = {}
    for match in header_pattern.finditer(data):
        nuclide = match.group(1)
        reaction = match.group(2)
        header_end = match.end()

        # Now record all available nuclide reactions
        if nuclide not in all_nuclide_reactions:
            all_nuclide_reactions[nuclide] = []
        all_nuclide_reactions[nuclide].append(reaction)

        # Check if the nuclide and reaction are in the nuclide_reaction_dict
        if nuclide in nuclide_reaction_dict and reaction in nuclide_reaction_dict[nuclide]:
            xs_data_match = xs_data_pattern.search(data, header_end)
            if xs_data_match:
                xs_data_text = xs_data_match.group(1)
                lines = xs_data_text.strip().split('\n')
                xs_data = [float(line.split()[1]) for line in lines[::2]]  # Extract second column and skip every other row
                if nuclide not in parsed_data_dict:
                    parsed_data_dict[nuclide] = {}
                parsed_data_dict[nuclide][reaction] = np.array(xs_data)

    # ========================================
    # Check for missing nuclides and reactions
    # ========================================

    nuclides_not_found = set(nuclide_reaction_dict.keys()) - set(parsed_data_dict.keys())
    reactions_not_found = {}
    for nuclide, reactions in nuclide_reaction_dict.items():
        # Skip the nuclides that aren't found
        if nuclide in nuclides_not_found:
            continue
        for reaction in reactions:
            if reaction not in parsed_data_dict[nuclide]:
                if nuclide not in reactions_not_found:
                    reactions_not_found[nuclide] = []
                reactions_not_found[nuclide].append(reaction)

    if len(nuclides_not_found) > 0 or reactions_not_found != {}:
        raise ValueError(f"Not all requested reactions were found in the data. Missing reactions: {reactions_not_found}. "
                         f"And missing nuclides: {nuclides_not_found}")


    if return_available_nuclide_reactions:
        # Remove the nuclide '0' from the list of available nuclides if it exists
        all_nuclide_reactions.pop('0', None)
        return parsed_data_dict, all_nuclide_reactions
    else:
        return parsed_data_dict

def _read_nuclide_reaction_from_multigroup_library(multigroup_library_path: Path, nuclide_zaid: str, reaction_mt: str, 
                                                   parsing_function: Callable[[Union[str, Path], bool, dict], Any]=_parse_nuclide_reaction, 
                                                   plot_option: str='plot', energy_boundaries: bool=False, **kwargs: dict) -> Any:
    """Uses SCALE to dump a binary multigroup library to a text file, and then calls the specified parsing function on the 
    output file.
    
    Parameters
    ----------
    multigroup_library_path
        The path to the SCALE multigroup library file.
    nuclide_zaid
        The ZAID of the nuclide.
    reaction_mt
        The reaction MT to read.
    parsing_function
        The function to call on the output file. This function takes the filename of the output file as its first argument, and
        whether or not energy boundaries should be returned as its second argument. Additional keyword arguments can be passed
        to the parsing function using the kwargs argument.
    plot_option
        The plot option to use when running the MG reader. For example ``'plot'`` or ``'fido'``.
    energy_boundaries
        If True, the energies at which the cross sections are defined are returned as well.
    
    Returns
    -------
        An output that is the result of the parsing function"""
    # Get the directory of the current file
    current_dir = Path(__file__).parent

    # Construct the path to the input file
    file_path = current_dir / 'input_files' / 'MG_reader.inp'

    # Create a tempfile for storing the output file of the MG reader dump.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name
    
    # Read the MG reader input template file
    with open(file_path, 'r') as f:
        template = Template(f.read())
        
    # Substitute the input file template variables
    input_file = template.safe_substitute(
        nuclide_zaid=nuclide_zaid, 
        reaction_mt=reaction_mt, 
        multigroup_library_path=multigroup_library_path,
        output_file_path=output_file_path,
        plot_option=plot_option
    )

    # Write the input file to another tempfile
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as input_temp_file:
        input_temp_file.write(input_file)
        input_temp_file_path = input_temp_file.name

    # Run the executable
    command = ['scalerte', input_temp_file_path]

    proc = subprocess.Popen(command)
    proc.wait()

    # Now delete the input file
    os.remove(input_temp_file_path)

    # Read the output file
    output = parsing_function(output_file_path, energy_boundaries=energy_boundaries, **kwargs)

    # Now delete the output file
    os.remove(output_file_path)
    return output

def _read_reactions_from_nuclide(multigroup_library_path: Path, nuclide_zaid: str, reaction_mts: List[str]
                                 ) -> Dict[str, np.ndarray]:
    """Function for reading a set of reactions from a given nuclide in a SCALE multigroup library.
    
    Parameters
    ----------
    multigroup_library_path
        The path to the SCALE multigroup library file.
    nuclide_zaid
        The ZAID of the nuclide.
    reaction_mts
        The list of reaction MTs to read.
    
    Returns
    -------
        A dictionary containing the cross sections for each reaction MT"""
    output = _read_nuclide_reaction_from_multigroup_library(multigroup_library_path, nuclide_zaid, reaction_mt='0', \
                                                           parsing_function=_parse_reactions_from_nuclide, \
                                                            reaction_mts=reaction_mts, plot_option='fido')
    return output

def read_multigroup_xs(multigroup_library_path: Path, nuclide_zaid_reaction_dict: Dict[str, Dict[str, str]], 
                       num_processes: int=multiprocessing.cpu_count(), return_available_nuclide_reactions: bool=False
                       ) -> Union[Dict[str, Dict[str, np.ndarray]], Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, list]]]:
    """Function for reading a set of reactions from a given nuclide in a SCALE multigroup library.
    
    Parameters
    ----------
    multigroup_library_path
        The path to the SCALE multigroup library file.
    nuclide_zaid_reaction_dict
        A dictionary mapping nuclide ZAIDs to a list of reaction MTs to read.
    num_processes
        The number of processes to use for reading the library. If None, the number of processes is set to the number of cores.
    return_available_nuclide_reactions
        If True, the available nuclide reactions (i.e. all of the nuclide reactions for which data exists in the given 
        multigroup library) are returned as well.

    Returns
    -------
        * If ``return_available_nuclide_reactions`` is ``False``, a dictionary containing the cross sections (as numpy arrays) 
          for each nuclide-reaction pair. Keyed first by nuclide, then by reaction.
        * If ``return_available_nuclide_reactions`` is ``True``, a tuple containing the above dictionary and a 
          dictionary containing all available nuclide reactions is returned.
    """

    NUCLIDE_THRESHOLD = 50 # Number of nuclides after which the large method is more performant and hence is used
    CORE_THRESHOLD = 2 # The large method is more performant if the number of cores is smaller than this number

    num_nuclides = len(list(nuclide_zaid_reaction_dict.keys()))
    use_small_method = ( num_nuclides < NUCLIDE_THRESHOLD ) and ( num_processes >= CORE_THRESHOLD )
    if use_small_method and not return_available_nuclide_reactions: # This method is slow but works for small amounts of nuclide reactions or a large amount of cores
        pool = multiprocessing.Pool(processes=num_processes)

        # Create a partial function with the common arguments
        read_reactions_partial = partial(_read_reactions_from_nuclide, multigroup_library_path)

        # Distribute the function calls among the processes
        results = pool.starmap(read_reactions_partial, nuclide_zaid_reaction_dict.items())

        # Close the pool and wait for the processes to finish
        pool.close()
        pool.join()

        # Convert the results to a dictionary
        output = dict(zip(nuclide_zaid_reaction_dict.keys(), results))

        return output
    else: # This method is faster (as in there's less scale run overhead) but requires the entire library to be read and is serial
        # If the user wants the available nuclide reactions, then we need to create a partial function which adds the appropriate
        # keyword argument to the parsing function
        if return_available_nuclide_reactions:
            parse_function = partial(_parse_from_total_library, return_available_nuclide_reactions=True)
        else:
            parse_function = _parse_from_total_library
        
        output = _read_nuclide_reaction_from_multigroup_library(multigroup_library_path, nuclide_zaid='0', reaction_mt='0', \
                                                        parsing_function=parse_function, \
                                                        nuclide_reaction_dict=nuclide_zaid_reaction_dict, plot_option='fido')
        return output