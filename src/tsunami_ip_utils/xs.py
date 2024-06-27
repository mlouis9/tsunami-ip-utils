from pyparsing import *
import numpy as np
from pathlib import Path
from string import Template
import tempfile
import subprocess, os

"""This module contains the functions necessary for processing multigroup cross sections and cross section covariance matrices."""

def parse_nuclide_reaction(filename, energy_boundaries=False):
    """Reads a multigroup cross section file produced by the extractor function and returns the energy-dependent cross sections 
    as a numpy array.
    
    Parameters
    ----------
    - filename: str The filename of the cross section file
    - energy_boundaries: bool If True, the energies at which the cross sections are defined are returned as well
    
    Returns
    -------"""
    xs = {}
    with open(filename, 'r') as f:
        data = f.read()

    # ---------------------------
    # Define grammar for parsing
    # ---------------------------

    xs_data_line = Suppress(pyparsing_common.sci_real) + pyparsing_common.sci_real + Suppress(LineEnd())

    # Note that the output is formatted such that the same cross section value is printed for both energy boundaries of the group
    # to avoid duplicating the cross section data, skip every other data line
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

def parse_reactions_from_nuclide(filename, **kwargs):
    """Reads a set of reactions (given by the list of reaction mt's) from a dump of all reactions for a single nuclide from
    a SCALE library. Note this function requires that the dump included headers
    
    Parameters
    ----------
    - filename: str The filename of the dump file
    - reaction_mts: list of str The list of reaction MTs to read (required kwarg)
    - energy_boundaries: bool If True, the energies at which the cross sections are defined are returned as well (optional kwarg)
    
    Returns
    -------
    - dict A dictionary containing the cross sections for each reaction MT"""

    if 'reaction_mts' not in kwargs:
        raise ValueError("Missing required keyword argument: reaction_mts")
    
    reaction_mts = kwargs['reaction_mts']
    energy_boundaries = kwargs.get('energy_boundaries', False)

    if energy_boundaries:
        raise NotImplementedError("Energy boundaries are not yet supported for this function")

    xs = {}
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

def read_nuclide_reaction_from_multigroup_library(multigroup_library_path: Path, nuclide_zaid, reaction_mt, \
                                                  parsing_function=parse_nuclide_reaction, plot_option='plot', \
                                                    energy_boundaries=False, **kwargs):
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

def read_reactions_from_nuclide(multigroup_library_path: Path, nuclide_zaid, reaction_mts):
    """Function for reading a set of reactions from a given nuclide in a SCALE multigroup library.
    
    Parameters
    ----------
    - multigroup_library_path: Path The path to the SCALE multigroup library file
    - nuclide_zaid: str The ZAID of the nuclide
    - reaction_mts: list The list of reaction MTs to read"""

    output = read_nuclide_reaction_from_multigroup_library(multigroup_library_path, nuclide_zaid, reaction_mt='0', \
                                                           parsing_function=parse_reactions_from_nuclide, \
                                                            reaction_mts=reaction_mts, plot_option='fido')
    return output