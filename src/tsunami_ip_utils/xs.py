from pyparsing import *
import numpy as np
from pathlib import Path
from string import Template
import tempfile
import subprocess, os

"""This module contains the functions necessary for processing multigroup cross sections and cross section covariance matrices."""

def read_multigroup_xs(filename, energy_boundaries=False):
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
    
def extract_multigroup_nuclide_reaction(multigroup_library_path: Path):
    pass

def run_multigroup_xs_reader(multigroup_library_path: Path, nuclide_zaid, reaction_mt, energy_boundaries=False):
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
        output_file_path=output_file_path
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
    output = read_multigroup_xs(output_file_path, energy_boundaries=energy_boundaries)

    # Now delete the output file
    os.remove(output_file_path)
    return output

