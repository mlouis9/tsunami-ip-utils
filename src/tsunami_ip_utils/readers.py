import numpy as np
from pyparsing import *
from uncertainties import unumpy, ufloat
import h5py
from pathlib import Path
from . import config
from typing import Tuple, List, Union, Dict

ParserElement.enablePackrat()

class SdfReader:
    """A class for reading TSUNAMI-B Sentitivity Data Files (SDFs, i.e. ``.sdf`` files produced by TSUNAMI-3D monte carlo
    transport simulations).
    
    The format for TSUNAMI-B SDF files is given `here <https://scale-manual.ornl.gov/tsunami-ip-appAB.html#format-of-tsunami-b-
    sensitivity-data-file>`_.

    Notes
    -----
    The SDF reader currently does not support TSUNAMI-A formatted SDF files (produced by deterministic transport simulations).
    """
    energy_boundaries: np.ndarray
    """Boundaries for the energy groups"""

    sdf_data: List[dict]
    """List of dictionaries containing the sensitivity profiles and other derived/descriptive data. The dictionary
            keys are given by ``SDF_DATA_NAMES`` = :globalparam:`SDF_DATA_NAMES`."""
    def __init__(self, filename: Union[str, Path]):
        """Create a TSUNAMI-B SDF reader object from the given filename
        
        Parameters
        ----------
        filename
            Path to the sdf file."""
        self.energy_boundaries, self.sdf_data = self._read_sdf(filename)
        
    def _read_sdf(self, filename: Union[str, Path]) -> Tuple[np.ndarray, List[dict]]:
        """Reads the SDF file and returns a dictionary of nuclide-reaction pairs and energy-dependent
        sensitivities (with uncertainties)
        
        Parameters
        ----------
        Filename
            Path to the sdf file.
        
        Returns
        -------
        energy_boundaries
            Energy boundaries for the energy groups.
        sdf_data
            List of dictionaries containing the sensitivity profiles and other derived/descriptive data. The dictionary
            keys are given by ``SDF_DATA_NAMES`` = :globalparam:`SDF_DATA_NAMES`."""
        with open(filename, 'r') as f:
            data = f.read()

        # ========================
        # Get sensitivity profiles
        # ========================

        # Get number of energy groups
        unused_lines = SkipTo(pyparsing_common.integer + "number of neutron groups")
        num_groups_parser = Suppress(unused_lines) + pyparsing_common.integer + Suppress("number of neutron groups")
        num_groups = num_groups_parser.parseString(data)[0]

        data_line = Group(OneOrMore(pyparsing_common.sci_real))
        data_block = OneOrMore(data_line)

        unused_lines = SkipTo("energy boundaries:")
        energy_boundaries = Suppress(unused_lines + "energy boundaries:") + data_block
        energy_boundaries = np.array(energy_boundaries.parseString(data)[0])

        # ------------------
        # SDF profile parser
        # ------------------
        atomic_number = Word(nums)
        element = Word(alphas.lower(), exact=1) 
        isotope_name = Combine(element + '-' + atomic_number)

        # Grammar for sdf header
        reaction_type = Word(alphanums + ',\'')
        zaid = Word(nums, max=6)
        reaction_mt = Word(nums, max=4)

        # Lines of the sdf header
        sdf_header_first_line = isotope_name + reaction_type + zaid + reaction_mt + Suppress(LineEnd())
        sdf_header_second_line = pyparsing_common.signed_integer + pyparsing_common.signed_integer + Suppress(LineEnd())
        sdf_header_third_line = pyparsing_common.sci_real + pyparsing_common.sci_real + pyparsing_common.signed_integer + \
                                    pyparsing_common.signed_integer + LineEnd()
        
        # This line contains the total energy integrated sensitivity data for the given profile, along with uncertainties, etc.
        sdf_data_first_line = Group(pyparsing_common.sci_real + pyparsing_common.sci_real) + \
                              pyparsing_common.sci_real + Group(pyparsing_common.sci_real + pyparsing_common.sci_real) + \
                              Suppress(LineEnd())
        
        # The total sdf header
        sdf_header = sdf_header_first_line + sdf_header_second_line + Suppress(sdf_header_third_line)

        # SDF profile data
        sdf_data_block = OneOrMore(data_line)
        sdf_data = sdf_header + sdf_data_first_line + sdf_data_block
        sdf_data = sdf_data.searchString(data)

        # Now break the data blocks into two parts: the first part is the sensitivities and the second is the uncertainties
        # NOTE: The sensitivities are read from largest to smallest energy group, so we need to reverse the order for them
        # to correspond to the cross section values
        sdf_data = [match[:-1] + [np.array(match[-1][:num_groups])[::-1], np.array(match[-1][num_groups:])[::-1]] for match in sdf_data]

        # -------------------------------------------------
        # Now parse each result into a readable dictionary
        # -------------------------------------------------

        # NOTE: sum_opposite_sign_groupwise_sensitivities refers to the groupwise sensitivities with opposite sign to the
        # integrated sensitivity coefficient
        names = config['SDF_DATA_NAMES']
        sdf_data = [dict(zip(names, match)) for match in sdf_data]

        # Convert the sensitivities and uncertainties to uncertainties.ufloat objects
        for match in sdf_data:
            match["sensitivities"] = unumpy.uarray(match['sensitivities'], match['uncertainties'])
            match["energy_integrated_sensitivity"] = \
                ufloat(match['energy_integrated_sensitivity'][0], match['energy_integrated_sensitivity'][1])
            match["sum_opposite_sign_groupwise_sensitivities"] = \
                ufloat(match['sum_opposite_sign_groupwise_sensitivities'][0], match['sum_opposite_sign_groupwise_sensitivities'][1])
            del match["uncertainties"]

        return energy_boundaries, sdf_data
        
class RegionIntegratedSdfReader(SdfReader):
    """Reads region integrated TSUNAMI-B sensitivity data files produced by TSUNAMI-3D. Useful when the spatial dependence of 
    sensitivty is not important."""
    filename: Union[str, Path]
    """Path to the sdf file."""
    sdf_data: Union[List[dict], Dict[str, Dict[str, dict]]]
    """Collection of region integrated sdf profiles. This can either be a list or a twice-nested dictionary (keyed by first by 
    nuclide and then reaction type) of dictionaries keyed by ``SDF_DATA_NAMES`` = :globalparam:`SDF_DATA_NAMES`."""
    def __init__(self, filename: Union[str, Path]):
        """Create a TSUNAMI-B region integrated SDF reader object from the given filename
        
        Parameters
        ----------
        filename:
            Path to the sdf file."""
        super().__init__(filename)
        
        # Now only return the region integrated sdf profiles
        # i.e. those with zone number and zone volume both equal to 0
        self.filename = filename
        self.sdf_data = [ match for match in self.sdf_data if match['zone_number'] == 0 and match['zone_volume'] == 0 ]
    
    def convert_to_dict(self, key: str='names'):
        """Converts the sdf data into a dictionary keyed by nuclide-reaction pair or by ZAID and reaction MT.
        
        Parameters
        ----------
        key
            The key to use for the dictionary. Default is ``'names'`` which uses the isotope name and reaction type
            if ``'numbers'`` is supplied instead then the ZAID and reaction MT are used."""
        # Since data is region and mixture integrated we can assume that there is only one entry for each nuclide-reaction pair
        if type(self.sdf_data) == dict:
            return self
        
        sdf_data_dict = {}
        for match in self.sdf_data:
            if key == 'names':
                nuclide = match['isotope']
                reaction_type = match['reaction_type']
            elif key == 'numbers':
                nuclide = match['zaid']
                reaction_type = match['reaction_mt']
            
            if nuclide not in sdf_data_dict:
                sdf_data_dict[nuclide] = {}

            sdf_data_dict[nuclide][reaction_type] = match
            
        self.sdf_data = sdf_data_dict
        return self
    
    def get_sensitivity_profiles(self, reaction_type: str='all') -> List[unumpy.uarray]:
        """Returns the sensitivity profiles for each nuclide-reaction pair in a list in the order they appear in the ``sdf_data``.
        
        Parameters
        ----------
        reaction_type
            The type of reaction to consider. Default is 'all' which considers all reactions.
        
        Returns
        -------
            List of sensitivity profiles for each nuclide-reaction pair."""
        if type(self.sdf_data) == list:
            if reaction_type == 'all':
                return [ data['sensitivities'] for data in RegionIntegratedSdfReader(self.filename).sdf_data ]
            else:
                return [ data['sensitivities'] for data in RegionIntegratedSdfReader(self.filename).sdf_data \
                        if data['reaction_type'] == reaction_type ]
        elif type(self.sdf_data) == dict:
            if reaction_type == 'all':
                return [ reaction['sensitivities'] for isotope in self.sdf_data.values() for reaction in isotope.values() ]
            else:
                return [ reaction['sensitivities'] for isotope in self.sdf_data.values() for reaction in isotope.values() \
                         if reaction['reaction_type'] == reaction_type ]     
        else:
            raise ValueError("Invalid data type for sdf_data. How did that happen?")

def read_covariance_matrix(filename: str):
    pass

def _read_ck_contributions(filename: str):
    pass

def read_uncertainty_contributions(filename: Union[str, Path]) -> Tuple[List[dict], List[dict]]:
    """Reads the output file from TSUNAMI and returns the uncertainty contributions for each nuclide-reaction
    covariance.
    
    Parameters
    ----------
    filename
        Path to the TSUNAMI output file.

    Returns
    -------
        * isotope_totals
            List of dictionaries containing the nuclide-wise contributions.
        * isotope_reaction
            List of dictionaries containing the nuclide-reaction pairs and the contributions.
    """
    with open(filename, 'r') as f:
        data = f.read()

    # ----------------------------------------------------
    # Define the formattting that precedes the data table
    # ----------------------------------------------------
    table_identifier = Literal("contributions to uncertainty in k-eff (% delta-k/k) by individual energy covariance matrices:")
    skipped_lines = SkipTo(table_identifier)
    pre_header = Literal("covariance matrix") + LineEnd()
    header = Word("nuclide-reaction") + Word("with") + Word("nuclide-reaction") + Word("% delta-k/k due to this matrix")
    dash_separator = OneOrMore(OneOrMore('-'))

    # ----------------------
    # Define the data lines
    # ----------------------

    # Define the grammar for the nuclide-reaction pair
    atomic_number = Word(nums)
    element = Word(alphas.lower(), exact=1) 
    isotope_name = Combine(element + '-' + atomic_number)
    reaction_type = Word(alphanums + ',\'')

    data_line = Group(isotope_name + reaction_type + isotope_name + reaction_type + \
                pyparsing_common.sci_real + Suppress(Literal("+/-")) + pyparsing_common.sci_real + Suppress(LineEnd()))
    data_block = OneOrMore(data_line)

    # -------------------------------------------
    # Define the total parser and parse the data
    # -------------------------------------------
    data_parser = Suppress(skipped_lines) + Suppress(table_identifier) + Suppress(pre_header) + Suppress(header) + \
                    Suppress(dash_separator) + data_block

    # -------------------------------------------------------------------------------
    # Now convert the data into isotope wise and isotope-reaction wise contributions
    # -------------------------------------------------------------------------------
    isotope_reaction = []
    for match in data_parser.parseString(data):
        isotope_reaction.append({
            'isotope': f'{match[0]} - {match[2]}',
            'reaction_type': f'{match[1]} - {match[3]}',
            'contribution': ufloat(match[4], match[5])
        })

    # Now calculate nuclide totals by summing the contributions for each nuclide via total = sqrt((pos)^2 - (neg)^2)
    isotope_totals = {}
    for data in isotope_reaction:
        # First add up squared sums of all reaction-wise contributions
        isotope = data['isotope']
        contribution = data['contribution']
        if isotope not in isotope_totals.keys():
            isotope_totals[isotope] = ufloat(0,0)

        if contribution < 0:
            isotope_totals[isotope] -= ( data['contribution'] )**2
        else:
            isotope_totals[isotope] += ( data['contribution'] )**2
        
    # Now take square root of all contributions
    for isotope, total in isotope_totals.items():
        isotope_totals[isotope] = total**0.5

    # Now convert into a list of dictionaries
    isotope_totals = [ {'isotope': isotope, 'contribution': total} for isotope, total in isotope_totals.items() ]

    return isotope_totals, isotope_reaction


def read_integral_indices(filename: Union[str, Path]) -> Dict[str, unumpy.uarray]:
    """Reads the output file from TSUNAMI-IP and returns the integral values for each application.

    Notes
    -----
    Currently, this function and only reads :math:`c_k`, :math:`E_{\\text{total}}`, 
    :math:`E_{\\text{fission}}`, :math:`E_{\\text{capture}}`, and :math:`E_{\\text{scatter}}`. If any of these are missing
    from the output file, the function will raise an error. To ensure these are present, please include at least
    
    ::

        read parameters
            e c
        end parameters

    in the TSUNAMI-IP input file.

    Parameters
    ----------
    filename
        Path to the TSUNAMI-IP output file.
    
    Returns
    -------
        Integral matrices for each integral index type. The shape of the matrices are ``(num_experiments, num_applications)``. 
        Keys are ``'C_k'``, ``'E_total'``, ``'E_fission'``, ``'E_capture'``, and ``'E_scatter'``."""

    with open(filename, 'r') as f:
        data = f.read()

    # Define the Integral Values parser
    dashed_line = OneOrMore("-")
    header = Literal("Integral Values for Application") + "#" + pyparsing_common.integer + LineEnd() + dashed_line
    table_header = Literal("Experiment") + Literal("Type") + Literal("Value") + Literal("s.d.") + \
                    Optional( Literal("xsec unc %") + Literal("s.d.") ) + Literal("c(k)") + \
                    Literal("s.d.") + Literal("E") + Literal("s.d.") + Literal("E(fis)") + Literal("s.d.") + Literal("E(cap)") + \
                    Literal("s.d.") + Literal("E(sct)") + Literal("s.d.") + LineEnd() + OneOrMore(dashed_line)
    
    # Define characters allowed in a filename (all printables except space)
    non_space_printables = ''.join(c for c in printables if c != ' ')
    sci_num = pyparsing_common.sci_real
    space_as_zero = White(' ', min=8).setParseAction(lambda: 0.0)  # Parsing nine spaces as zero
    value_or_space = MatchFirst([ space_as_zero, sci_num ])

    data_line = Group(Suppress(pyparsing_common.integer + Word(non_space_printables) + Opt(Literal('-') + Word(alphas)) + \
                               Word(alphas)) + \
                        Opt( Suppress( sci_num + value_or_space ) ) + \
                        Group( sci_num + value_or_space ) + Group( sci_num + value_or_space ) + Group( sci_num + value_or_space ) + \
                        Group( sci_num + value_or_space ) + Group( sci_num + value_or_space ) + Group( sci_num + value_or_space ))
    data_block = OneOrMore(data_line)
    integral_values = Suppress(header + table_header) + data_block
    parsed_integral_values = integral_values.searchString(data)

    # Parse the integral value tables into a uarray
    num_applications = len(parsed_integral_values)
    num_experiments = len(parsed_integral_values[0]) - 1 # First row seems to be a repeat, i.e. in the output it's "experiment 0"
    
    integral_matrices = {}
    integral_matrix = unumpy.umatrix( np.zeros( (num_experiments, num_applications) ), 
                                      np.zeros( (num_experiments, num_applications) ) )
    
    # Initialize the integral matrices
    C_k       = np.copy(integral_matrix)
    E_total   = np.copy(integral_matrix)
    E_fission = np.copy(integral_matrix)
    E_capture = np.copy(integral_matrix)
    E_scatter = np.copy(integral_matrix)

    # Now populate the integral matrices from the parsed output
    for match_index, match in enumerate(parsed_integral_values):
        for row_index, row in enumerate(match[1:]):
            C_k[row_index, match_index]       = ufloat(row[1][0], row[1][1])
            E_total[row_index, match_index]   = ufloat(row[2][0], row[2][1])
            E_fission[row_index, match_index] = ufloat(row[3][0], row[3][1])
            E_capture[row_index, match_index] = ufloat(row[4][0], row[4][1])
            E_scatter[row_index, match_index] = ufloat(row[5][0], row[5][1])

    integral_matrices.update({
        "C_k": C_k,
        "E_total": E_total,
        "E_fission": E_fission,
        "E_capture": E_capture,
        "E_scatter": E_scatter
    })

    return integral_matrices

def read_region_integrated_h5_sdf(filename: Path) -> Dict[unumpy.uarray]:
    """Reads all region integrated SDFs from a HDF5 (``.h5``) formatted TSUNAMI-B sdf file and returns a dictionary of 
    the data
    
    Parameters
    ----------
    filename
        Path to the .h5 SDF file (e.g. ``my_model.sdf.h5``)
        
    Returns
    -------
    sdf_data
        Dictionary of the region integrated SDF data. The dictionary is twice-nested, and keyed first by nuclide, then by
        reaction type. The values are the sensitivity profiles with uncertainties."""
    with h5py.File(filename, 'r') as f:
        # Get the region integrated sdf's, i.e. where unit=0
        region_integrated_indices = np.where(f['unit'][:] == 0)[0]
        mts = f['mt'][region_integrated_indices]
        zaids = f['nuclide_id'][region_integrated_indices]
        sensitivities = f['profile_values'][region_integrated_indices, :]
        uncertainties = f['profile_sigmas'][region_integrated_indices, :]

        # Now convert the data into a dictionary
        unique_zaids = set(zaids)
        sdf_data = {str(zaid): {} for zaid in unique_zaids}
        for index in region_integrated_indices:
            sdf_data[str(zaids[index])][str(mts[index])] = unumpy.uarray(sensitivities[index][::-1], uncertainties[index][::-1])
    return sdf_data