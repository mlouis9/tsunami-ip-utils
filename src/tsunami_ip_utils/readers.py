import numpy as np
from pyparsing import *
from uncertainties import unumpy, ufloat
import h5py
from pathlib import Path
from . import config
from typing import Tuple, List, Union, Dict
from string import Template
from tempfile import NamedTemporaryFile
import subprocess
import os

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
        element = Word(alphas.lower(), max=2) 
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
    """Collection of region integrated sdf profiles. Only includes SDF profiles with ``zone_number == 0`` and ``zone_volume == 0`` 
    This can either be a list or a twice-nested dictionary (keyed by first by 
    nuclide and then reaction type) of dictionaries keyed by ``SDF_DATA_NAMES`` = :globalparam:`SDF_DATA_NAMES`."""
    
    def __init__(self, filename: Union[str, Path]):
        """Create a TSUNAMI-B region integrated SDF reader object from the given filename
        
        Parameters
        ----------
        filename:
            Path to the sdf file.
            
        Examples
        --------
        >>> reader = RegionIntegratedSdfReader('tests/example_files/sphere_model_1.sdf')
        >>> reader.sdf_data[0]
        {'isotope': 'u-234', 'reaction_type': 'total', 'zaid': '92234', 'reaction_mt': '1', 'zone_number': 0, 'zone_volume': 0, 'energy_integrated_sensitivity': 0.008984201+/-0.0002456359, 'abs_sum_groupwise_sensitivities': 0.009319556, 'sum_opposite_sign_groupwise_sensitivities': -0.0001676775+/-4.959921e-05, 'sensitivities': array([0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -4.393147e-10+/-4.39204e-10, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -1.0871e-09+/-1.086821e-09, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, -8.666447e-10+/-8.664239e-10, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, -2.021382e-09+/-2.020862e-09, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -2.351273e-11+/-2.350687e-11, 1.707797e-08+/-1.707371e-08,
               -7.172835e-10+/-7.171047e-10, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -1.89505e-09+/-1.894573e-09, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, -5.624288e-09+/-5.622858e-09,
               -1.859026e-09+/-1.858552e-09, -3.250635e-12+/-3.249808e-12,
               -5.068914e-09+/-5.067634e-09, -4.871941e-09+/-4.87071e-09, 0.0+/-0,
               0.0+/-0, -6.871712e-09+/-4.952527e-09,
               -1.545212e-09+/-1.544828e-09, -6.40573e-10+/-6.404087e-10, 0.0+/-0,
               -4.935283e-09+/-2.901975e-09, -1.843172e-08+/-1.251576e-08,
               -6.71445e-09+/-6.71273e-09, -1.388548e-07+/-4.274043e-08,
               -4.241146e-08+/-1.762258e-08, -2.999543e-08+/-2.998778e-08,
               -1.520626e-07+/-5.006299e-07, -3.266883e-07+/-1.449028e-07,
               7.846369e-07+/-8.803142e-07, -1.519928e-07+/-1.07218e-07,
               -4.081559e-07+/-1.873273e-07, -3.434015e-08+/-1.027987e-06,
               -2.516069e-07+/-1.73851e-07, -8.333603e-07+/-2.853665e-07,
               1.439084e-06+/-1.77756e-06, -2.836471e-06+/-6.774013e-07,
               -7.692097e-07+/-2.280938e-07, -3.207971e-06+/-5.799841e-06,
               -8.805729e-06+/-3.719169e-06, -3.810922e-06+/-3.29137e-06,
               -1.671083e-05+/-5.398771e-06, -8.131938e-06+/-8.703075e-06,
               -1.213469e-05+/-6.712768e-06, 2.241294e-05+/-1.816346e-05,
               -3.46584e-05+/-2.47274e-05, -2.504294e-05+/-1.322211e-05,
               2.988232e-06+/-9.999848e-06, 9.248787e-06+/-2.051111e-05,
               8.657588e-06+/-2.534945e-05, -9.951473e-06+/-9.817888e-06,
               3.696361e-06+/-2.029254e-05, 8.702776e-06+/-1.343905e-05,
               1.237458e-05+/-3.052268e-05, 9.421495e-05+/-4.640119e-05,
               -3.918389e-05+/-3.697521e-05, 0.0001051884+/-6.277818e-05,
               0.0001998479+/-6.489579e-05, 0.0001497279+/-5.563772e-05,
               0.0001947749+/-5.449071e-05, 7.17222e-05+/-2.914425e-05,
               0.000134732+/-2.9629e-05, 0.0001109158+/-3.34169e-05,
               0.0001089335+/-2.818967e-05, 0.0003405998+/-4.422492e-05,
               0.0001581347+/-2.790862e-05, 0.0001519637+/-2.893234e-05,
               0.0002925294+/-4.266703e-05, 3.390917e-05+/-1.328072e-05,
               0.0003311925+/-4.068107e-05, 0.0003964351+/-3.973347e-05,
               0.0002032989+/-2.820316e-05, 4.749484e-05+/-1.391542e-05,
               7.308714e-05+/-1.919152e-05, 8.477883e-05+/-1.880586e-05,
               0.0003783144+/-3.910496e-05, 0.0002999691+/-3.663134e-05,
               0.0002810895+/-3.531906e-05, 0.0001691159+/-2.466477e-05,
               0.0001964518+/-2.750338e-05, 0.0001036065+/-1.975991e-05,
               7.118574e-05+/-1.79444e-05, 0.0003142769+/-3.554042e-05,
               0.0006935364+/-5.715178e-05, 0.0008185451+/-5.864059e-05,
               0.0001756991+/-2.74122e-05, 0.0005926753+/-4.821636e-05,
               0.001032059+/-6.463625e-05, 0.0001732703+/-2.667443e-05,
               0.0003230355+/-3.82033e-05, 0.0001058208+/-2.194266e-05,
               4.386201e-05+/-1.377494e-05, 2.231817e-05+/-7.640867e-06,
               3.150126e-06+/-3.827818e-06, 9.240951e-08+/-1.080839e-07, 0.0+/-0,
               0.0+/-0, 0.0+/-0], dtype=object)}
        >>> reader.energy_boundaries
        array([2.000e+07, 1.733e+07, 1.568e+07, 1.455e+07, 1.384e+07, 1.284e+07,
               1.000e+07, 8.187e+06, 6.434e+06, 4.800e+06, 4.304e+06, 3.000e+06,
               2.479e+06, 2.354e+06, 1.850e+06, 1.500e+06, 1.400e+06, 1.356e+06,
               1.317e+06, 1.250e+06, 1.200e+06, 1.100e+06, 1.010e+06, 9.200e+05,
               9.000e+05, 8.750e+05, 8.611e+05, 8.200e+05, 7.500e+05, 6.790e+05,
               6.700e+05, 6.000e+05, 5.730e+05, 5.500e+05, 4.920e+05, 4.700e+05,
               4.400e+05, 4.200e+05, 4.000e+05, 3.300e+05, 2.700e+05, 2.000e+05,
               1.490e+05, 1.283e+05, 1.000e+05, 8.500e+04, 8.200e+04, 7.500e+04,
               7.300e+04, 6.000e+04, 5.200e+04, 5.000e+04, 4.500e+04, 3.000e+04,
               2.000e+04, 1.700e+04, 1.300e+04, 9.500e+03, 8.030e+03, 5.700e+03,
               3.900e+03, 3.740e+03, 3.000e+03, 2.500e+03, 2.250e+03, 2.200e+03,
               1.800e+03, 1.550e+03, 1.500e+03, 1.150e+03, 9.500e+02, 6.830e+02,
               6.700e+02, 5.500e+02, 3.050e+02, 2.850e+02, 2.400e+02, 2.200e+02,
               2.095e+02, 2.074e+02, 2.020e+02, 1.930e+02, 1.915e+02, 1.885e+02,
               1.877e+02, 1.800e+02, 1.700e+02, 1.430e+02, 1.220e+02, 1.190e+02,
               1.175e+02, 1.160e+02, 1.130e+02, 1.080e+02, 1.050e+02, 1.012e+02,
               9.700e+01, 9.000e+01, 8.170e+01, 8.000e+01, 7.600e+01, 7.200e+01,
               6.750e+01, 6.500e+01, 6.300e+01, 6.100e+01, 5.800e+01, 5.340e+01,
               5.060e+01, 4.830e+01, 4.520e+01, 4.400e+01, 4.240e+01, 4.100e+01,
               3.960e+01, 3.910e+01, 3.800e+01, 3.763e+01, 3.727e+01, 3.713e+01,
               3.700e+01, 3.600e+01, 3.550e+01, 3.500e+01, 3.375e+01, 3.325e+01,
               3.175e+01, 3.125e+01, 3.000e+01, 2.750e+01, 2.500e+01, 2.250e+01,
               2.175e+01, 2.120e+01, 2.050e+01, 2.000e+01, 1.940e+01, 1.850e+01,
               1.700e+01, 1.600e+01, 1.440e+01, 1.290e+01, 1.190e+01, 1.150e+01,
               1.000e+01, 9.100e+00, 8.100e+00, 7.150e+00, 7.000e+00, 6.875e+00,
               6.750e+00, 6.500e+00, 6.250e+00, 6.000e+00, 5.400e+00, 5.000e+00,
               4.700e+00, 4.100e+00, 3.730e+00, 3.500e+00, 3.200e+00, 3.100e+00,
               3.000e+00, 2.970e+00, 2.870e+00, 2.770e+00, 2.670e+00, 2.570e+00,
               2.470e+00, 2.380e+00, 2.300e+00, 2.210e+00, 2.120e+00, 2.000e+00,
               1.940e+00, 1.860e+00, 1.770e+00, 1.680e+00, 1.590e+00, 1.500e+00,
               1.450e+00, 1.400e+00, 1.350e+00, 1.300e+00, 1.250e+00, 1.225e+00,
               1.200e+00, 1.175e+00, 1.150e+00, 1.140e+00, 1.130e+00, 1.120e+00,
               1.110e+00, 1.100e+00, 1.090e+00, 1.080e+00, 1.070e+00, 1.060e+00,
               1.050e+00, 1.040e+00, 1.030e+00, 1.020e+00, 1.010e+00, 1.000e+00,
               9.750e-01, 9.500e-01, 9.250e-01, 9.000e-01, 8.500e-01, 8.000e-01,
               7.500e-01, 7.000e-01, 6.500e-01, 6.250e-01, 6.000e-01, 5.500e-01,
               5.000e-01, 4.500e-01, 4.000e-01, 3.750e-01, 3.500e-01, 3.250e-01,
               3.000e-01, 2.750e-01, 2.500e-01, 2.250e-01, 2.000e-01, 1.750e-01,
               1.500e-01, 1.250e-01, 1.000e-01, 9.000e-02, 8.000e-02, 7.000e-02,
               6.000e-02, 5.000e-02, 4.000e-02, 3.000e-02, 2.530e-02, 1.000e-02,
               7.500e-03, 5.000e-03, 4.000e-03, 3.000e-03, 2.500e-03, 2.000e-03,
               1.500e-03, 1.200e-03, 1.000e-03, 7.500e-04, 5.000e-04, 1.000e-04,
               1.000e-05])
        """
        super().__init__(filename)
        
        # Now only return the region integrated sdf profiles
        # i.e. those with zone number and zone volume both equal to 0
        self.filename = filename
        self.sdf_data = [ match for match in self.sdf_data if match['zone_number'] == 0 and match['zone_volume'] == 0 ]
    
    def __repr__(self):
        return f"<tsunami_ip_utils.readers.RegionIntegratedSdfReader object from {self.filename}>"

    def convert_to_dict(self, key: str='names'):
        """Converts the sdf data into a dictionary keyed by nuclide-reaction pair or by ZAID and reaction MT.
        
        Parameters
        ----------
        key
            The key to use for the dictionary. Default is ``'names'`` which uses the isotope name and reaction type.
            If ``'numbers'`` is supplied instead then the ZAID and reaction MT are used.
            
        Examples
        --------
        >>> reader.convert_to_dict()
        <tsunami_ip_utils.readers.RegionIntegratedSdfReader object from tests/example_files/sphere_model_1.sdf>
        >>> reader.sdf_data.keys()
        dict_keys(['u-234', 'u-235', 'u-238'])
        >>> reader.sdf_data['u-234'].keys()
        dict_keys(['total', 'elastic', "n,n'", 'n,2n', 'fission', 'capture', 'n,gamma', 'chi', 'nubar'])
        >>> reader.sdf_data['u-234']['total']
        {'isotope': 'u-234', 'reaction_type': 'total', 'zaid': '92234', 'reaction_mt': '1', 'zone_number': 0, 'zone_volume': 0, 'energy_integrated_sensitivity': 0.008984201+/-0.0002456359, 'abs_sum_groupwise_sensitivities': 0.009319556, 'sum_opposite_sign_groupwise_sensitivities': -0.0001676775+/-4.959921e-05, 'sensitivities': array([0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -4.393147e-10+/-4.39204e-10, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -1.0871e-09+/-1.086821e-09, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, -8.666447e-10+/-8.664239e-10, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, -2.021382e-09+/-2.020862e-09, 0.0+/-0, 0.0+/-0,
               0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -2.351273e-11+/-2.350687e-11, 1.707797e-08+/-1.707371e-08,
               -7.172835e-10+/-7.171047e-10, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               -1.89505e-09+/-1.894573e-09, 0.0+/-0, 0.0+/-0, 0.0+/-0, 0.0+/-0,
               0.0+/-0, -5.624288e-09+/-5.622858e-09,
               -1.859026e-09+/-1.858552e-09, -3.250635e-12+/-3.249808e-12,
               -5.068914e-09+/-5.067634e-09, -4.871941e-09+/-4.87071e-09, 0.0+/-0,
               0.0+/-0, -6.871712e-09+/-4.952527e-09,
               -1.545212e-09+/-1.544828e-09, -6.40573e-10+/-6.404087e-10, 0.0+/-0,
               -4.935283e-09+/-2.901975e-09, -1.843172e-08+/-1.251576e-08,
               -6.71445e-09+/-6.71273e-09, -1.388548e-07+/-4.274043e-08,
               -4.241146e-08+/-1.762258e-08, -2.999543e-08+/-2.998778e-08,
               -1.520626e-07+/-5.006299e-07, -3.266883e-07+/-1.449028e-07,
               7.846369e-07+/-8.803142e-07, -1.519928e-07+/-1.07218e-07,
               -4.081559e-07+/-1.873273e-07, -3.434015e-08+/-1.027987e-06,
               -2.516069e-07+/-1.73851e-07, -8.333603e-07+/-2.853665e-07,
               1.439084e-06+/-1.77756e-06, -2.836471e-06+/-6.774013e-07,
               -7.692097e-07+/-2.280938e-07, -3.207971e-06+/-5.799841e-06,
               -8.805729e-06+/-3.719169e-06, -3.810922e-06+/-3.29137e-06,
               -1.671083e-05+/-5.398771e-06, -8.131938e-06+/-8.703075e-06,
               -1.213469e-05+/-6.712768e-06, 2.241294e-05+/-1.816346e-05,
               -3.46584e-05+/-2.47274e-05, -2.504294e-05+/-1.322211e-05,
               2.988232e-06+/-9.999848e-06, 9.248787e-06+/-2.051111e-05,
               8.657588e-06+/-2.534945e-05, -9.951473e-06+/-9.817888e-06,
               3.696361e-06+/-2.029254e-05, 8.702776e-06+/-1.343905e-05,
               1.237458e-05+/-3.052268e-05, 9.421495e-05+/-4.640119e-05,
               -3.918389e-05+/-3.697521e-05, 0.0001051884+/-6.277818e-05,
               0.0001998479+/-6.489579e-05, 0.0001497279+/-5.563772e-05,
               0.0001947749+/-5.449071e-05, 7.17222e-05+/-2.914425e-05,
               0.000134732+/-2.9629e-05, 0.0001109158+/-3.34169e-05,
               0.0001089335+/-2.818967e-05, 0.0003405998+/-4.422492e-05,
               0.0001581347+/-2.790862e-05, 0.0001519637+/-2.893234e-05,
               0.0002925294+/-4.266703e-05, 3.390917e-05+/-1.328072e-05,
               0.0003311925+/-4.068107e-05, 0.0003964351+/-3.973347e-05,
               0.0002032989+/-2.820316e-05, 4.749484e-05+/-1.391542e-05,
               7.308714e-05+/-1.919152e-05, 8.477883e-05+/-1.880586e-05,
               0.0003783144+/-3.910496e-05, 0.0002999691+/-3.663134e-05,
               0.0002810895+/-3.531906e-05, 0.0001691159+/-2.466477e-05,
               0.0001964518+/-2.750338e-05, 0.0001036065+/-1.975991e-05,
               7.118574e-05+/-1.79444e-05, 0.0003142769+/-3.554042e-05,
               0.0006935364+/-5.715178e-05, 0.0008185451+/-5.864059e-05,
               0.0001756991+/-2.74122e-05, 0.0005926753+/-4.821636e-05,
               0.001032059+/-6.463625e-05, 0.0001732703+/-2.667443e-05,
               0.0003230355+/-3.82033e-05, 0.0001058208+/-2.194266e-05,
               4.386201e-05+/-1.377494e-05, 2.231817e-05+/-7.640867e-06,
               3.150126e-06+/-3.827818e-06, 9.240951e-08+/-1.080839e-07, 0.0+/-0,
               0.0+/-0, 0.0+/-0], dtype=object)}"""
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

def read_uncertainty_contributions_out(filename: Union[str, Path]) -> Tuple[List[dict], List[dict]]:
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
    element = Word(alphas.lower(), max=2) 
    isotope_name = Combine(element + Optional('-' + atomic_number)) # To handle the case of carbon in ENDF-7.1 libraries
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

def read_uncertainty_contributions_sdf(filenames: List[Path]):
    """Reads the uncertainty contributions from a list of TSUNAMI-B SDF files and returns the contributions for each nuclide-
    reaction covariance by first running a TSUNAMI-IP calculation to generate the extended uncertainty edit.
    
    Parameters
    ----------
    filenames
        List of paths to the SDF files.
        
    Returns
    -------"""

    # ===============================
    # Generate the Uncertainty Edits
    # ===============================

    # First generate the extended uncertainty edits
    current_dir = Path(__file__).parent
    with open(current_dir / "input_files" / "tsunami_ip_uncertainty_contributions.inp", 'r') as f:
        tsunami_ip_template = Template(f.read())

    # Generate a string containing a list of all filenames
    filenames = [ str(filename.absolute()) for filename in filenames ]

    # Now template the input file
    tsunami_ip_input = tsunami_ip_template.substitute(
        filenames='\n'.join(filenames),
        first_file=filenames[0]
    )

    # Now write the template file to a temporary file and run it
    with NamedTemporaryFile('w', delete=False) as f:
        f.write(tsunami_ip_input)
        input_filename = f.name

    # Run the TSUNAMI-IP calculation
    process = subprocess.Popen( ['scalerte', input_filename], cwd=str( Path( input_filename ).parent ) )
    process.wait()

    # ========================
    # Process the Output File
    # ========================
    with open(f"{input_filename}.out", 'r') as f:
        data = f.read()
    
    # ----------------------------------------------------
    # Define the formattting that precedes the data table
    # ----------------------------------------------------
    table_identifier = Literal("contributions to uncertainty in keff ( % dk/k ) by individual energy covariance matrices:")
    pre_header = Literal("covariance matrix") + LineEnd()
    header = Word("nuclide-reaction") + Word("with") + Word("nuclide-reaction") + Word("% delta-k/k due to this matrix")
    dash_separator = OneOrMore(OneOrMore('-'))

    # ----------------------
    # Define the data lines
    # ----------------------

    # Define the grammar for the nuclide-reaction pair
    atomic_number = Word(nums)
    element = Word(alphas.lower(), max=2) 
    isotope_name = Combine(element + Optional('-' + atomic_number)) # To handle the case of carbon in ENDF-7.1 libraries
    reaction_type = Word(alphanums + ',\'')

    data_line = Group(isotope_name + reaction_type + isotope_name + reaction_type + \
                pyparsing_common.sci_real + Suppress(Literal("+/-")) + pyparsing_common.sci_real + Suppress(LineEnd()))
    data_block = OneOrMore(data_line)

    # -------------------------------------------
    # Define the total parser and parse the data
    # -------------------------------------------
    data_parser = Suppress(table_identifier) + Suppress(pre_header) + Suppress(header) + \
                    Suppress(dash_separator) + data_block

    # -------------------------------------------------------------------------------
    # Now convert the data into isotope wise and isotope-reaction wise contributions
    # -------------------------------------------------------------------------------
    parsed_data = data_parser.searchString(data)
    num_sdfs = len(parsed_data)
    isotope_reaction = [ [] for _ in range(num_sdfs) ]
    for i, match in enumerate(parsed_data):
        for data_element in match:
            isotope_reaction[i].append({
                'isotope': f'{data_element[0]} - {data_element[2]}',
                'reaction_type': f'{data_element[1]} - {data_element[3]}',
                'contribution': ufloat(data_element[4], data_element[5])
            })

    # Now calculate nuclide totals by summing the contributions for each nuclide via total = sqrt((pos)^2 - (neg)^2)
    isotope_totals = [ {} for _ in range(num_sdfs) ]
    for i, application in enumerate(isotope_reaction):
        for data in application:
            # First add up squared sums of all reaction-wise contributions
            isotope = data['isotope']
            contribution = data['contribution']
            if isotope not in isotope_totals[i].keys():
                isotope_totals[i][isotope] = ufloat(0,0)

            if contribution < 0:
                isotope_totals[i][isotope] -= ( data['contribution'] )**2
            else:
                isotope_totals[i][isotope] += ( data['contribution'] )**2
        
        # Now take square root of all contributions
        for isotope, total in isotope_totals[i].items():
            isotope_totals[i][isotope] = total**0.5 if total > 0 else -(-total)**0.5 # Allow for negative isotope totals

        # Now convert into a list of dictionaries
        isotope_totals[i] = [ {'isotope': isotope, 'contribution': total} for isotope, total in isotope_totals[i].items() ]

    # Remove the temporary input and output files
    os.remove(input_filename)
    os.remove(f"{input_filename}.out")

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
        Integral matrices for each integral index type. The shape of the matrices are ``(num_applications, num_experiments)``. 
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
        "C_k":       C_k.transpose(),
        "E_total":   E_total.transpose(),
        "E_fission": E_fission.transpose(),
        "E_capture": E_capture.transpose(),
        "E_scatter": E_scatter.transpose()
    })

    return integral_matrices

def read_region_integrated_h5_sdf(filename: Path) -> Dict[str, unumpy.uarray]:
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