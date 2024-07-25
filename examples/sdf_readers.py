"""
SDF Readers
===========

This demonstrates how to use the tools for reading `TSUNAMI-B <https://scale-manual.ornl.gov/tsunami-ip-appAB.html#format-of-tsunami-b-sensitivity-data-file>`_
formatted Sensitivity Data Files (SDFs).
"""

# %%
# Reading .txt formatted .sdf files
# ---------------------------------
# 
# Reading all SDF Profiles
# ^^^^^^^^^^^^^^^^^^^^^^^^
# There is a general class (:class:`tsunami_ip_utils.readers.SdfReader`) which reads all SDF profiles 
# (including region-wise and region-integrated) profiles.

from tsunami_ip_utils.readers import RegionIntegratedSdfReader, SdfReader
from paths import EXAMPLES

data = SdfReader(f'{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-001.sdf')

# %%
# This class has two attributes, :attr:`tsunami_ip_utils.readers.SdfReader.sdf_data` and 
# :attr:`tsunami_ip_utils.readers.SdfReader.energy_boundaries` which represent the sdf data and the energy boundaries of the 
# multigroup structure (used to tally the senstivities), respectively. To access the data, use the ``sdf_data`` attribute.

print(data.sdf_data)

# %%
# The data is stored as a list of dictionaries representing each entry in the SDF file. The dictionary has keys corresponding
# to the isotope, reaction, region number, integral value for the sdf profile (i.e. the energy-integrated sensitivity), and
# a key ``sensitivities`` which contains the full energy-dependent sensitivity profile stored as a :func:`uncertainties.unumpy.uarray`
# objects (essentially just numpy arrays with :func:`uncertainties.ufloat` objects). Detailed documentation about the keys
# of the dicionaries can be found in the documentation for the :attr:`tsunami_ip_utils.readers.SdfReader.sdf_data` attribute.
# 
# For example, the sensitivity profile for the first entry can be accessed by

print(data.sdf_data[0]['isotope'], data.sdf_data[0]['reaction_type'])
print(data.sdf_data[0]['sensitivities'])

# %%
# Reading Region-Integrated SDF profiles
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# There is a separate class (:class:`tsunami_ip_utils.readers.RegionIntegratedSdfReader`) which reads only the region-integrated
# SDF profiles. Note that this data could be obtained using :class:`tsunami_ip_utils.readers.SdfReader` by simply filtering
# the data for the region-integrated profiles. However, this class is provided for convenience, and the fact that only region
# integrated data exists guarantees that there is one and only one SDF profile for each nuclide-reaction pair, which allows
# the data to be stored as a nested dictionary. Sdf data is retreived in the same way as the general SDF reader

data = RegionIntegratedSdfReader(f'{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-001.sdf')
print(data.sdf_data)

# %%
# The data can be converted to a nested dictionary by calling the 
# :meth:`tsunami_ip_utils.readers.RegionIntegratedSdfReader.convert_to_dict` metod

data.convert_to_dict()
print(data.sdf_data)

# %%
# Reading .h5 formatted .sdf files
# ---------------------------------
# SCALE provides utilities (like `tao <https://scale-manual.ornl.gov/tsunami-ip-appAB.html#format-of-hdf5-based-sensitivity-data-file>`_)
# for converting .sdf files to .h5 files. These files can be read easily using `h5py <https://docs.h5py.org/en/stable>`_, but
# the function :func:`tsunami_ip_utils.readers.read_region_integrated_h5_sdf` provides an interface that outputs a dictionary
# of a similar format as that produced by the :meth:`tsunami_ip_utils.readers.RegionIntegratedSdfReader.convert_to_dict` method,
# except that only the sensitivity data is included, and the dictionary is keyed by nuclide zaid and reaction mt rather than
# the names.

from tsunami_ip_utils.readers import read_region_integrated_h5_sdf

data = read_region_integrated_h5_sdf(f'{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-001.sdf.h5')
print(data)