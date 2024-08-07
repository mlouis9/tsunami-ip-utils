"""
Reading Integral Indices
========================
This example demonstrates how to read/get TSUNAMI-IP integral indices (like :math:`c_k` and :math:`E`) computed from a set of 
experiments and applications.
"""

# %%
# Reading Integral Indices From a TSUNAMI-IP ``.out`` File
# --------------------------------------------------------
# The first way to obtain integral indices is to first create a TSUNAMI-IP input file (see `the SCALE manual 
# <https://scale-manual.ornl.gov/tsunami-ip.html#user-input>`_), run it with SCALE, and then read the integral indices from the
# resulting ``.out`` file. There's a built in function to do this in ``tsunami_ip_utils``, provided the TSUNAMI-IP input file
# follows the :func:`expected format <tsunami_ip_utils.readers.read_integral_indices>`: 
# :func:`tsnami_ip_utils.readers.read_integral_indices`. This function just takes a path to the ``.out`` file and returns a
# dictionary of integral index matrices. Here's an example of how to use it:

from tsunami_ip_utils.readers import read_integral_indices
from paths import EXAMPLES

# Path to the .out file
out_file = EXAMPLES / 'data' / 'tsunami_ip_hmf.out'

# Read the integral indices
integral_indices = read_integral_indices(out_file)

# Print the integral indices
for key, value in integral_indices.items():
    print(f'{key}:\n{value}\n')

# %%
# From the output, we can see that the integral indices are stored in a dictionary where the keys are the names of the integral
# indices and the values are the matrices of the integral indices. The integral indices are stored as numpy arrays, and have shapes
# ``(num_applications, num_experiments)``.

# %%
# Getting Integral Indices From a Set of Application and Experiment SDFs
# -----------------------------------------------------------------------
# Another way to get the integral indices from a given set of applications and experiments that requires less manual effort is to use
# the :func:`tsunami_ip_utils.integral_indices.get_integral_indices` function. This function takes a list of application SDFs and
# experiment SDFs, and the covariance library (which the user would have to manually input when creating their TSUNAMI-IP input
# file in the previous section anyways). This function just creates a TSUNAMI-IP input file, runs it, then extracts the results, simply
# abstracting away the manual effort required in the previous method.

from tsunami_ip_utils.integral_indices import get_integral_indices

application_sdfs = [ EXAMPLES / 'data' / 'example_sdfs' / 'HMF' / f'HEU-MET-FAST-003-00{i}.sdf' for i in range(1, 13) ]
experiment_sdfs = application_sdfs
coverx_library = '252groupcov7.1'
integral_indices = get_integral_indices(application_sdfs, experiment_sdfs, coverx_library=coverx_library)

# Print the integral indices
for key, value in integral_indices.items():
    print(f'{key}:\n{value}\n')

# %%
# Spaces in Annotated Names
# -------------------------
# In each SDF file, the first line contains an annotated name for the SDF file, that is used to identify the SDF in the TSUNAMI-IP
# output tables. If these names contain spaces (as some of those provided in the MCT directory do), they can cause the
# :func:`tsnami_ip_utils.readers.read_integral_indices` function to fail. For example, the following code will raise an error:

try:
    integral_indices = read_integral_indices(EXAMPLES / 'data' / 'tsunami_ip_mct.out')
    print(integral_indices)
except Exception as e:
    print(e)

# %%
# This is because of the spaces in the following SDF annotated names (that were used to generate the TSUNAMI-IP output file):

sdf_paths = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-002-00{i}S.sdf' for i in range(1, 7) ]

for sdf_path in sdf_paths:
    with open(sdf_path, 'r') as f:
        print(f.readline())

# %%
# So if you have SDFs with names that contain spaces, you should either use the 
# :func:`tsunami_ip_utils.integral_indices.get_integral_indices`, which will handle these names by replacing whitespace with ``'-'``,
# or you can manually replace the whitespaces in the SDF annotated names using the :func:`tsunami_ip_utils.utils.modify_sdf_names`
# function. Below is an example showing the correct output using the :func:`tsunami_ip_utils.integral_indices.get_integral_indices`
# function:

from tsunami_ip_utils.integral_indices import get_integral_indices

application_sdfs = [ EXAMPLES / 'data' / 'example_sdfs' / 'MCT' / f'MIX-COMP-THERM-002-00{i}S.sdf' for i in range(1, 7) ]
experiment_sdfs = application_sdfs
coverx_library = '252groupcov7.1'
integral_indices = get_integral_indices(application_sdfs, experiment_sdfs, coverx_library)

# Print the integral indices
for key, value in integral_indices.items():
    print(f'{key}:\n{value}\n')