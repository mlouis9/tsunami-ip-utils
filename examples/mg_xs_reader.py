"""
SCALE Multigroup Cross Section Reader
=====================================

This example demonstrates how to read SCALE multigroup cross section libraries. This functionality relies on AMPX modules
for creating text dumps of COVERX formatted libraries. More details are given here: :ref:`mg_reader.inp`. 
"""

# %%
# A Disclaimer
# ------------
# The COVERX formatted cross section library referenced in this example ``dummy_52_v7.1`` was created with AMPX using ExSite
# and publicly available cross section data from `ENDF VII.1 <https://www.nndc.bnl.gov/endf-b7.1/download.html>`_.
# This library is **NOT** intended for use in simulations, and only contains a small subset of nuclides and reactions. This
# library is included purely for demonstration purposes, and when using functionality in this package that requires nuclear
# data, please use the provided SCALE libraries.

# %%
# Reading Multigroup Cross Sections
# ---------------------------------
# The following code snippet demonstrates how to read a SCALE multigroup cross section library. The function 
# :func:`tsunami_ip_utils.xs.read_multigroup_xs` is used to read COVER formatted SCALE cross section libraries, and takes
# a :class:`pathlib.Path` object to the multigroup library, and a dictionary of nuclide-reaction pairs to read. The function
# returns a nested dictionary of :class:`numpy.ndarray` objects, where the outer dictionary keys are nuclides, and the inner
# dictionary keys are reactions. The values are the multigroup cross sections for the corresponding nuclide-reaction pair. 

from tsunami_ip_utils.xs import read_multigroup_xs
from paths import EXAMPLES

nuclide_reaction_dict = {'92235': ['1', '18'], '5011': ['1', '27'], '94239': ['1', '18']}
multigroup_library_path = EXAMPLES / 'data' / 'dummy_52_v7.1'
out = read_multigroup_xs(multigroup_library_path, nuclide_reaction_dict)
print(out)

# %%
# This function is parallel, and reads cross section libraries on multiple cores, which can be useful for large libraries.
# This function isn't the most user-friendly, and requires the user to input a nuclide_reaction_dict in terms of nuclide
# ZAID and reaction MT numbers, but it is primarily used by the :mod:`tsunami_ip_utils.perturbations` module for
# reading cross sections for perturbation calculations.

# %%
# Caching a Multigroup Cross Section Library
# ------------------------------------------
# If working with python-based applications that deal with SCALE multigroup cross section data, it may be useful to cache the
# cross section libary in a convenient format for reading into python, like a ``.pkl`` file. To avoid having to manually supply
# the list of all nuclides and reactions in the multigroup library, the :func:`tsunami_ip_utils.xs.read_multigroup_xs` function
# can be run with the ``return_available_nuclide_reactions`` flag set to ``True``. This will return a dictionary of all nuclides
# and reactions in the library, which can then be used to read the library and cache it.

out, available_nuclide_reactions = read_multigroup_xs(
    multigroup_library_path, 
    nuclide_reaction_dict, 
    return_available_nuclide_reactions=True
)
print(available_nuclide_reactions)

# %%
# The list of all available nuclide reactions can then be supplied to the function to read the entire library.

out = read_multigroup_xs(multigroup_library_path, available_nuclide_reactions)

# %%
# And the cross sections can be easily cached

import pickle

with open(EXAMPLES / 'data' / 'dummy_52_v7.1.pkl', 'wb') as f:
    pickle.dump(out, f)

# Now compare the dump to the gold standard
import filecmp

assert filecmp.cmp(EXAMPLES / 'data' / 'dummy_52_v7.1.pkl', EXAMPLES / 'gold' / 'dummy_52_v7.1.pkl')