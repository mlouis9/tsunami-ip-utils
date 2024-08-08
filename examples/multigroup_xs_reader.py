"""
SCALE Multigroup Cross Section Reader
=====================================

This example demonstrates how to read SCALE multigroup cross section libraries. This functionality relies on AMPX modules
for creating text dumps of COVERX formatted libraries. More details are given here: :ref:`mg_reader.inp`. 
"""

# %%
# A Disclaimer
# ------------
# The COVERX formatted cross section library referenced in this example ``dummy_56_v7.1`` was created with AMPX using ExSite
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
multigroup_library_path = EXAMPLES / 'data' / 'dummy_56_v7.1'
out = read_multigroup_xs(multigroup_library_path, nuclide_reaction_dict)
print(len(out['92235']['1']))

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

with open(EXAMPLES / 'data' / 'dummy_56_v7.1.pkl', 'wb') as f:
    pickle.dump(out, f)

# %%
# Getting the Energy Structure
# ----------------------------
# By default the reader doesn't return the energy group structure (this is something that could be improved in the future), but if
# you're examining SCALE built-in cross section libraries, the energy group structure can be obtained from the 
# :func:`tsunami_ip_utils.xs.get_scale_multigroup_structure` function. This function scrapes the SCALE manual for the energy
# group structure of a given SCALE library.

from tsunami_ip_utils.xs import get_scale_multigroup_structure

num_groups = 56
multigroup_structure = get_scale_multigroup_structure(num_groups)

print(multigroup_structure)

# %%
# We can easily plot cross section data using the energy group structure. The following code snippet demonstrates how to plot
# the u-235 fission cross section

import matplotlib.pyplot as plt
import numpy as np

u235_fission_xs = out['92235']['18']

# To get the cross sections to be visually constant within each energy group, we need to repeat the cross sections and add in the
# implied lower energy bound of 1E-05 eV
u235_fission_xs = u235_fission_xs.repeat(2)
modified_energies = np.zeros_like(u235_fission_xs)

energies = multigroup_structure[:, 1]
modified_energies[0] = 1E-05
modified_energies[-1] = energies[-1]
modified_energies[1:] = np.repeat(energies, 2)[:-1]

plt.plot(modified_energies, u235_fission_xs)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy (eV)')
plt.ylabel('Cross Section (b)')
plt.grid()
plt.show()

# %%
# A Future Improvement
# ---------------------
# It is unfortunate that to cache a cross section library, the library must be read twice (the most consuming part is making
# the text dump again). This is a limitation of the current implementation, and may be improved in the future. To avoid this,
# the text dump just needs to be saved so that it can be read by the second function call. This could all be implemented by
# adding an additional flag to the :func:`tsunami_ip_utils.xs.read_multigroup_xs` that does this under the hood.