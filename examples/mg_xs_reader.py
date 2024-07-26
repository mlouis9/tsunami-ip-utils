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
# This function is not the most user-friendly, and requires teh user to input a nuclide_reaction_dict in terms of nuclide
# ZAID and reaction MT numbers. This function is primarily used by the :mod:`tsunami_ip_utils.perturbations` module for
# reading cross sections for perturbation calculations.