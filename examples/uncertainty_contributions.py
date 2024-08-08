"""
Uncertainty/Variance Contributions
==================================
This example demonstrates how to read (nuclear-data induced) uncertainty/variance contributions from a set of TSUNAMI-B SDF files.
"""

# Getting the Contributions
# -------------------------
# The first step is to get the contributions to the uncertainty. This can be done using the 
# :func:`tsunami_ip_utils.integral_indices.get_uncertainty_contributions` function.

from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from paths import EXAMPLES

application_filenames = [ f"{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-00{i}.sdf" for i in range(1, 3) ]
experiment_filenames  = [ f"{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-00{i}.sdf" for i in range(3, 6) ]
uncertainty_contributions_nuclide, uncertainty_contributions_nuclide_reaction = \
      get_uncertainty_contributions(application_filenames, experiment_filenames)
# %%
# As explained in the API documentation, the above function returns two dictionaries. The first dictionary,
# ``uncertainty_contributions_nuclide``, contains the contributions to the uncertainty for each nuclide. The second dictionary,
# ``uncertainty_contributions_nuclide_reaction``, contains the contributions to the uncertainty for each nuclide-reaction pair.
# Each dictionary has keys ``'application'`` and ``'experiment'``.

print(uncertainty_contributions_nuclide.keys())
application_nuclide_contributions = uncertainty_contributions_nuclide['application']

# %%
# Isotope-Wise Contributions
# --------------------------
# The contents of the dictionary is a list of contributions for each application/experiment. I.e. the (nuclide-wise) 
# contributions to the nuclear data induced uncertainty for application 1 (i.e. ``application_filenames[0]``) can be 
# accessed as follows:

print(application_nuclide_contributions[1])

# %%
# The ouput is a list of dictionaries with the keys ``'isotope'`` and ``'contribution'``. The contributions are the nuclide-wise
# contributions to the nuclear data induced uncertainty. For TSUNAMI-B formatted SDFs, these values have an associated
# (monte carlo) uncertainty, and so are represented as :func:`uncertainties.ufloat` objects.

print(application_nuclide_contributions[1][0]['contribution'])

# %%
# These objects automatically handle uncertainty propagation, and there use is further documented in the :mod:`uncertainties`
# package.

# %%
# Isotope-Reaction-Wise Contributions
# -----------------------------------
# Isotope-reaction-wise contributions are accessed similarly.

application_nuclide_contributions = uncertainty_contributions_nuclide_reaction['application']
print(application_nuclide_contributions[1])

# %%
# The output is a list of dictionaries with the keys ``'isotope'``, ``'reaction'``, and ``'contribution'``. The contributions
# are (like before) :func:`uncertainties.ufloat` objects. A specific contribution can be accessed via

print(application_nuclide_contributions[1][0])
print(application_nuclide_contributions[1][0]['contribution'])

# %%
# Getting Contributions for a Single File
# ---------------------------------------
# As discussed in the API documentation, the application_filenames, and experiment_filenames arugments are both optional - to
# streamline the process of getting contributions for a single file, the function can be called with only one argument. For
# example, to get the contributions for the first application file, we can call the function as follows:

application_filenames = [ f"{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-001.sdf" ]
uncertainty_contributions_nuclide, uncertainty_contributions_nuclide_reaction = \
      get_uncertainty_contributions(application_filenames)
print(uncertainty_contributions_nuclide['application'])

# %%
# Note that since we only passed application filenames, the function only returns the contributions for the application. The
# experiment contributions are an empty list

print(uncertainty_contributions_nuclide['experiment'])

# %%
# Variance Contributions
# ----------------------
# The contributions to the nuclear data induced variance can also be obtained by passing the ``variance=True`` argument to the
# function. This will return the contributions to the variance, rather than the uncertainty.

variance_contributions_nuclide, variance_contributions_nuclide_reaction = \
      get_uncertainty_contributions(application_filenames, variance=True)

print(variance_contributions_nuclide['application'])

# %%
# The variances are just the squares of the uncertainties, except when the are negative, in which case the variance is the
# negative of the square of the uncertainty.

variance_contribution = variance_contributions_nuclide['application'][0][0]['contribution']
uncertainty_contribution = uncertainty_contributions_nuclide['application'][0][0]['contribution']
squared_uncertainty = uncertainty_contribution**2 if uncertainty_contribution >= 0 else  -(uncertainty_contribution)**2

print(variance_contribution, squared_uncertainty)
assert ( variance_contribution.n == squared_uncertainty.n ) and ( variance_contribution.s == squared_uncertainty.s )

# %%
# Note that even though the variance contribution and the squared uncertainty contribution have the same nominal value
# and standard deviation, they are not equal, i.e. ``variance_contribution == squared_uncertainty`` evaluates to ``False``. 
# This is because the uncertainties are represented as :func:`uncertainties.ufloat` which represent distinct random variables, 
# as discussed in `the uncertainties package documentation <https://pythonhosted.org/uncertainties/user_guide.html#comparison-operators>`_.