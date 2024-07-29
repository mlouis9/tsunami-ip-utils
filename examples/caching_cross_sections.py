"""
Cache Base and Perturbed Cross Section Libraries
================================================
This demonstrates how to cache the base and perturbed cross section libraries for a given SCALE multigroup library. This is useful for
preventing the need to regenerate the perturbed cross section libraries each time a perturbation correlation plot is created.
"""

# %%
# This module uses clarolplus along with the paths to the base cross section library and the perturbation factors to create the
# perturbed cross section libraries, then cache them as ``.pkl`` files which can be readily loaded into python for future use.

from tsunami_ip_utils.perturbations import cache_all_libraries
from paths import EXAMPLES
from pathlib import Path

multigroup_library = EXAMPLES / 'data' / 'dummy_56_v7.1'
perturbation_factors = Path("~/codes/SCALE-6.3.1/data/perturb/56n.v7.1")
cache_all_libraries(multigroup_library, perturbation_factors)

# %%
# Note that the reported runtime is not accurate, and actually caching all of the perturbations of the dummy 56 group library takes
# around 30 minutes (depending on the number of available cores, since this caching is multi-threaded by default). As mentioned in the
# API documentation, this can take up to ~5 hours for the 252 group library, and will consume a significant amount of disk space. Now
# to cache one of the SCALE libraries for production use, just change ``multigroup_library`` to be a path to one of SCALE's default
# libraries (e.g. ``Path("~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1")``), and change the ``perturbation_factors`` to the path to
# the perturbation factors for the library you are using (e.g. ``Path("~/codes/SCALE-6.3.1/data/perturb/252n.v7.1")``) and run the
# the lines of code above.