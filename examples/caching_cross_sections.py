"""
Cache Base and Perturbed Cross Section Libraries
================================================
This demonstrates how to cache the base and perturbed cross section libraries for a given SCALE multigroup library. This is useful for
preventing the need to regenerate the perturbed cross section libraries each time they are needed.
"""
from tsunami_ip_utils.perturbations import cache_all_libraries
from paths import EXAMPLES
from pathlib import Path

multigroup_library = EXAMPLES / 'data' / 'dummy_56_v7.1'
perturbation_factors = Path("~/codes/SCALE-6.3.1/data/perturb/56n.v7.1")
cache_all_libraries(multigroup_library, perturbation_factors)