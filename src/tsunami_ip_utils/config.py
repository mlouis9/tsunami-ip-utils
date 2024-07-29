"""A package level configuration module."""
from pathlib import Path

NUM_SAMPLES = 1000
"""Number of cross section perturbation factors available in SCALE."""

SDF_DATA_NAMES = [
    "isotope",
    "reaction_type",
    "zaid",
    "reaction_mt",
    "zone_number",
    "zone_volume",
    "energy_integrated_sensitivity",
    "abs_sum_groupwise_sensitivities",
    "sum_opposite_sign_groupwise_sensitivities",
    "sensitivities",
    "uncertainties"
]
"""Names of the data fields parsed by the SDF reader for TSUNAMI-B formatted SDF files."""

COMPARISON_HEATMAP_LABELS = {
    'Calculated': 'Calculated Integral Index',
    'TSUNAMI-IP': 'TSUNAMI-IP Integral Index',
    'Percent Difference': 'Percent Difference (%)'
}
"""Labels for comparison heatmaps."""

generating_docs = False
"""Whether or not to kill interactive legend (flask/dash applications) plots after a short amount of time. This is
not intended for use by users, but is necessary for generating documentation properly."""

cache_dir = Path("~/.tsunami_ip_utils_cache").expanduser()
"""Directory to store cached cross section libraries and perturbations. This is also where the package will look for already cached
data, so be sure it corresponds with where your cached data actually is, if you have manually changed this."""