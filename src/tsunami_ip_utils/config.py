"""A package level configuration module."""

# config.py

# Number of cross section perturbation factors available in SCALE
NUM_SAMPLES = 1000

# Names of the data fields parsed by the SDF reader for TSUNAMI-B formatted SDF files
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

# Labels for comparison heatmaps
COMPARISON_HEATMAP_LABELS = {
    'Calculated': 'Calculated Integral Index',
    'TSUNAMI-IP': 'TSUNAMI-IP Integral Index',
    'Percent Difference': 'Percent Difference (%)'
}

# Whether or not to kill interactive legend (flask/dash applications) plots after a short amount of time. This is
# not intended for use by users, but is necessary for generating documentation properly.
generating_docs = False