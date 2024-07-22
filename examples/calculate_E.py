"""
Calculating the Similarity Index E
==================================

This example demonstrates the use of ``tsunami_ip_utils`` for calculating the similarity index :math:`E` for a set of TSUNAMI-IP 
SDF files. This example uses a set of 12 SDF files generated from the highly enriched uranium fast metal (HMF) series of 
critical experiments. :func:`tsunami_ip_utils.comparisons.E_calculation_comparison` is used to calculate the similarity index
"""
from tsunami_ip_utils.comparisons import E_calculation_comparison
from tsunami_ip_utils.viz.scatter_plot import InteractiveScatterLegend

# List of all SDF files
all_sdfs = [ f"u-235-test-case/sphere_model_{i}.sdf" for i in range(1, 9) ]
data = E_calculation_comparison("u-235-test-case/tsunami_ip.out", all_sdfs, all_sdfs)

# Save data as excel spreadsheets
for E_type in data.keys():
    data[E_type].to_excel(f"{E_type}_comparison.xlsx")

