"""
Calculating the Similarity Index E
==================================

This example demonstrates the use of ``tsunami_ip_utils`` for calculating the similarity index :math:`E` for a set of TSUNAMI-IP 
SDF files. This example uses a set of 12 SDF files generated from the highly enriched uranium fast metal (HMF) series of 
critical experiments [Bess2019]_.
"""
from tsunami_ip_utils.comparisons import E_calculation_comparison

# List of all SDF files
all_sdfs = [ f"data/example_sdfs/sphere_model_{i}.sdf" for i in range(1, 4) ]
data = E_calculation_comparison(all_sdfs, all_sdfs)

# Save data as excel spreadsheets
for E_type in data.keys():
    data[E_type].to_excel(f"{E_type}_comparison.xlsx")

# %%
# References
# ==========
# .. [Bess2019] J. Bess, "International handbook of evaluated criticality safety benchmark experiments (ICSBEP)," 
#   Organization for Economic Co-operation and Development-Nuclear Energy Agency Report NEA/NSC/DOC (95) 3 (2019).