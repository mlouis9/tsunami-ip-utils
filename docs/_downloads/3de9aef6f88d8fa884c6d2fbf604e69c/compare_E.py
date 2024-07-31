"""
Comparing Calculated and TSUNAMI-IP E for a Set of SDF Files
============================================================

This example demonstrates the use of ``tsunami_ip_utils`` for calculating the similarity index :math:`E` for a set of TSUNAMI-IP 
SDF files and comparing to the TSUNAMI-IP calculated value. This example uses a set of 12 SDF files generated from the highly 
enriched uranium fast metal (HMF) series of critical experiments [Bess2019]_.
"""
from tsunami_ip_utils.comparisons import E_calculation_comparison
import pandas as pd
import pickle
from paths import EXAMPLES

# Create a list of application and experiment SDF files (same set of applications and experiments)
applications = [ f"{EXAMPLES}/data/example_sdfs/HMF/HEU-MET-FAST-003-00{i}.sdf" for i in range(1, 3) ]
experiments = applications
data = E_calculation_comparison(applications, experiments)
E_total_comparison = data['total']

# %%
# The output is a pandas dataframe. The application is the first index, and can be accessed via passing the application
# number as an index. For example, to access the first application:

print(E_total_comparison[1])

# %%
# From the output, you can see that each application has a few attributes that are calculated. The attributes are:

print(E_total_comparison[1].keys())

# %%
# To access, for example, the calculated E value, you can pass the appropriate key

print(E_total_comparison[1]['Calculated'])

# %%
# This gives the calculated E value for application 1 compared to all other experiments. To select a specific application
# experiement pair, you can pass the experiment number as the second index. For example, to access the calculated E value
# for application 1 compared to experiment 2:

print(E_total_comparison[1]['Calculated'][2])

# %%
# Now the pandas dataframe comparison for each of the :math:`E` types (i.e. :math:`E_{\text{total}}`, :math:`E_{\text{fission}}`, 
# :math:`E_{\text{capture}}`, and :math:`E_{\text{scatter}}`) can be written to an excel file. The gold standard data is also 
# provided in the ``gold`` directory.

# Save output data as xlsx
for E_type in data.keys():
    data[E_type].to_excel(f"{E_type}_comparison.xlsx")

# Test that the calculated E comparison is the same as the gold standard
for E_type in data.keys():
    dumped_data = pd.read_excel( f"{E_type}_comparison.xlsx"  )
    gold_data   = pd.read_excel( f"gold/{E_type}_comparison.xlsx" )

    assert dumped_data.equals( gold_data )

# %%
# References
# ----------
# .. [Bess2019] J. Bess, "International handbook of evaluated criticality safety benchmark experiments (ICSBEP)," 
#   Organization for Economic Co-operation and Development-Nuclear Energy Agency Report NEA/NSC/DOC (95) 3 (2019).