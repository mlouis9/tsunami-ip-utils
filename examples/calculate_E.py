"""
Placeholder
===========

This is a placeholder description
"""
# Script for recreating the global similarity parameter E
from tsunami_ip_utils.comparisons import E_calculation_comparison

if __name__ == '__main__':
    # List of all SDF files
    all_sdfs = [ f"u-235-test-case/sphere_model_{i}.sdf" for i in range(1, 9) ]
    data = E_calculation_comparison("u-235-test-case/tsunami_ip.out", all_sdfs, all_sdfs)

    # Save data as excel spreadsheets
    for E_type in data.keys():
        data[E_type].to_excel(f"{E_type}_comparison.xlsx")

