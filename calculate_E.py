# Script for recreating the global similarity parameter E
from tsunami_ip_utils import comparison

if __name__ == '__main__':
    all_sdfs = [ f"sphere_model_{i}.sdf" for i in range(1, 9) ]
    data = comparison("tsunami_ip.out", all_sdfs, all_sdfs)
    for E_type in data.keys():
        data[E_type].to_excel(f"{E_type}_comparison.xlsx")

