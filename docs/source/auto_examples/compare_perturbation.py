"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.comparisons import correlation_comparison
from tsunami_ip_utils.readers import read_integral_indices
from pathlib import Path
import numpy as np

# Define the filenames for the applications and experiments
filenames = [
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-001-001.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-001-002.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-001-003.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-001-004.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-001S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-002S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-003S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-004S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-005S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-002-006S.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-001.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-002.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-003.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-004.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-005.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-006.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-007.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-008.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-009.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-010.sdf"),
    Path("/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/MIX-COMP-THERM-004-011.sdf"),
]

indices = range(21)
experiment_filenames = [filenames[i] for i in indices]
application_filenames = [filenames[i] for i in indices]

# Get the TSUNAMI-IP integral indices for the applications and experiments
integral_indices = read_integral_indices('/home/mlouis9/projects/suli2024/main/Matthew_VALID_sdfs/MCT/tsunami_ip.out')
C_k = integral_indices['C_k']
C_k = C_k[np.ix_(indices, indices)] # Slice the C_k matrix to only include the indices of interest

# Define paths used for generating cross section perturbations
multigroup_library = Path("~/codes/SCALE-6.3.1/data/scale.rev01.xn252v8.0")
perturbation_factors = Path("/home/mlouis9/codes/SCALE-6.3.1/data/perturb/252n.v7.1")

# Perform the comparison
for num_perturbations in [200, 500, 1000]:
    comparisons = correlation_comparison(
        integral_index_matrix=C_k,
        integral_index_name='c_k',
        application_files=application_filenames, 
        experiment_files=experiment_filenames, 
        method='perturbation',
        base_library=multigroup_library,
        perturbation_factors=perturbation_factors,
        num_perturbations=num_perturbations,
        make_plot=False,
        num_cores=2
    )

    # fig.save_state(f'results/allHMFcases/perturbation_comparison_{num_perturbations}.pkl')
    comparisons.to_excel(f'results/MCT/perturbation_comparison_{num_perturbations}.xlsx')
    # fig.save_state(f'results/MCT/perturbation_comparison_{num_perturbations}.pkl')
