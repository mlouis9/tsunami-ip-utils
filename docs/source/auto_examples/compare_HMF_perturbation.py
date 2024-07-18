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
filenames = [ Path(f'../3d-sphere/sphere_model_{i}.sdf') for i in range(1,13) ]

indices = range(12)
experiment_filenames = [filenames[i] for i in indices]
application_filenames = [filenames[i] for i in indices]

# Get the TSUNAMI-IP integral indices for the applications and experiments
integral_indices = read_integral_indices('../3d-sphere/tsunami_ip.out')
C_k = integral_indices['C_k']
C_k = C_k[np.ix_(indices, indices)] # Slice the C_k matrix to only include the indices of interest

# Define paths used for generating cross section perturbations
multigroup_library = Path("~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1")
perturbation_factors = Path("/home/mlouis9/codes/SCALE-6.3.1/data/perturb/252n.v7.1")

# Perform the comparison
for num_perturbations in [50, 200, 500, 1000]:
    comparisons = correlation_comparison(
        integral_index_matrix=C_k,
        integral_index_name='c_k',
        application_files=application_filenames, 
        experiment_files=experiment_filenames, 
        method='perturbation',
        base_library=multigroup_library,
        perturbation_factors=perturbation_factors,
        num_perturbations=num_perturbations,
        make_plot=False
    )

    # fig.save_state(f'results/allHMFcases/perturbation_comparison_{num_perturbations}.pkl')
    comparisons.to_excel(f'results/allHMFcases/perturbation_comparison_{num_perturbations}.xlsx')
