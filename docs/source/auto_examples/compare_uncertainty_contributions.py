"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.comparisons import correlation_comparison
from tsunami_ip_utils.readers import read_integral_indices
import numpy as np
from pathlib import Path

filenames = [ Path(f'../3d-sphere/sphere_model_{i}.sdf') for i in range(1,13) ]

indices = range(12)
experiment_filenames = [filenames[i] for i in indices]
application_filenames = [filenames[i] for i in indices]

# Get the TSUNAMI-IP integral indices for the applications and experiments
integral_indices = read_integral_indices('../3d-sphere/tsunami_ip.out')
C_k = integral_indices['C_k']
C_k = C_k[np.ix_(indices, indices)] # Slice the C_k matrix to only include the indices of interest

comparisons, fig = correlation_comparison(
    integral_index_matrix=C_k,
    integral_index_name='c_k',
    application_files=application_filenames, 
    experiment_files=experiment_filenames, 
    method='uncertainty_contributions_nuclide'
)

fig.save_state('results/uncertainty_comparisons/HMF_comparisons_isotope.pkl')
comparisons.to_excel('results/uncertainty_comparisons/HMF_comparisons_isotope.xlsx')
fig.show()