"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.integral_indices import get_uncertainty_contributions
from pathlib import Path

sdfs = [ Path(f'../3d-sphere/sphere_model_{i}.sdf') for i in range(1, 13) ]

dk_over_k_nuclide_wise, dk_over_k_nuclide_reaction_wise = get_uncertainty_contributions(sdfs, sdfs)
print(dk_over_k_nuclide_wise['application'])