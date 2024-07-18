"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.readers import read_uncertainty_contributions_out
from tsunami_ip_utils.readers import read_uncertainty_contributions_sdf
from pathlib import Path

isotope_totals, isotope_reaction = read_uncertainty_contributions_out('../3d-sphere/sphere_model_9.out')
print(isotope_reaction)

isotope_totals, isotope_reaction = read_uncertainty_contributions_sdf([Path('../3d-sphere/sphere_model_9.sdf')])
print(isotope_totals[0][0])