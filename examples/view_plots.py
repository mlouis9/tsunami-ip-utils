"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.viz import load_interactive_matrix_plot

# This loads the plot matrix for the 1000 sample points comparison
fig = load_interactive_matrix_plot('results/HMF/perturbation_comparisons_200.pkl')
fig.show()

# Load the example matrix plot of %Δk/k contributions
# fig = load_interactive_matrix_plot('results/uncertainty_comparisons/MCT_comparisons_isotope_reaction.pkl')
# fig.show()