"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.viz import generate_heatmap_from_comparison

# Load the Excel file into a DataFrame with MultiIndex columns
file_path = 'results/uncertainty_comparisons/HMF_comparison.xlsx'

plot_dict = generate_heatmap_from_comparison(file_path)

fig, axs = plot_dict['Percent Difference']
axs.set_title('test')
fig.show()