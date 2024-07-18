"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.viz import generate_heatmap_from_comparison

suffix = ['calculated', 'tsunami_ip', 'percent_diff']
keys = ['Calculated', 'Tsunami-IP', 'Percent Difference']
experiment_types = ['HMF', 'MCT']
fontsizes = [6, 4]
contribution_types = ['isotope', 'isotope_reaction']

# Generate heatmaps for uncertainty contributions
for contribution_type in contribution_types:
    for experiment_type, fontsize in zip(experiment_types, fontsizes):
        plot_dict = generate_heatmap_from_comparison(
            f"./results/uncertainty_comparisons/{experiment_type}_comparisons_{contribution_type}.xlsx",
            base_fontsize=fontsize)
        for index, value in enumerate(plot_dict.values()):
            fig, ax = value
            if fontsize == 4:
                ax.tick_params(axis='x', labelsize=6)
                ax.tick_params(axis='y', labelsize=6)
            fig.savefig(f'./results/heatmaps/{contribution_type}_wise_{experiment_type}_{suffix[index]}.png', dpi=300, bbox_inches='tight')

# ------------------------------------
# Generate heatmaps for perturbations
# ------------------------------------

# MCT cases
plot_dict = generate_heatmap_from_comparison(
    f"./results/MCT/perturbation_comparisons_500.xlsx",
    base_fontsize=4)
for index, value in enumerate(plot_dict.values()):
    fig, ax = value
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    fig.savefig(f'./results/heatmaps/perturbation_MCT_{suffix[index]}.png', dpi=300, bbox_inches='tight')

# HMF cases
plot_dict = generate_heatmap_from_comparison(
    f"./results/HMF/perturbation_comparisons_500.xlsx",
    base_fontsize=4)
for index, value in enumerate(plot_dict.values()):
    fig, ax = value
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    fig.savefig(f'./results/heatmaps/perturbation_HMF_{suffix[index]}.png', dpi=300, bbox_inches='tight')
# points_array = [50, 200, 500, 1000]
# for points in points_array:
#     plot_dict = generate_heatmap_from_comparison(
#         f"./results/HMF/perturbation_comparisons_{points}.xlsx",
#         base_fontsize=4)
#     ax.tick_params(axis='x', labelsize=6)
#     ax.tick_params(axis='y', labelsize=6)
#     fig.savefig(f'./results/heatmaps/perturbation_HMF_{points}_percent_diff.png', dpi=300, bbox_inches='tight')