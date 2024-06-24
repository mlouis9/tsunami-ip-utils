from tsunami_ip_utils.viz import correlation_plot, contribution_plot
from tsunami_ip_utils.integral_indices import calculate_E_contributions, add_missing_reactions_and_nuclides
import matplotlib.pyplot as plt
from tsunami_ip_utils.readers import RegionIntegratedSdfReader

all_sdfs = [ f"3d-sphere/sphere_model_{i}.sdf" for i in range(1, 13) ]

indices = [0, 9]
application_filenames = [ all_sdfs[index] for index in indices ]
experiment_filenames = [ all_sdfs[index] for index in indices ]

# E = calculate_E(application_filenames, experiment_filenames)
# print(E)
E_contributions_nuclide, E_contributions_nuclide_reaction = calculate_E_contributions(application_filenames, experiment_filenames)

# fig = contribution_plot(E_contributions_nuclide_reaction[0, 0], integral_index_name='E', plot_type='interactive_pie', plot_redundant_reactions=True)
# fig.show()

# fig = contribution_plot(E_contributions_nuclide_reaction[1, 1], integral_index_name='E', plot_type='interactive_pie', plot_redundant_reactions=True)
# fig.show()
print(E_contributions_nuclide_reaction[0, 0])
fig = correlation_plot(E_contributions_nuclide_reaction[0, 0], E_contributions_nuclide_reaction[1, 1], plot_type='interactive_scatter', \
                 integral_index_name='E', plot_redundant_reactions=True, interactive_legend=True)
fig.write_html('correlation_plot.html')