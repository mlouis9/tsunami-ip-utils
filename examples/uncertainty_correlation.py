from tsunami_ip_utils.integral_indices import calculate_uncertainty_contributions
from tsunami_ip_utils.viz import correlation_plot
from tsunami_ip_utils.readers import read_integral_indices

filenames = [ f'3d-sphere/sphere_model_{i}.out' for i in range(1, 13) ]

indices = [0, 9]
application_filenames = [ filenames[index] for index in indices ]
experiment_filenames = [ filenames[index] for index in indices ]

isotope_total, isotope_reaction = calculate_uncertainty_contributions(application_filenames, experiment_filenames)

integral_indices = read_integral_indices('3d-sphere/tsunami_ip.out')
C_k = integral_indices['C_k']
print(C_k[1, 0])
fig = correlation_plot(isotope_reaction['application'][0], isotope_reaction['experiment'][1], plot_type='interactive_scatter', \
                 integral_index_name='%Δk/k', plot_redundant_reactions=True, interactive_legend=True)
fig.save_state('test.pkl')
# fig.write_html('uncertainty_correlation.html')