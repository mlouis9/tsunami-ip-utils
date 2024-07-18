"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.perturbations import generate_points, cache_all_libraries
from pathlib import Path
from tsunami_ip_utils.viz import perturbation_plot
from tsunami_ip_utils.viz.scatter_plot import EnhancedPlotlyFigure
from tsunami_ip_utils.readers import read_integral_indices, RegionIntegratedSdfReader
import pickle
import matplotlib.pyplot as plt
from uncertainties import unumpy
import h5py
from tsunami_ip_utils.comparisons import _update_annotation

multigroup_library = Path("~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1")
perturbation_factors = Path("/home/mlouis9/codes/SCALE-6.3.1/data/perturb/252n.v7.1")

# cache_all_libraries(multigroup_library, perturbation_factors)
# points = generate_points(Path("../3d-sphere/sphere_model_1.sdf"), Path("../3d-sphere/sphere_model_8.sdf"), multigroup_library, \
#             perturbation_factors, 50)

# with open("points.pkl", "wb") as f:
#     pickle.dump(points, f)

with open("points.pkl", "rb") as f:
    points = pickle.load(f)

integral_indices = read_integral_indices('../3d-sphere/tsunami_ip.out')
C_k = integral_indices['C_k']
E = integral_indices['E_total']

experiment_index = 7
print(C_k[experiment_index, 0], E[experiment_index, 0])

fig = perturbation_plot(points)

_update_annotation(fig, C_k[experiment_index, 0], 'c_k')
fig.show()