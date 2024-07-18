"""
Placeholder
===========

This is a placeholder description
"""
from tsunami_ip_utils.xs import read_multigroup_xs
import numpy as np
from pathlib import Path

# xs, energies = read_multigroup_xs('data/u235-mg-252g.txt', energy_boundaries=True)
# print(np.shape(energies))

current_dir = Path(__file__).parent
# xs = read_nuclide_reaction_from_multigroup_library(Path('~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1'), '6000', 102)
# print(xs)

# out = parse_reactions_from_nuclide(str(current_dir / 'data' / 'u235-252g-all-rxns.txt'), reaction_mts = ['102', '18'], energy_boundaries=False)

# out = read_reactions_from_nuclide(Path('~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1'), '92235', ['102', '18'])
nuclide_reaction_dict = {'92235': ['102', '18'], '92238': ['102', '18'], '92234': ['102', '18'], '92232': ['102', '18'], '94239': ['102', '18']}
out = read_multigroup_xs(Path('~/codes/SCALE-6.3.1/data/scale.rev05.xn252v7.1'), nuclide_reaction_dict, method='small')
print(out)