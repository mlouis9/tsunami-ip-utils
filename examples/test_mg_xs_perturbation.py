from tsunami_ip_utils.xs import perturb_multigroup_xs_dump
from pathlib import Path
script_dir = Path(__file__).parent.absolute()

input_filename = script_dir / 'data' / 'u235-252g-all-rxns.txt'

perturbed_xs = perturb_multigroup_xs_dump(input_filename, 0.1, output_file='./data/perturbed_library.txt')