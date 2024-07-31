import os

def remove_header_from_rst(filepath, header_to_replace, replacement_header):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    with open(filepath, 'w') as file:
        replace_index = 0
        for index, line in enumerate(lines):
            if header_to_replace in line:
                replace_index = index
                break
        
        lines[replace_index] = replacement_header + '\n'
        file.writelines(lines)            

remove_header_from_rst('source/auto_examples/matrix_plot/index.rst', ':pseudoheader:`Matrix Plots`', 'Matrix Plots\n------------')
remove_header_from_rst('source/auto_examples/correlation_plot/index.rst', ':pseudoheader:`Correlation Plots`', 'Correlation Plots\n-----------------')