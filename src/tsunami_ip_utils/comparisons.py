from tsunami_ip_utils.readers import read_integral_indices
from tsunami_ip_utils.integral_indices import calculate_E
import numpy as np
from uncertainties import unumpy
import pandas as pd

def comparison(tsunami_ip_output_filename, application_filenames: list, experiment_filenames: list):
    """Function that compares the calculated similarity parameter E with the TSUNAMI-IP output for each application with each
    experiment. The comparison is done for the nominal values and the uncertainties of the E values. In addition, the
    difference between manually calculated uncertainties and automatically calculated uncertainties (i.e. via the uncertainties
    package) is also calculated. The results are returned as a pandas DataFrame.
    
    Parameters
    ----------
    - tsunami_ip_output_filename: str, path to the TSUNAMI-IP output file
    - application_filenames: list of str, paths to the application sdf files
    - experiment_filenames: list of str, paths to the experiment sdf files
    
    Returns
    -------
    - data: dict, dictionary of pandas DataFrames for each type of E index. The DataFrames contain the calculated E values,
        the manual uncertainties, the TSUNAMI-IP values, the relative difference in the mean, and the relative difference
        in the manual uncertainty. The DataFrames are indexed by the experiment number and the columns are a MultiIndex
        with the application number as the main index and the attributes as the subindex."""
    
    # First perform the manual calculations for each type of E index
    E_types = ['total', 'fission', 'capture', 'scatter']
    E = {}
    for E_type in E_types + ['total_manual', 'fission_manual', 'capture_manual', 'scatter_manual']:
        if 'manual' in E_type: # Manual uncertainty propagation
            if E_type.replace('_manual', '') == 'total':
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type='all', uncertainties='manual')
            elif E_type.replace('_manual', '') == 'scatter':
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type='elastic', uncertainties='manual')
            else:
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type=E_type.replace('_manual', ''), uncertainties='manual')
        else: # Automatic uncertainty propagation
            if E_type == 'total':
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type='all')
            elif E_type == 'scatter':
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type='elastic')
            else:
                E[E_type] = calculate_E(application_filenames, experiment_filenames, reaction_type=E_type)

    print("Done with calculations")

    # Now read the tsunami_ip output
    tsunami_ip_output = read_integral_indices(tsunami_ip_output_filename)

    # Compare the nominal values
    E_diff = {}
    for E_type in E_types:
        E_diff[E_type] = np.abs(unumpy.nominal_values(E[E_type]) - unumpy.nominal_values(tsunami_ip_output[f"E_{E_type}"])) \
                            / unumpy.nominal_values(tsunami_ip_output[f"E_{E_type}"])

    # Compare the calculated (manual) uncertainty with the TSUNAMI-IP uncertainty
    E_diff_unc = {}
    for E_type in E_types:
        E_diff_unc[E_type] = np.abs( unumpy.std_devs(E[E_type + '_manual']) - unumpy.std_devs(tsunami_ip_output[f"E_{E_type}"]) ) \
                                / unumpy.std_devs(tsunami_ip_output[f"E_{E_type}"])

    # -----------------------------------------
    # Format the results as a pandas DataFrame
    # -----------------------------------------

    # Create a MultiIndex for columns
    num_experiments, num_applications = np.shape(E['total'])

    columns = pd.MultiIndex.from_product([
        np.arange(1, num_applications + 1),  # Main column indices
        ['Calculated', 'Manual Uncertainty', 'TSUNAMI-IP', 'Relative Difference in Mean', 'Relative Difference in Manual Uncertainty']  # Subcolumns
    ], names=['Application Number', 'Attribute'])

    # Initialize DataFrame
    data = {}

    print("Creating pandas dataframes")

    # Create a pandas DataFrame for each type of E index
    for E_type in E_types:
        data[E_type] = pd.DataFrame(index=pd.Index(np.arange(1, num_experiments + 1), name='Experiment Number'), columns=columns)

        # Populate DataFrame
        for application_index in range(num_applications):
            for experiment_index in range(num_experiments):
                # Now write the data to the DataFrame
                data[E_type].loc[experiment_index + 1, (application_index + 1, 'Calculated')] = \
                    f"{E[E_type][experiment_index, application_index].n:1.3E}+/-{E[E_type][experiment_index, application_index].s:1.2E}"
                data[E_type].loc[experiment_index + 1, (application_index + 1, 'Manual Uncertainty')] = \
                    f"{E[E_type + '_manual'][experiment_index, application_index].s:1.2E}"
                data[E_type].loc[experiment_index + 1, (application_index + 1, 'TSUNAMI-IP')] = \
                    f"{tsunami_ip_output[f'E_{E_type}'][experiment_index, application_index].n:1.3E}+/-{tsunami_ip_output[f'E_{E_type}'][experiment_index, application_index].s:1.2E}"
                data[E_type].loc[experiment_index + 1, (application_index + 1, 'Relative Difference in Mean')] = f"{E_diff[E_type][experiment_index, application_index]:1.4E}"
                data[E_type].loc[experiment_index + 1, (application_index + 1, 'Relative Difference in Manual Uncertainty')] = f"{E_diff_unc[E_type][experiment_index, application_index]:1.4E}"

    return data