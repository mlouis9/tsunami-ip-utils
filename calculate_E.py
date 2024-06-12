# Script for recreating the global similarity parameter E
from pyparsing import *
from uncertainties import ufloat
from uncertainties import unumpy, umath
import numpy as np
import math
import pandas as pd

ParserElement.enablePackrat()

def read_region_integrated_sdf(filename: str):
    """Function that reads the sdf file and returns a dictionary of nuclide-reaction pairs and energy-dependent
    sensitivities (with uncertainties)
    
    Parameters
    ----------
    - Filename: str, path to the sdf file"""
    with open(filename, 'r') as f:
        data = f.read()

    # ===========================================
    # Get region integrated sensitivity profiles
    # ===========================================

    # Number of columns in SDF profiles
    NUM_COLUMNS = 5

    # Get number of energy groups
    unused_lines = SkipTo(pyparsing_common.integer + "number of neutron groups")
    num_groups_parser = Suppress(unused_lines) + pyparsing_common.integer + Suppress("number of neutron groups")
    num_groups = num_groups_parser.parseString(data)[0]

    data_line = Group(OneOrMore(pyparsing_common.sci_real))
    data_block = OneOrMore(data_line)

    unused_lines = SkipTo("energy boundaries:")
    energy_boundaries = Suppress(unused_lines + "energy boundaries:") + data_block
    energy_boundaries = np.array(energy_boundaries.parseString(data)[0])

    # ------------------
    # SDF profile parser
    # ------------------
    atomic_number = Word(nums)
    element = Word(alphas.lower(), exact=1) 
    isotope_name = Combine(element + '-' + atomic_number)

    # Grammar for sdf header
    reaction_type = Word(alphanums + ',\'')
    zaid = Word(nums, max=6)
    reaction_mt = Word(nums, max=4)

    # Lines of the sdf header
    sdf_header_first_line = isotope_name + reaction_type + zaid + reaction_mt + Suppress(LineEnd())
    sdf_header_second_line = pyparsing_common.signed_integer + pyparsing_common.signed_integer + Suppress(LineEnd())
    sdf_header_third_line = Suppress(pyparsing_common.sci_real + pyparsing_common.sci_real + pyparsing_common.signed_integer + \
          pyparsing_common.signed_integer + LineEnd())
    
    # This line contains the total energy integrated sensitivity data for the given profile, along with uncertainties, etc.
    sdf_data_first_line = pyparsing_common.sci_real + pyparsing_common.sci_real + pyparsing_common.sci_real + \
                            pyparsing_common.sci_real + pyparsing_common.sci_real + LineEnd()
    
    # sdf_headers = sdf_header.searchString(data)
    # The total sdf header
    sdf_header = sdf_header_first_line + sdf_header_second_line + sdf_header_third_line

    # SDF profile data

    sdf_data_block = OneOrMore(data_line)
    sdf_data = sdf_header + Suppress(sdf_data_first_line) + sdf_data_block
    sdf_data = sdf_data.searchString(data)
    # Now break the data blocks into two parts: the first part is the sensitivities and the second is the uncertainties
    sdf_data = [match[:-1] + [np.array(match[-1][:num_groups]), np.array(match[-1][num_groups:])] for match in sdf_data]

    # --------------------------------------------------
    # Now only return the region integrated sdf profiles
    # --------------------------------------------------
    # i.e. those with zone number and zone volume both equal to 0

    sdf_data = [ match for match in sdf_data if match[4] == 0 and match[5] == 0 ]

    # Now parse each result into a readable dictionary
    names = ["isotope", "reaction_type", "zaid", "reaction_mt", "zone_number", "zone_volume", "sensitivities", "uncertainties"]
    sdf_data = [dict(zip(names, match)) for match in sdf_data]

    # Convert the sensitivities and uncertainties to uncertainties.ufloat objects
    for match in sdf_data:
        match["sensitivities"] = unumpy.uarray(match['sensitivities'], match['uncertainties'])
        del match["uncertainties"]

    return energy_boundaries, sdf_data

def calculate_E(application_filenames: list, experiment_filenames: list, reaction_type='all', uncertainties='automatic'):
    """Calculates the similarity parameter, E for each application with each available experiment given the application 
    and experiment sdf files
    
    Parameters
    ----------
    - application_filenames: list of str, paths to the application sdf files
    - experiment_filenames: list of str, paths to the experiment sdf files
    - reaction_type: str, the type of reaction to consider. Default is 'all' which considers all reactions
    - uncertainties: str, the type of uncertainty propagation to use. Default is 'automatic' which uses the uncertainties
        package for error propagation. If set to 'manual', then manual error propagation is used
    
    Returns
    -------
    - E: np.ndarray, similarity parameter for each application with each experiment (experiment x application)"""

    # Read the application and experiment sdf files

    if reaction_type == "all":
        application_sdfs = [ [ data['sensitivities'] for data in read_region_integrated_sdf(filename)[1] ] for filename in application_filenames ]
        experiment_sdfs = [ [ data['sensitivities'] for data in read_region_integrated_sdf(filename)[1] ] for filename in experiment_filenames ]
    else:
        application_sdfs = [ [ data['sensitivities'] for data in read_region_integrated_sdf(filename)[1] if data['reaction_type'] == reaction_type ] \
                            for filename in application_filenames ]
        experiment_sdfs = [ [ data['sensitivities'] for data in read_region_integrated_sdf(filename)[1] if data['reaction_type'] == reaction_type ] \
                            for filename in experiment_filenames ]

    def create_sensitivity_vector(sdfs):
        """Creates a senstivity vector from all of the sensitivity profiles from a specific application or experiment
        
        Parameters:
        -----------
        - sdfs: list of unumpy.uarrays of sensitivities"""
        uncertainties = np.concatenate([unumpy.std_devs(sdf) for sdf in sdfs])
        senstivities = np.concatenate([unumpy.nominal_values(sdf) for sdf in sdfs])

        return unumpy.uarray(senstivities, uncertainties)

    E_vals = unumpy.umatrix(np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ), \
                            np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ))
    
    for i, experiment in enumerate(experiment_sdfs):
        for j, application in enumerate(application_sdfs):
            application_vector = create_sensitivity_vector(application)
            experiment_vector = create_sensitivity_vector(experiment)

            if uncertainties == 'automatic':
                application_norm = umath.sqrt(np.sum(application_vector**2))
                experiment_norm = umath.sqrt(np.sum(experiment_vector**2))
                
                application_unit_vector = application_vector / application_norm
                experiment_unit_vector = experiment_vector / experiment_norm

                # For some reason, the above introduces a correlation between the application and experiment vectors
                # which should only be the case if the application and the experiment are the same, so we manually 
                # break this correlation otherwise
                
                if application_filenames[j] != experiment_filenames[i]:
                    # Break dependency to treat as independent
                    application_unit_vector = np.array([ufloat(v.n, v.s) for v in application_unit_vector])
                    experiment_unit_vector = np.array([ufloat(v.n, v.s) for v in experiment_unit_vector])

                dot_product = np.dot(application_unit_vector, experiment_unit_vector)

                E = dot_product
                E_vals[i, j] = E

            else:
                # Manual error propagation

                # Same calculation for E
                application_norm = umath.sqrt(np.sum(application_vector**2))
                experiment_norm = umath.sqrt(np.sum(experiment_vector**2))

                application_unit_vector = application_vector / application_norm
                experiment_unit_vector = experiment_vector / experiment_norm

                dot_product = np.dot(application_vector, experiment_vector)
                E = ( dot_product / \
                    (application_norm * experiment_norm)).n
                
                # ------------------------------------------
                # Now manually perform the error propagation
                # ------------------------------------------

                """Idea: propagate uncertainty of components of each sensitivity vector to the normalized sensitivty vector
                i.e. u / ||u|| = u^, v / ||v|| = v^
                
                E = u^ . v^

                Since u^ and v^ are completely uncorrelated, we may first calculate the uncertainty of u^ and v^ separately then
                use those to calculate the uncertainty of E. To caluculate the uncertainty of u^ and v^, note

                σ_{u^} = √( ∑_j ( ∂u_i / ∂u_j * σ_{u_j} )^2 )

                Since u^ = u_i / √( u_1^2 + ... + u_n^2 ), we have

                For i=j:

                    ∂u_i / ∂u_j = ( ||u||^2 - u_i^2 ) / ||u||^3

                For i!=j:

                    ∂u_i / ∂u_j = -u_i * u_j / ||u||^3
                """

                def unit_vector_uncertainty_propagation(vector):
                    """Does error propagation for the components of a vector v that is normalized to a unit vector via 
                    v^ = v / ||v||"""
                    # Calculate norm of the vector
                    vector_norm = np.sqrt(np.sum(unumpy.nominal_values(vector)**2))

                    # Extract the uncertainties of the vector components
                    vector_uncertainties = unumpy.std_devs(vector)

                    # Compute the derivative matrix for uncertainty propagation
                    derivative_matrix = -unumpy.nominal_values(vector)[:, np.newaxis] * unumpy.nominal_values(vector)[np.newaxis, :] / vector_norm**3
                    np.fill_diagonal(derivative_matrix, (vector_norm**2 - unumpy.nominal_values(vector)**2) / vector_norm**3)

                    # Calculate the uncertainties of the unit vector components
                    unit_vector_uncertainties = np.sqrt(np.sum((derivative_matrix**2 * vector_uncertainties**2), axis=1))

                    # Return the unit vector with uncertainties
                    return unit_vector_uncertainties

                def dot_product_uncertainty_propagation(vect1, vect2):
                    """Calculates the uncertainty in the dot product of two vectors with uncertainties"""
                    dot_product_uncertainty = 0
                    for i in range(len(vect1)):
                        if vect1[i].n == 0 or vect2[i].n == 0:
                            product_uncertainty = 0
                        else:
                            product_uncertainty = vect1[i].n * vect2[i].n * np.sqrt((vect1[i].s/vect1[i].n)**2 + \
                                                                                    (vect2[i].s/vect2[i].n)**2)
                        dot_product_uncertainty += product_uncertainty**2
                    dot_product_uncertainty = np.sqrt(dot_product_uncertainty)

                    return dot_product_uncertainty

                application_unit_vector_error = unit_vector_uncertainty_propagation(application_vector)
                experiment_unit_vector_error = unit_vector_uncertainty_propagation(experiment_vector)

                application_unit_vector = unumpy.uarray(unumpy.nominal_values(application_unit_vector), \
                                                         application_unit_vector_error)
                experiment_unit_vector = unumpy.uarray(unumpy.nominal_values(experiment_unit_vector), experiment_unit_vector_error)

                # Now calculate error in dot product to get the uncertainty in E
                E_uncertainty = dot_product_uncertainty_propagation(application_unit_vector, experiment_unit_vector)

                E_vals[i, j] = ufloat(E, E_uncertainty)

    return E_vals

def read_tsunami_ip_ouput(filename):
    with open(filename, 'r') as f:
        data = f.read()

    # Define the Integral Values parser
    dashed_line = OneOrMore("-")
    header = Literal("Integral Values for Application") + "#" + pyparsing_common.integer + LineEnd() + dashed_line
    table_header = Literal("Experiment") + Literal("Type") + Literal("Value") + Literal("s.d.") + Literal("c(k)") + \
                    Literal("s.d.") + Literal("E") + Literal("s.d.") + Literal("E(fis)") + Literal("s.d.") + Literal("E(cap)") + \
                    Literal("s.d.") + Literal("E(sct)") + Literal("s.d.") + LineEnd() + OneOrMore(dashed_line)
    
    # Define characters allowed in a filename (all printables except space)
    non_space_printables = ''.join(c for c in printables if c != ' ')
    sci_num = pyparsing_common.sci_real
    data_line = Group(Suppress(pyparsing_common.integer + Word(non_space_printables) + Word(alphas)) + \
                        Group(sci_num + sci_num) + Group(sci_num + sci_num) + Group(sci_num + sci_num) + \
                        Group(sci_num + sci_num) + Group(sci_num + sci_num) + Group(sci_num + sci_num))
    data_block = OneOrMore(data_line)
    integral_values = Suppress(header + table_header) + data_block
    parsed_integral_values = integral_values.searchString(data)

    # Parse the integral value tables into a uarray
    num_applications = len(parsed_integral_values)
    num_experiments = len(parsed_integral_values[0]) - 1 # First row seems to be a repeat
    integral_matrices = {}
    
    integral_matrix = unumpy.umatrix(np.zeros((num_experiments, num_applications)), np.zeros((num_experiments, num_applications)))
    C_k       = np.copy(integral_matrix)
    E_total   = np.copy(integral_matrix)
    E_fission = np.copy(integral_matrix)
    E_capture = np.copy(integral_matrix)
    E_scatter = np.copy(integral_matrix)

    for match_index, match in enumerate(parsed_integral_values):
        for row_index, row in enumerate(match[1:]):
            C_k[row_index, match_index]       = ufloat(row[1][0], row[1][1])
            E_total[row_index, match_index]   = ufloat(row[2][0], row[2][1])
            E_fission[row_index, match_index] = ufloat(row[3][0], row[3][1])
            E_capture[row_index, match_index] = ufloat(row[4][0], row[4][1])
            E_scatter[row_index, match_index] = ufloat(row[5][0], row[5][1])

    integral_matrices.update({
        "C_k": C_k,
        "E_total": E_total,
        "E_fission": E_fission,
        "E_capture": E_capture,
        "E_scatter": E_scatter
    })

    return integral_matrices

def comparison(tsunami_ip_output_filename, application_filenames: list, experiment_filenames: list):
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
    tsunami_ip_output = read_tsunami_ip_ouput(tsunami_ip_output_filename)

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

if __name__ == '__main__':
    all_sdfs = [ f"sphere_model_{i}.sdf" for i in range(1, 9) ]
    data = comparison("tsunami_ip.out", all_sdfs, all_sdfs)
    for E_type in data.keys():
        data[E_type].to_excel(f"{E_type}_comparison.xlsx")

