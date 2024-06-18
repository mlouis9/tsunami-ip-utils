from pyparsing import *
from uncertainties import ufloat, unumpy, umath
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
from abc import ABC, abstractmethod
from threading import Timer
import webbrowser
import os, sys, signal
from flask import Flask, render_template_string
import threading

ParserElement.enablePackrat()
plt.rcParams['hatch.linewidth'] = 0.5

# ==================
# Reading Utilities
# ==================


class SdfReader:
    def __init__(self, filename):
        self.energy_boundaries, self.sdf_data = self._read_sdf(filename)
        
    def _read_sdf(self, filename):
        """Function that reads the sdf file and returns a dictionary of nuclide-reaction pairs and energy-dependent
        sensitivities (with uncertainties)
        
        Parameters
        ----------
        - Filename: str, path to the sdf file
        
        Returns
        -------
        - energy_boundaries: np.ndarray, energy boundaries for the energy groups
        - sdf_data: list of dict, list of dictionaries containing the nuclide-reaction pairs and the sensitivities
            and uncertainties for the region-integrated sensitivity profile for each nuclide-reaction pair. The dictionary
            keys are 'isotope', 'reaction_type', 'zaid', 'reaction_mt', 'zone_number', 'zone_volume', 
            'energy_integrated_sensitivity', 'abs_sum_groupwise_sensitivities', 'sum_opposite_sign_groupwise_sensitivities', 
            'sensitivities', and 'uncertainties'. The sensitivities and uncertainties are stored as unumpy.uarray objects."""
        with open(filename, 'r') as f:
            data = f.read()

        # ========================
        # Get sensitivity profiles
        # ========================

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
        sdf_header_third_line = pyparsing_common.sci_real + pyparsing_common.sci_real + pyparsing_common.signed_integer + \
                                    pyparsing_common.signed_integer + LineEnd()
        
        # This line contains the total energy integrated sensitivity data for the given profile, along with uncertainties, etc.
        sdf_data_first_line = Group(pyparsing_common.sci_real + pyparsing_common.sci_real) + \
                              pyparsing_common.sci_real + Group(pyparsing_common.sci_real + pyparsing_common.sci_real) + \
                              Suppress(LineEnd())
        
        # The total sdf header
        sdf_header = sdf_header_first_line + sdf_header_second_line + Suppress(sdf_header_third_line)

        # SDF profile data

        sdf_data_block = OneOrMore(data_line)
        sdf_data = sdf_header + sdf_data_first_line + sdf_data_block
        sdf_data = sdf_data.searchString(data)

        # Now break the data blocks into two parts: the first part is the sensitivities and the second is the uncertainties
        sdf_data = [match[:-1] + [np.array(match[-1][:num_groups]), np.array(match[-1][num_groups:])] for match in sdf_data]

        # -------------------------------------------------
        # Now parse each result into a readable dictionary
        # -------------------------------------------------

        # NOTE sum_opposite_sign_groupwise_sensitivities referrs to the groupwise sensitivities with opposite sign to the
        # integrated sensitivity coefficient
        names = ["isotope", "reaction_type", "zaid", "reaction_mt", "zone_number", "zone_volume", \
                 "energy_integrated_sensitivity", "abs_sum_groupwise_sensitivities", \
                 "sum_opposite_sign_groupwise_sensitivities", "sensitivities", "uncertainties"]
        sdf_data = [dict(zip(names, match)) for match in sdf_data]

        # Convert the sensitivities and uncertainties to uncertainties.ufloat objects
        for match in sdf_data:
            match["sensitivities"] = unumpy.uarray(match['sensitivities'], match['uncertainties'])
            match["energy_integrated_sensitivity"] = \
                ufloat(match['energy_integrated_sensitivity'][0], match['energy_integrated_sensitivity'][1])
            match["sum_opposite_sign_groupwise_sensitivities"] = \
                ufloat(match['sum_opposite_sign_groupwise_sensitivities'][0], match['sum_opposite_sign_groupwise_sensitivities'][1])
            del match["uncertainties"]

        return energy_boundaries, sdf_data
        
class RegionIntegratedSdfReader(SdfReader):
    def __init__(self, filename):
        super().__init__(filename)
        
        # Now only return the region integrated sdf profiles
        # i.e. those with zone number and zone volume both equal to 0
        self.filename = filename
        self.sdf_data = [ match for match in self.sdf_data if match['zone_number'] == 0 and match['zone_volume'] == 0 ]
    
    def convert_to_dict(self):
        # Transform the data into a dictionary keyed by nuclide and reaction type. Since data is region and mixture integrated
        # we can assume that there is only one entry for each nuclide-reaction pair
        sdf_data_dict = {}
        for match in self.sdf_data:
            nuclide = match['isotope']
            reaction_type = match['reaction_type']
            
            if nuclide not in sdf_data_dict:
                sdf_data_dict[nuclide] = {}

            sdf_data_dict[nuclide][reaction_type] = match
            
        self.sdf_data = sdf_data_dict
        return self
    
    def get_sensitivity_profiles(self, reaction_type='all'):
        """Returns the sensitivity profiles for each nuclide-reaction pair in a list
        
        Parameters
        ----------
        - reaction_type: str, the type of reaction to consider. Default is 'all' which considers all reactions
        
        Returns
        -------
        - sensitivity_profiles: list of unumpy.uarrays, list of sensitivity profiles for each nuclide-reaction pair"""
        if reaction_type == 'all':
            return [ data['sensitivities'] for data in RegionIntegratedSdfReader(self.filename).sdf_data ]
        else:
            return [ data['sensitivities'] for data in RegionIntegratedSdfReader(self.filename).sdf_data \
                    if data['reaction_type'] == reaction_type ]

def read_covariance_matrix(filename: str):
    pass

def read_ck_contributions(filename: str):
    pass

def read_integral_indices(filename):
    """Reads the output file from TSUNAMI-IP and returns the integral values for each application
    
    Parameters
    ----------
    - filename: str, path to the TSUNAMI-IP output file
    
    Returns
    -------
    - integral_matrices: dict, dictionary of integral matrices for each integral index type. The dimensions
        of the matrices are (num_experiments x num_applications) where num_experiments is the number of experiments"""

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
    num_experiments = len(parsed_integral_values[0]) - 1 # First row seems to be a repeat, i.e. in the output it's "experiment 0"
    
    integral_matrices = {}
    integral_matrix = unumpy.umatrix( np.zeros( (num_experiments, num_applications) ), 
                                      np.zeros( (num_experiments, num_applications) ) )
    
    # Initialize the integral matrices
    C_k       = np.copy(integral_matrix)
    E_total   = np.copy(integral_matrix)
    E_fission = np.copy(integral_matrix)
    E_capture = np.copy(integral_matrix)
    E_scatter = np.copy(integral_matrix)

    # Now populate the integral matrices from the parsed output
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


# ============================
# Error Propagation Functions
# ============================


def unit_vector_uncertainty_propagation(vector):
    """Does error propagation for the components of a vector v that is normalized to a unit vector via 
    v^ = v / ||v||
    
    Parameters
    ----------
    - vector: unumpy.uarray, vector to normalize and calculate uncertainties for
    
    Returns
    -------
    - unit_vector_uncertainties: np.ndarray, uncertainties of the unit vector components"""


    """To calculate the uncertainty of u^ and v^, note

    σ_{u^} = √( ∑_j ( ∂u_i / ∂u_j * σ_{u_j} )^2 )

    Since u^ = u_i / √( u_1^2 + ... + u_n^2 ), we have

    For i=j:

        ∂u_i / ∂u_j = ( ||u||^2 - u_i^2 ) / ||u||^3

    For i!=j:

        ∂u_i / ∂u_j = -u_i * u_j / ||u||^3"""

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
    """Calculates the uncertainty in the dot product of two vectors with uncertainties
    
    Parameters
    ----------
    - vect1: unumpy.uarray, first vector in dot product
    - vect2: unumpy.uarray, second vector in dot product

    Returns
    -------
    - dot_product_uncertainty: float, uncertainty in the dot product of the two vectors"""

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


# ============================
# Format Conversion Utilities
# ============================


def isotope_reaction_list_to_nested_dict(isotope_reaction_list, field_of_interest):
    """Converts a list of dictionaries containing isotope-reaction pairs (and some other key that represents a value of
    interest, e.g. an sdf profile or a contribution) to a nested dictionary
    
    Parameters
    ----------
    - isotope_reaction_list: list of dict, list of dictionaries containing isotope-reaction pairs and some other key
    - field_of_interest: str, the key in the dictionary that represents the value of interest

    Returns
    -------
    - nested_dict: dict, nested dictionary containing the isotope-reaction pairs and the value of interest"""

    isotope_reaction_dict = {}

    def get_atomic_number(isotope):
        return int(re.findall(r'\d+', isotope)[0])
    
    # Sort isotopes by atomic number so plots will have similar colors across different calls
    all_isotopes = list(set([isotope_reaction['isotope'] for isotope_reaction in isotope_reaction_list]))
    all_isotopes.sort(key=get_atomic_number)
    isotope_reaction_dict = { isotope: {} for isotope in all_isotopes }

    for isotope_reaction in isotope_reaction_list:
        isotope = isotope_reaction['isotope']
        reaction = isotope_reaction['reaction_type']
        value = isotope_reaction[field_of_interest]

        isotope_reaction_dict[isotope][reaction] = value

    return isotope_reaction_dict


# =====================================
# Integral Index Calculation Utilities
# =====================================

def calculate_E_from_sensitivity_vecs(application_vector, experiment_vector, application_filename=None, \
                                      experiment_filename=None, uncertainties='automatic', experiment_norm=None, \
                                        application_norm=None):
    """Calculates E given the sensitivity vectors for an application and an experiment. 
    
    NOTE the application and experiment
    filenames are used to break the correlation between the application and experiment vectors (that should not exist). This 
    is because the vectors are only correlated if the application and experiment are the same.
    
    Parameters
    ----------
    - application_vector: unumpy.uarray, sensitivity vector for the application
    - experiment_vector: unumpy.uarray, sensitivity vector for the experiment
    - application_filename: str, filename of the application sdf file (only needed for automatic uncertianty propagation)
    - experiment_filename: str, filename of the experiment sdf file (only needed for automatic uncertianty propagation)
    - uncertainties: str, type of error propagation to use. Default is 'automatic' which uses the uncertainties package.
        If set to 'manual', then manual error propagation is used which is generally faster
    - experiment_norm: ufloat, norm of the experiment vector. If not provided, it is calculated. This is mainly used for
        calculating E contributions, where the denominator is not actually the norm of the application and experiment
        vectors.
    - application_norm: ufloat, norm of the application vector. If not provided, it is calculated. This is mainly used for
        calculating E contributions, where the denominator is not actually the norm of the application and experiment
        vectors.
        
    Returns
    -------
    - E: ufloat, similarity parameter between the application and the experiment"""
    
    norms_not_provided = ( experiment_norm == None ) and ( application_norm == None )
    if uncertainties == 'automatic': # Automatic error propagation with uncertainties package
        if application_filename == None or experiment_filename == None:
            raise ValueError("Application and experiment filenames must be provided for automatic error propagation")
        
        if norms_not_provided:
            application_norm = umath.sqrt(np.sum(application_vector**2))
            experiment_norm = umath.sqrt(np.sum(experiment_vector**2))
        
        application_unit_vector = application_vector / application_norm
        experiment_unit_vector = experiment_vector / experiment_norm

        # For some reason, the above introduces a correlation between the application and experiment vectors
        # which should only be the case if the application and the experiment are the same, so we manually 
        # break this correlation otherwise
        
        if application_filename != experiment_filename:
            # Break dependency to treat as independent
            application_unit_vector = np.array([ufloat(v.n, v.s) for v in application_unit_vector])
            experiment_unit_vector = np.array([ufloat(v.n, v.s) for v in experiment_unit_vector])

        dot_product = np.dot(application_unit_vector, experiment_unit_vector)

        E = dot_product
        return E

    else:
        # Manual error propagation

        # Same calculation for E
        if norms_not_provided:
            application_norm = umath.sqrt(np.sum(application_vector**2))
            experiment_norm = umath.sqrt(np.sum(experiment_vector**2))

        application_unit_vector = application_vector / application_norm
        experiment_unit_vector = experiment_vector / experiment_norm

        dot_product = np.dot(application_unit_vector, experiment_unit_vector)
        E = dot_product.n
        
        # ------------------------------------------
        # Now manually perform the error propagation
        # ------------------------------------------

        """Idea: propagate uncertainty of components of each sensitivity vector to the normalized sensitivty vector
        i.e. u / ||u|| = u^, v / ||v|| = v^
        
        E = u^ . v^

        Since u^ and v^ are completely uncorrelated, we may first calculate the uncertainty of u^ and v^ separately then
        use those to calculate the uncertainty of E."""

        # Calculate the uncertainties in the unit vectors
        application_unit_vector_error = unit_vector_uncertainty_propagation(application_vector)
        experiment_unit_vector_error = unit_vector_uncertainty_propagation(experiment_vector)

        # Construct application and experiment unit vectors as unumpy.uarray objects
        application_unit_vector = unumpy.uarray(unumpy.nominal_values(application_unit_vector), \
                                                    application_unit_vector_error)
        experiment_unit_vector = unumpy.uarray(unumpy.nominal_values(experiment_unit_vector), experiment_unit_vector_error)

        # Now calculate error in dot product to get the uncertainty in E
        E_uncertainty = dot_product_uncertainty_propagation(application_unit_vector, experiment_unit_vector)

        return ufloat(E, E_uncertainty)


def create_sensitivity_vector(sdfs):
    """Creates a senstivity vector from all of the sensitivity profiles from a specific application or experiment
    
    Parameters:
    -----------
    - sdfs: list of unumpy.uarrays of sensitivities for each nuclide-reaction pair in consideration
    
    Returns
    -------
    - sensitivity_vector: unumpy.uarray of sensitivities from all of the sensitivity profiles"""
    uncertainties = np.concatenate([unumpy.std_devs(sdf) for sdf in sdfs])
    senstivities = np.concatenate([unumpy.nominal_values(sdf) for sdf in sdfs])

    return unumpy.uarray(senstivities, uncertainties)


def calculate_E(application_filenames: list, experiment_filenames: list, reaction_type='all', uncertainties='manual'):
    """Calculates the similarity parameter, E for each application with each available experiment given the application 
    and experiment sdf files
    
    Parameters
    ----------
    - application_filenames: list of str, paths to the application sdf files
    - experiment_filenames: list of str, paths to the experiment sdf files
    - reaction_type: str, the type of reaction to consider in teh calculation of E. Default is 'all' which considers all 
        reactions
    - uncertainties: str, the type of uncertainty propagation to use. Default is 'automatic' which uses the uncertainties
        package for error propagation. If set to 'manual', then manual error propagation is used
    
    Returns
    -------
    - E: np.ndarray, similarity parameter for each application with each experiment (experiment x application)"""

    # Read the application and experiment sdf files

    application_sdfs = [ RegionIntegratedSdfReader(filename).get_sensitivity_profiles(reaction_type) \
                        for filename in application_filenames ]
    experiment_sdfs  = [ RegionIntegratedSdfReader(filename).get_sensitivity_profiles(reaction_type) \
                        for filename in experiment_filenames ]

    # Create a matrix to store the similarity parameter E for each application with each experiment
    E_vals = unumpy.umatrix(np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ), \
                            np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ))
    
    # Now calculate the similarity parameter E for each application with each experiment
    for i, experiment in enumerate(experiment_sdfs):
        for j, application in enumerate(application_sdfs):
            application_vector = create_sensitivity_vector(application)
            experiment_vector = create_sensitivity_vector(experiment)

            E_vals[i, j] = calculate_E_from_sensitivity_vecs(application_vector, experiment_vector, \
                                                             application_filenames[j], experiment_filenames[i], uncertainties)

    return E_vals

def get_reaction_wise_E_contributions(application, experiment, isotope, all_reactions, application_norm, experiment_norm):
    """Calculate contributions to the similarity parameter E for each reaction type for a given isotope
    
    Parameters
    ----------
    - application: dict, dictionary of application sensitivity profiles
    - experiment: dict, dictionary of experiment sensitivity profiles
    - isotope: str, isotope to consider
    - all_reactions: list of str, list of all possible reaction types
    - application_norm: ufloat, norm of the application sensitivity vector
    - experiment_norm: ufloat, norm of the experiment sensitivity vector
    
    Returns
    -------
    - E_contributions: list of dict, list of dictionaries containing the contribution to the similarity parameter E for each
        reaction type"""
    
    E_contributions = []
    for reaction in all_reactions:
        application_vector = application[isotope][reaction]['sensitivities']
        experiment_vector = experiment[isotope][reaction]['sensitivities']

        E_contribution = calculate_E_from_sensitivity_vecs(application_vector, experiment_vector, uncertainties='manual', \
                                                           application_norm=application_norm, experiment_norm=experiment_norm)
        E_contributions.append({
            "isotope": isotope,
            "reaction_type": reaction,
            "contribution": E_contribution
        })

    return E_contributions


def get_nuclide_and_reaction_wise_E_contributions(application: RegionIntegratedSdfReader, experiment: RegionIntegratedSdfReader):
    """Calculate the contributions to the similarity parameter E for each nuclide and for each reaction type for a given
    application and experiment
    
    Parameters
    ----------
    - application: RegionIntegratedSdfReader, contains application sensitivity profile dictionaries
    - experiment: RegionIntegratedSdfReader, contains experiment sensitivity profile dictionaries
    
    Returns
    -------
    - nuclide_wise_contributions: list of dict, list of dictionaries containing the contribution to the similarity parameter E
        for each nuclide
    - nuclide_reaction_wise_contributions: list of dict, list of dictionaries containing the contribution to the similarity"""

    # First, extract the sensitivity vectors for the application and experiment
    
    # Calculate |S_A| and |S_E| to normalize the E contributions properly

    application_vector = create_sensitivity_vector(application.get_sensitivity_profiles())
    experiment_vector = create_sensitivity_vector(experiment.get_sensitivity_profiles())

    application_norm = umath.sqrt(np.sum(application_vector**2))
    experiment_norm = umath.sqrt(np.sum(experiment_vector**2))

    # Now convert the application and experiment sdf's to dictionaries keyed by nuclide and reaction type
    application = application.convert_to_dict().sdf_data
    experiment  = experiment.convert_to_dict().sdf_data

    nuclide_wise_contributions = []
    nuclide_reaction_wise_contributions = []

    # All isotopes in the application and experiment
    all_isotopes = set(application.keys()).union(set(experiment.keys()))

    # Assuming each nuclide has the same reactions, choose an arbitrary isotope in the application to get the reactions
    arbitrary_nuclide = list(application.keys())[0]
    all_reactions = application[arbitrary_nuclide].keys()

    for isotope in all_isotopes:
        isotope_does_not_contribute = isotope not in application.keys() or isotope not in experiment.keys()
        if isotope_does_not_contribute:
            nuclide_wise_contributions.append({
                "isotope": isotope,
                "contribution": ufloat(0, 0)
            })
            nuclide_reaction_wise_contributions += [ {
                "isotope": isotope,
                "reaction_type": reaction_type,
                "contribution": ufloat(0, 0)
            } for reaction_type in all_reactions ]
            continue

        # For isotope-wise contribution, the sensitivity vector is all of the reaction sensitivities concatenated together
        application_vector = create_sensitivity_vector([ application[isotope][reaction]['sensitivities'] \
                                                        for reaction in all_reactions] )
        experiment_vector  = create_sensitivity_vector([ experiment[isotope][reaction]['sensitivities'] \
                                                        for reaction in all_reactions ])

        E_isotope_contribution = calculate_E_from_sensitivity_vecs(application_vector, experiment_vector, uncertainties='manual',
                                                                   application_norm=application_norm, \
                                                                    experiment_norm=experiment_norm)
        nuclide_wise_contributions.append({
            "isotope": isotope,
            "contribution": E_isotope_contribution
        })

        # For nuclide-reaction-wise contribution, we need to consider each reaction type
        nuclide_reaction_wise_contributions += \
            get_reaction_wise_E_contributions(application, experiment, isotope, all_reactions, \
                                              application_norm, experiment_norm)
        
    return nuclide_wise_contributions, nuclide_reaction_wise_contributions


def calculate_E_contributions(application_filenames: list, experiment_filenames: list):
    """Calculates the contributions to the similarity parameter E for each application with each available experiment 
    on a nuclide basis and on a nuclide-reaction basis
    
    Parameters
    ----------
    - application_filenames: list of str, paths to the application sdf files
    - experiment_filenames: list of str, paths to the experiment sdf files
    
    Returns
    -------
    - E_contributions_nuclide: unumpy.uarray of contributions to the similarity parameter E for each application with each
        experiment on a nuclide basis.
    - E_contributions_nuclide_reaction: unumpy.uarray of contributions to the similarity parameter E for each application with
        each experiment on a nuclide-reaction basis"""
    
    application_sdfs = [ RegionIntegratedSdfReader(filename) for filename in application_filenames ]
    experiment_sdfs  = [ RegionIntegratedSdfReader(filename) for filename in experiment_filenames ]
    
    # Initialize np object arrays to store the E contributions
    E_nuclide_wise          = np.empty( ( len(experiment_sdfs), len(application_sdfs) ), dtype=object )
    E_nuclide_reaction_wise = np.empty( ( len(experiment_sdfs), len(application_sdfs) ), dtype=object )

    for i, experiment in enumerate(experiment_sdfs):
        for j, application in enumerate(application_sdfs):
            nuclide_wise_contributions, nuclide_reaction_wise_contributions = \
                get_nuclide_and_reaction_wise_E_contributions(application, experiment)

            E_nuclide_wise[i, j] = nuclide_wise_contributions
            E_nuclide_reaction_wise[i, j] = nuclide_reaction_wise_contributions

    return E_nuclide_wise, E_nuclide_reaction_wise


# ========================
# Visualization Utilities
# ========================

class Plotter(ABC):
    @abstractmethod
    def create_plot(self, data, nested):
        pass

    @abstractmethod
    def add_to_subplot(self, fig, position):
        pass

    @abstractmethod
    def get_plot(self):
        pass

    @abstractmethod
    def style(self):
        pass


class BarPlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.fig, self.axs = plt.subplots()
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested):
        if nested:
            self.nested_barchart(contributions)
        else:
            self.barchart(contributions)

        self.style()

    def get_plot(self):
        return self.fig, self.axs
        
    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.axs, sharey=self.axs)
        
    def barchart(self, contributions):
        self.axs.bar(contributions.keys(), [contribution.n for contribution in contributions.values()],
            yerr=[contribution.s for contribution in contributions.values()], capsize=5, error_kw={'elinewidth': 0.5})

    def nested_barchart(self, contributions):
        # Colors for each reaction type
        num_reactions = len(next(iter(contributions.values())))
        cmap = plt.get_cmap('Set1')
        colors = cmap(np.linspace(0, 1, num_reactions))

        # Variables to hold the bar positions and labels
        indices = range(len(contributions))
        labels = list(contributions.keys())

        # Bottom offset for each stack
        bottoms_pos = [0] * len(contributions)
        bottoms_neg = [0] * len(contributions)

        color_index = 0
        for reaction in next(iter(contributions.values())).keys():
            values = [contributions[nuclide][reaction].n for nuclide in contributions]
            errs = [contributions[nuclide][reaction].s for nuclide in contributions]
            # Stacking positive values
            pos_values = [max(0, v) for v in values]
            neg_values = [min(0, v) for v in values]
            self.axs.bar(indices, pos_values, label=reaction, bottom=bottoms_pos, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            self.axs.bar(indices, neg_values, bottom=bottoms_neg, color=colors[color_index % len(colors)],
                    yerr=errs, capsize=5, error_kw={'capthick': 0.5})
            # Update the bottom positions
            bottoms_pos = [bottoms_pos[i] + pos_values[i] for i in range(len(bottoms_pos))]
            bottoms_neg = [bottoms_neg[i] + neg_values[i] for i in range(len(bottoms_neg))]
            color_index += 1

        # Adding 'effective' box with dashed border
        total_values = [sum(contributions[label][r].n for r in contributions[label]) for label in labels]
        for idx, val in zip(indices, total_values):
            self.axs.bar(idx, abs(val), bottom=0 if val > 0 else val, color='none', edgecolor='black', hatch='///', linewidth=0.5)

        self.axs.set_xticks(indices)
        self.axs.set_xticklabels(labels)
        self.axs.legend()

    def style(self):
        if self.plot_redundant:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.set_ylabel(f"Contribution to {self.index_name}")
        self.axs.set_xlabel("Isotope")
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)


class PiePlotter(Plotter):
    def __init__(self, integral_index_name, plot_redudant=False, **kwargs):
        self.fig, self.axs = plt.subplots()
        self.index_name = integral_index_name
        self.plot_redundant = plot_redudant
    
    def create_plot(self, contributions, nested):
        if nested:
            self.nested_pie_chart(contributions)
        else:
            self.pie_chart(contributions)

        self.style()

    def add_to_subplot(self, fig, position):
        return fig.add_subplot(position, sharex=self.ax, sharey=self.ax)

    def get_plot(self):
        return self.fig, self.axs

    def nested_pie_chart(self, contributions):
        # Create a nested ring chart
        num_reactions = len(next(iter(contributions.values())))
        nuclide_colors = plt.get_cmap('rainbow')(np.linspace(0, 1, len(contributions.keys())))
        nuclide_totals = { nuclide: sum(contribution.n for contribution in contributions[nuclide].values()) \
                        for nuclide in contributions }
        nuclide_labels = list(nuclide_totals.keys())

        # Now, deal with negative values

        nuclides_with_opposite_sign_contributions = []
        for nuclide, contribution in contributions.items():
            contribution_values = [contribution[reaction].n for reaction in contribution]
            if not (all(v >= 0 for v in contribution_values) or all(v <= 0 for v in contribution_values)):
                nuclides_with_opposite_sign_contributions.append(nuclide)
            
        # For nuclides with opposite sign contributions, we distinguish the positive and negative contributions
        # by coloring some of the inner ring a lighter color to indicate the negative contributions in the outer ring
        wedge_widths = list(nuclide_totals.values())
        inner_wedge_hatches = [None] * len(wedge_widths)

        def blend_colors(color1, color2, alpha):
            return np.array( [ alpha * c1 + (1 - alpha) * c2 for c1, c2 in zip(color1, color2 ) ] )

        if len(nuclides_with_opposite_sign_contributions) > 0:
            for nuclide in nuclides_with_opposite_sign_contributions:
                # First, determine the fraction of the contributions that are opposite (in sign) to the total
                total_sign = np.sign(nuclide_totals[nuclide])
                
                # Now, we want to plot the "lost" wedge width in white, i.e. the width lost from cancellations between the
                # positive and negative contributions. This will be colored a lighter color. The absolute sum of the
                # contributions represents the wedge width if there were no cancellations, so the total wedge width
                # minus the absolute sum of the contributions is "lost" wedge width.

                absolute_sum_of_contributions = sum(np.abs(contribution.n) for contribution in contributions[nuclide].values())
                
                # NOTE the sign function is needed to handle the case when the nuclide total is negative
                lost_wedge_width = absolute_sum_of_contributions - total_sign * nuclide_totals[nuclide]

                # Now, insert the lost wedge width into the wedge widths list right after the nuclide
                nuclide_index = list(nuclide_totals.keys()).index(nuclide)
                wedge_widths.insert(nuclide_index + 1, lost_wedge_width)
                nuclide_labels.insert(nuclide_index + 1, '')
                
                # The color of the lost wedge width will be a blend of the nuclide color and white
                white_color = np.array([1, 1, 1, 1])
                opacity = 0.8
                blended_color = blend_colors(white_color, nuclide_colors[nuclide_index], opacity)
                nuclide_colors = np.insert(nuclide_colors, nuclide_index + 1, blended_color, axis=0)
                
                # Add hatches to the negative total sum wedge
                if nuclide_totals[nuclide] < 0:
                    inner_wedge_hatches[nuclide_index] = '//'

        # Now make everything positive for the pie chart
        wedge_widths = np.abs(wedge_widths)

        # Plot the inner ring for nuclide totals
        inner_ring, _ = self.axs.pie(wedge_widths, radius=0.7, labels=nuclide_labels, \
                                colors=nuclide_colors, labeldistance=0.6, textprops={'fontsize': 8}, \
                                    wedgeprops=dict(width=0.3, edgecolor='w'))

        # Add hatches to the negative total sum wedges
        for wedge, hatch in zip(inner_ring, inner_wedge_hatches):
            if hatch:
                wedge.set_hatch(hatch)

        # Get colors for reactions from the "rainbow" colormap
        reaction_colors = plt.get_cmap('Set1')(np.linspace(0, 1, num_reactions))

        # Plot the outer ring for reaction-specific contributions
        outer_labels = []
        outer_colors = []
        outer_sizes = []
        outer_hatches = []
        for i, (nuclide, reactions) in enumerate(contributions.items()):
            for j, (reaction, contribution) in enumerate(list(reactions.items())):
                outer_labels.append(f"{nuclide} - {reaction}")
                
                outer_colors.append(reaction_colors[j])
                outer_sizes.append(np.abs(contribution.n))
                
                if contribution.n < 0:
                    outer_hatches.append('//')
                else:
                    outer_hatches.append(None)

        outer_ring, _ = self.axs.pie(outer_sizes, radius=1, labels=outer_labels, labeldistance=0.9, colors=outer_colors, \
                textprops={'fontsize': 6}, startangle=inner_ring[0].theta1, counterclock=True, \
                    wedgeprops=dict(width=0.3, edgecolor='w'))

        # Add hatches to the negative contribution wedges
        for wedge, hatch in zip(outer_ring, outer_hatches):
            if hatch:
                wedge.set_hatch(hatch)
        
    def pie_chart(self, contributions):
        labels = list(contributions.keys())
        values = [abs(contributions[key].n) for key in labels]

        # Determining hatching patterns: empty string for positive, cross-hatch for negative
        hatches = ['//' if contributions[key].n < 0 else '' for key in labels]

        # Creating the pie chart
        wedges, _ = self.axs.pie(values, labels=labels, startangle=90)

        # Applying hatching patterns to the wedges
        for wedge, hatch in zip(wedges, hatches):
            wedge.set_hatch(hatch)

    def style(self):
        if self.plot_redundant:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.axs.grid(True, which='both', axis='y', color='gray', linestyle='-', linewidth=0.5)
        self.axs.set_title(title_text)


class InteractiveLegend:
    def __init__(self, fig, df):
        """Return a flask webapp that will display an interactive legend for the sunburst chart"""
        self.fig = fig

        self.app = Flask(__name__)

        @self.app.route('/shutdown', methods=['POST'])
        def shutdown():
            os.kill(os.getpid(), signal.SIGINT)  # Send the SIGINT signal to the current process
            return 'Server shutting down...'

        @self.app.route('/')
        def show_sunburst():
            # Extract root nodes (nodes without parents)
            root_nodes = df[df['parents'] == '']

            # Generate legend HTML with a title
            legend_html = '<div id="legend" style="margin-left: 20px; border: 2px solid black; padding: 10px;"><h3 style="margin-top: 0;">Legend</h3>\n'
            for _, row in root_nodes.iterrows():
                legend_html += f'    <div class="legend-item" style="cursor: pointer; margin-bottom: 5px;" data-target="{row["ids"]}">{row["ids"]}: {row["values"]}</div>\n'
            legend_html += '</div>\n'

            # JavaScript for interactivity and shutdown
            script_html = """
            <script>
            window.addEventListener('beforeunload', (event) => {
                navigator.sendBeacon('/shutdown');
            });
            document.addEventListener('DOMContentLoaded', function () {
                const legendItems = document.querySelectorAll('.legend-item');
                legendItems.forEach(item => {
                    item.addEventListener('mouseenter', function() {
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {
                                path.style.opacity = 0.5; // Highlight by changing opacity
                            }
                        });
                    });
                    item.addEventListener('mouseleave', function() {
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            path.style.opacity = 1; // Reset opacity
                        });
                    });
                    item.addEventListener('click', function() {
                        const target = this.getAttribute('data-target');
                        const paths = document.querySelectorAll('path.surface');
                        paths.forEach(path => {
                            const labelText = path.nextElementSibling ? path.nextElementSibling.textContent : "";
                            if (labelText.includes(target)) {
                                path.dispatchEvent(new MouseEvent('click', { 'view': window, 'bubbles': true, 'cancelable': true }));
                            }
                        });
                    });
                });
            });
            </script>
            """

            # Save the chart with interactivity and layout adjustments
            fig_html = self.fig.to_html(full_html=False, include_plotlyjs='cdn')
            full_html = f"""
            <!DOCTYPE html>
            <html>
            <head>
            <title>Interactive Sunburst Chart</title>
            <style>
                body, html {{
                    height: 100%;
                    margin: 0;
                    font-family: Arial, sans-serif;
                }}
                #container {{
                    display: flex;
                    justify-content: center;
                    align-items: center;
                    height: 100%;
                    width: 100%;
                }}
                #chart {{
                    min-width: 0;    /* Prevent flex item from overflowing its container */
                }}
                #legend {{
                    flex: 0 1 auto;  /* Do not grow, allow shrink */
                    margin-left: 20px;
                    padding: 10px;
                    border: 2px solid black;
                    max-height: 100%; /* Make sure it does not overflow vertically */
                    overflow: auto;   /* Enable scrolling if content is too large */
                }}
            </style>
            </head>
            <body>
            <div id="container">
                <div id="chart">{fig_html}</div>
                {legend_html}
            </div>
            {script_html}
            </body>
            </html>
            """

            return render_template_string(full_html)
        
    def open_browser(self):
        # webbrowser.open_new("http://localhost:5000/")
        pass

    def show(self):
        # Suppress Flask's startup and runtime messages by redirecting them to dev null
        log = open(os.devnull, 'w')
        # sys.stdout = log
        # sys.stderr = log

        threading.Timer(1, self.open_browser).start()
        self.app.run()


class InteractivePiePlotter(Plotter):
    def __init__(self, integral_index_name, plot_redundant=False, **kwargs):
        self.fig = make_subplots()

        # Check if the user wants an interactive legend
        if 'interactive_legend' in kwargs.keys():
            self.interactive_legend = kwargs['interactive_legend']
        else:
            self.interactive_legend = True
        
        self.index_name = integral_index_name
        self.plot_redundant = plot_redundant

    def create_plot(self, contributions, nested=True):
        # Prepare data for the sunburst chart
        if nested:
            df = self._create_nested_sunburst_data(contributions)
        else:
            df = self._create_sunburst_data(contributions)
        
        # Create a sunburst chart
        self.fig = px.sunburst(
            data_frame=df,
            names='labels',
            parents='parents',
            ids='ids',
            values='normalized_values',
            custom_data=['values', 'uncertainties']
        )

        # Update hovertemplate with correct syntax
        self.fig.update_traces(
            hovertemplate=(
                "<b>%{label}</b><br>"
                "Value: %{customdata[0]:1.4E} +/- %{customdata[1]:1.4E}"  # Corrected format specifiers
                "<extra></extra>"  # This hides the trace info
            )
        )

        # Now style the plot
        self.style()

        if self.interactive_legend:
            self.fig = InteractiveLegend(self.fig, df)


    
    def add_to_subplot(self, fig, position):
        if self.interactive_legend:
            raise ValueError("Interactive legends are not supported when adding to a subplot")
        else:
            for trace in self.fig.data:
                fig.add_trace(trace, row=position[0], col=position[1])
            return fig

    def get_plot(self):
        return self.fig

    def _create_sunburst_data(self, contributions):
        data = {
            'labels': [], 
            'ids': [], 
            'parents': [], 
            'values': [], 
            'uncertainties': [],
            'normalized_values': [],
            'nuclide': []
        }

        abs_sum_of_nuclide_totals = sum( abs(contribution.n) for contribution in contributions.values())

        for nuclide, nuclide_total in contributions.items():
            # Caclulate the nuclide total, and the positive and negative contributions
            norm_nuclide_total = abs(nuclide_total) / abs_sum_of_nuclide_totals

            # Add the nuclide as a parent
            data['labels'].append(nuclide)
            data['ids'].append(nuclide)
            data['parents'].append('')
            data['values'].append(nuclide_total.n)
            data['uncertainties'].append(nuclide_total.s)
            data['normalized_values'].append(norm_nuclide_total.n)
            data['nuclide'].append(nuclide)

        return pd.DataFrame(data)

    def _create_nested_sunburst_data(self, contributions):
        data = {
            'labels': [], 
            'ids': [], 
            'parents': [], 
            'values': [], 
            'uncertainties': [],
            'normalized_values': [],
            'nuclide': []
        }

        abs_sum_of_nuclide_totals = sum(sum(abs(contribution.n) for contribution in reactions.values()) \
                                    for reactions in contributions.values())

        for nuclide, reactions in contributions.items():
            # Caclulate the nuclide total, and the positive and negative contributions
            nuclide_total = sum(contribution for contribution in reactions.values())
            norm_nuclide_total = abs(nuclide_total) / abs_sum_of_nuclide_totals

            positive_contributions = { reaction: contribution for reaction, contribution in reactions.items() \
                                      if contribution.n >= 0 }
            negative_contributions = { reaction: contribution for reaction, contribution in reactions.items() \
                                      if contribution.n < 0 }
            positive_total = sum(contribution for contribution in positive_contributions.values())
            negative_total = sum(contribution for contribution in negative_contributions.values())

            # Add the nuclide as a parent
            data['labels'].append(nuclide)
            data['ids'].append(nuclide)
            data['parents'].append('')
            data['values'].append(nuclide_total.n)
            data['uncertainties'].append(nuclide_total.s)
            data['normalized_values'].append(norm_nuclide_total.n)
            data['nuclide'].append(nuclide)
    
            # --------------------------------------------------------
            # Add the positive and negative contributions as children
            # --------------------------------------------------------

            # Normalize the contributions by the absolute value of the nuclide total 
            absolute_sum = positive_total + abs(negative_total)
            normalization_factor = abs(norm_nuclide_total) / absolute_sum

            # Positive contributions
            if positive_total != 0:
                norm_positive_total = positive_total * normalization_factor
                data['labels'].append('Positive')
                data['ids'].append(f"{nuclide}-Positive")
                data['parents'].append(nuclide)
                data['values'].append(positive_total.n)
                data['uncertainties'].append(positive_total.s)
                data['normalized_values'].append( norm_positive_total.n )
                data['nuclide'].append(nuclide)

            # Negative contributions
            if negative_total != 0:
                norm_negative_total = abs(negative_total) * normalization_factor
                data['labels'].append('Negative')
                data['ids'].append(f"{nuclide}-Negative")
                data['parents'].append(nuclide)
                data['values'].append(negative_total.n)
                data['uncertainties'].append(negative_total.s)
                data['normalized_values'].append( norm_negative_total.n )
                data['nuclide'].append(nuclide)

            # -------------------------------
            # Add the reaction contributions
            # -------------------------------
            # NOTE: Plotly express apparently has issues dealing with small numbers, so unless the contribution is
            # multiplied by a sufficiently large scale factor, the data won't be displayed correctly
            scale_factor = 10000
            for reaction, contribution in positive_contributions.items():
                # Now normalize contributions so they sum to the "normalized_positive_total"
                normalization_factor = norm_positive_total / positive_total
                norm_reaction_contribution = contribution.n * normalization_factor
                
                if contribution.n != 0:
                    data['labels'].append(reaction)
                    data['ids'].append(f"{nuclide}-{reaction}")
                    data['parents'].append(f"{nuclide}-Positive")
                    data['values'].append(contribution.n)
                    data['uncertainties'].append(contribution.s)
                    data['normalized_values'].append(scale_factor*norm_reaction_contribution.n)
                    data['nuclide'].append(nuclide)

            for reaction, contribution in negative_contributions.items():
                # Now normalize contributions so they sum to the "normalized_negative_total"
                normalization_factor = norm_negative_total / abs(negative_total)
                norm_reaction_contribution = abs(contribution.n) * normalization_factor

                if contribution.n != 0:
                    data['labels'].append(reaction)
                    data['ids'].append(f"{nuclide}-{reaction}")
                    data['parents'].append(f"{nuclide}-Negative")
                    data['values'].append(contribution.n)
                    data['uncertainties'].append(contribution.s)
                    data['normalized_values'].append(scale_factor*norm_reaction_contribution.n)
                    data['nuclide'].append(nuclide)


        return pd.DataFrame(data)

    def style(self):
        if self.plot_redundant:
            title_text = f'Contributions to {self.index_name} (including redundant/irrelvant reactions)'
        else:
            title_text = f'Contributions to {self.index_name}'
        self.fig.update_layout(title_text=title_text, title_x=0.5)  # 'title_x=0.5' centers the title


def plot_contributions(contributions, plot_type='bar', integral_index_name='E', plot_redundant_reactions=False, **kwargs):
    """Plots the contributions to an arbitrary similarity parameter for a single experiment application pair
    
    Parameters
    ----------
    - contributions: list of dict, list of dictionaries containing the contributions to the similarity parameter for each
        nuclide or nuclide-reaction pair
    - plot_type: str, type of plot to create. Default is 'bar' which creates a bar plot. Other option is 'pie' which creates
        a pie chart
    - integral_index_name: str, name of the integral index being plotted. Default is 'E'
    - plot_redundant_reactions: bool, whether to plot redundant reactions (or irrelevant reactions) when considering
        nuclide-reaction-wise contributions. Default is False
    - kwargs: additional keyword arguments to pass to the plotting function"""

    # Determine if the contributions are nuclide-wise or nuclide-reaction-wise
    if 'reaction_type' in contributions[0]:
        # Nuclide-reaction-wise contributions
        nested_plot = True # Nested plot by nuclide then by reaction type

        # Create a dictionary of contributions keyed by isotope then by reaction type
        contributions = isotope_reaction_list_to_nested_dict(contributions, 'contribution')

        # If viewing nuclide-reaction wise contributions, it's important (at least for the visualizations in this function)
        # that if viewing the true contributions to the nuclide total, that redundant interactions (e.g. capture and fission
        # + (n, g)) and irrelevant interactions (e.g. chi and nubar) are not plotted.

        if not plot_redundant_reactions:
            # Remove redundant interactions
            redundant_interactions = ['chi', 'capture', 'nubar', 'total']
            contributions = { isotope: { reaction: contributions[isotope][reaction] for reaction in contributions[isotope] \
                                if reaction not in redundant_interactions } for isotope in contributions }
    else:
        # Nuclide-wise contributions
        nested_plot = False
        contributions = { contribution['isotope']: contribution['contribution'] for contribution in contributions }

    plotters = {
        'bar': BarPlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'pie': PiePlotter(integral_index_name, plot_redundant_reactions, **kwargs),
        'interactive_pie': InteractivePiePlotter(integral_index_name, plot_redundant_reactions, **kwargs)
    }
    
    # Get the requested plotter
    plotter = plotters.get(plot_type)
    if plotter is None:
        raise ValueError("Unsupported plot type")

    # Create the plot and style it
    plotter.create_plot(contributions, nested_plot)

    return plotter.get_plot()

# =====
# Misc
# =====

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