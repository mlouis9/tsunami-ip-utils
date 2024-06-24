from uncertainties import ufloat, umath, unumpy
import numpy as np
from tsunami_ip_utils.readers import RegionIntegratedSdfReader, read_uncertainty_contributions
from tsunami_ip_utils.error import unit_vector_uncertainty_propagation, dot_product_uncertainty_propagation
from copy import deepcopy

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

    application_sdfs = [ RegionIntegratedSdfReader(filename).convert_to_dict() for filename in application_filenames ]
    experiment_sdfs  = [ RegionIntegratedSdfReader(filename).convert_to_dict() for filename in experiment_filenames ]

    # Create a matrix to store the similarity parameter E for each application with each experiment
    E_vals = unumpy.umatrix(np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ), \
                            np.zeros( ( len(experiment_sdfs), len(application_sdfs) ) ))
    
    # Now calculate the similarity parameter E for each application with each experiment
    for i, experiment in enumerate(experiment_sdfs):
        for j, application in enumerate(application_sdfs):
            # Now add missing data to the application and experiment dictionaries
            all_isotopes = set(application.sdf_data.keys()).union(set(experiment.sdf_data.keys()))
            add_missing_reactions_and_nuclides(application.sdf_data, experiment.sdf_data, all_isotopes)
            
            # Sometimes the application and experiment dictionaries have different orders, which causes their sensitivity
            # profiles and hence sensitivity vectors to be different. To fix this, we need to sort the dictionaries

            # Sort the application by the experiment keys
            application.sdf_data = { isotope: { reaction: application.sdf_data[isotope][reaction] \
                                               for reaction in experiment.sdf_data[isotope] } \
                                                for isotope in experiment.sdf_data.keys() }

            application_profiles = application.get_sensitivity_profiles()
            experiment_profiles  = experiment.get_sensitivity_profiles()

            application_vector = create_sensitivity_vector(application_profiles)
            experiment_vector  = create_sensitivity_vector(experiment_profiles)

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

def add_missing_reactions_and_nuclides(application, experiment, all_isotopes, mode='sdfs'):
    """Add missing reactions and nuclides to the application and experiment dictionaries with an sdf profile of all zeros.
    NOTE: Since dictionaries are passed by reference, this function does not return anything, but modifies the
    application and experiment dictionaries in place.
    
    Parameters
    ----------
    - application: dict, dictionary of application sensitivity profiles
    - experiment: dict, dictionary of experiment sensitivity profiles
    - all_isotopes: list of str, list of all isotopes in the application and experiment dictionaries
    - mode: str, the mode to use for adding missing reactions and nuclides. Default is 'sdfs' which adds sdf profiles to missing
        reactions and nuclides. If set to 'contribution', then the contributions to the similarity parameter E are set to zero
        for missing nuclides and reactions
    
    Returns
    -------
    - all_isotopes: set of str, set of all isotopes in the application and experiment dictionaries"""
    application_reactions = set([ key for isotope in application.keys() for key in application[isotope].keys() ])
    experiment_reactions = set([ key for isotope in experiment.keys() for key in experiment[isotope].keys() ])
    all_reactions = application_reactions.union(experiment_reactions)

    # Now, for any reactions that are in experiment but not application (and vice versa), they need to be added with an sdf
    # profile of all zeros

    # Get an arbitrary sdf profile to get the shape from
    first_application_nuclide = list(application.keys())[0]
    first_application_reaction = list(application[first_application_nuclide].keys())[0]
    if mode == 'sdfs':
        zero_data = {
            'sensitivities': application[first_application_nuclide][first_application_reaction]['sensitivities']
        }
    elif mode == 'contribution':
        zero_data = ufloat(0,0)

    for isotope in application.keys():
        # If reaction is missing for this isotope, add it with an sdf profile of all zeros
        for reaction in all_reactions:
            if reaction not in application[isotope].keys():
                application[isotope][reaction] = deepcopy(zero_data)
    for isotope in experiment.keys():
        # If reaction is missing for this isotope, add it with an sdf profile of all zeros
        for reaction in all_reactions:
            if reaction not in experiment[isotope].keys():
                experiment[isotope][reaction] = deepcopy(zero_data)

    # Now zero out nuclides that are not in the application or experiment
    for isotope in all_isotopes:
        if isotope not in application.keys():
            application[isotope] = { reaction: deepcopy(zero_data) for reaction in all_reactions }
            
        if isotope not in experiment.keys():
            experiment[isotope] = { reaction: deepcopy(zero_data) for reaction in all_reactions }

    return all_reactions

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

    # Since different nuclides can have different reactions, we need to consider all reactions for each nuclide (e.g. 
    # fissile isotopes will have fission reactions, while non-fissile isotopes will not)
    all_reactions = add_missing_reactions_and_nuclides(application, experiment, all_isotopes)

    for isotope in all_isotopes:

        # For isotope-wise contribution, the sensitivity vector is all of the reaction sensitivities concatenated together
        application_vector = create_sensitivity_vector([ application[isotope][reaction]['sensitivities'] \
                                                        for reaction in all_reactions ] )
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


def calculate_uncertainty_contributions(application_filenames: list, experiment_filenames: list):
    """Calculates the contributions to the uncertainty in k (i.e. dk/k) for each application with each available experiment
    on a nuclide basis and on a nuclide-reaction basis

    Parameters
    ----------
    - application_filenames: list of str, paths to the application sdf files
    - experiment_filenames: list of str, paths to the experiment sdf files

    Returns
    -------
    - uncertainty_contributions_nuclide: dict of unumpy.uarray of contributions to the uncertainty in k for each application and 
        each experiment on a nuclide basis. Keyed by 'application' and 'experiment'
    - uncertainty_contributions_nuclide_reaction: dict of unumpy.uarray of contributions to the uncertainty in k
        for each application and each experiment on a nuclide-reaction basis. Keyed by 'application' and 'experiment'"""
    
    dk_over_k_nuclide_wise = {
        'application': np.empty( len(application_filenames), dtype=object ),
        'experiment': np.empty( len(experiment_filenames), dtype=object )    
    }
    dk_over_k_nuclide_reaction_wise = {
        'application': np.empty( len(application_filenames), dtype=object ),
        'experiment': np.empty( len(experiment_filenames), dtype=object )    
    }

    for i, application_filename in enumerate(application_filenames):
        dk_over_k_nuclide_wise['application'][i], dk_over_k_nuclide_reaction_wise['application'][i] = \
            read_uncertainty_contributions(application_filename)

    for i, experiment_filename in enumerate(experiment_filenames):
        dk_over_k_nuclide_wise['experiment'][i], dk_over_k_nuclide_reaction_wise['experiment'][i] = \
            read_uncertainty_contributions(experiment_filename)
            
    return dk_over_k_nuclide_wise, dk_over_k_nuclide_reaction_wise