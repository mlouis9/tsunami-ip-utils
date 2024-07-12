import re

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
        matches = re.findall(r'\d+', isotope)
        return int(matches[0] if matches else -1) # Return -1 if no atomic number is found, these isotopes will be sorted last
                                                  # Should only be applicable to carbon in ENDF 7.1
    
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

def filter_redundant_reactions(data_dict, redundant_reactions=['chi', 'capture', 'nubar', 'total']):
    """Filters out redundant reactions from a nested isotope-reaction dictionary
    
    Parameters
    ----------
    - data_dict: dict, nested dictionary containing isotope-reaction pairs
    - redundant_reactions: list of str, list of reactions to filter out"""
    return { isotope: { reaction: data_dict[isotope][reaction] for reaction in data_dict[isotope] \
                        if reaction not in redundant_reactions } for isotope in data_dict }

def filter_by_nuclie_reaction_dict(data_dict, nuclide_reactions):
    """Filters out isotopes that are not in the nuclide_reactions dictionary
    
    Parameters
    ----------
    - data_dict: dict, nested dictionary containing isotope-reaction pairs
    - nuclide_reactions: dict, dictionary containing isotopes and their reactions"""
    return {nuclide: {reaction: xs for reaction, xs in reactions.items() if reaction in nuclide_reactions[nuclide]} \
                        for nuclide, reactions in data_dict.items() if nuclide in nuclide_reactions.keys()}