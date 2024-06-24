from uncertainties import unumpy
import numpy as np

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
    if vector_norm != 0:
        derivative_matrix = -unumpy.nominal_values(vector)[:, np.newaxis] * unumpy.nominal_values(vector)[np.newaxis, :] / vector_norm**3
        np.fill_diagonal(derivative_matrix, (vector_norm**2 - unumpy.nominal_values(vector)**2) / vector_norm**3)
    else:
        derivative_matrix = np.zeros((len(vector), len(vector)))

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
