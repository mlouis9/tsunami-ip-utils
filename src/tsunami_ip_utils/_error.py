from uncertainties import unumpy
import numpy as np

def _unit_vector_uncertainty_propagation(vector: unumpy.uarray) -> np.ndarray:
    """Performs error propagation for the components of a vector :math:`\\boldsymbol{v}` that is normalized to a unit vector 
    via :math:`\\hat{\\boldsymbol{v}} = \\frac{\\boldsymbol{v}}{\\lVert\\boldsymbol{v}\\rVert}`. To calculate the uncertainty 
    of :math:`\\hat{\\boldsymbol{v}}\\equiv \\boldsymbol{u}`, recall from linear error propagation theory 
    (see `this <https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Simplification>`_ equation)

    .. math::
        \\sigma_{u_i} = \\sqrt{ \\sum_j \\left( \\frac{\\partial u_i}{\\partial v_j} \\sigma_{v_j}  \\right)^2 }

    since
    
    .. math::
        u_i = \\frac{v_i}{\\sqrt{\\sum_j v_j^2}}
        
    we have for :math:`i=j`:

    .. math::
        \\frac{\\partial u_i}{\\partial v_j} &= 
        \\frac{1}{\\sqrt{\\sum_j v_j^2}} - \\frac{1}{2}\\frac{v_i\\cdot 2v_i}{\\left(\\sum_j v_j^2\\right)^{\\frac{3}{2}}}\\\\
        &=\\frac{1}{\\lVert \\boldsymbol{v}\\rVert} - \\frac{v_i^2}{\\lVert \\boldsymbol{v}\\rVert^3}\\\\
        &=\\frac{\\lVert \\boldsymbol{v}\\rVert^2 - v_i^2}{\\lVert \\boldsymbol{v}\\rVert^3}

    For :math:`i\\neq j`:

    .. math::
        \\frac{\\partial u_i}{\\partial v_j} = -\\frac{v_iv_j}{\\lVert\\boldsymbol{v}\\rVert^3}

    Examples
    --------
    >>> unit_vector_uncertainty_propagation(unumpy.uarray([1, 2, 3], [0.1, 0.2, 0.3]))
    array([0.03113499, 0.05150788, 0.03711537])

    >>> unit_vector_uncertainty_propagation(unumpy.uarray([0, 0, 0], [0.1, 0.2, 0.3]))
    array([0., 0., 0.])

    >>> unit_vector_uncertainty_propagation(unumpy.uarray([0, 0, 0], [0.0, 0.0, 0.0]))
    array([0., 0., 0.])
    
    Parameters
    ----------
    vector
        The vector :math:`\\boldsymbol{v}` (with uncertainties) whose normalization the error propagation is
        calculated for.
    
    Returns
    -------
        Uncertainties of the unit vector components."""

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

def _dot_product_uncertainty_propagation(vect1: unumpy.uarray, vect2: unumpy.uarray) -> float:
    """Calculates the uncertainty in the dot product of two vectors with uncertainties
    
    Parameters
    ----------
    vect1
        First vector in dot product.
    vect2
        Second vector in dot product.

    Returns
    -------
        Uncertainty in the dot product of the two vectors."""

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
