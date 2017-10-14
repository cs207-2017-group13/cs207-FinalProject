import numpy as np


def calculate_progress_rate(concentrations, reactant_coeffs, reaction_coeffs):
    """For a system of reactions, determine reaction progress rate(s).

    Parameters
    ----------
    concentrations : array_like, 1D
        An ``N``-dimensional vector of the concentrations of the
        reacting species.
    reactant_coeffs : array_like
        An ``NxM`` matrix representing stoichiometric coefficients for
        the reactants in the system of reactions, where ``N``is the
        number of reactants and ``M`` is the number of elementary
        reactions.
    reaction_coeffs : array_like, 1D
        An ``M``-dimensinal vector of the reaction rate coefficients,
        each corresponding to one of the elementary reactions in the
        system of reactions.

    Returns
    -------
    progress_rates : array_like, 1D
        An ``M``-dimensional vector of progress rates, one for each of
        the elementary reactions in the system of reactions.

    Raises
    ------
    ValueError
        Non-numeric values in function arguments.
    TypeError
        Function arguments cannot be converted to numpy arrays.
    AssertionError
        Incompatible dimensions between arguments
    IndexError
        Incompatible dimensions between arguments

    Examples
    --------
    >>> calculate_progress_rate([1., 2., 3.], [2, 1, 0], 10)
    array([ 20.])

    >>> reactant_coeffs = np.rollaxis(np.array([[1,2,0],[2,0,2]]), -1)
    >>> calculate_progress_rate([1., 2., 1.], reactant_coeffs, [10,10])
    array([ 40.,  10.])

    """
    # It's easier to convert to np arrays first.
    # Will also catch if data is not numeric
    try:
        concentrations = np.atleast_1d(concentrations).astype(float)
        reactant_coeffs = np.atleast_1d(reactant_coeffs).astype(int)
        if reactant_coeffs.ndim == 1:
            reactant_coeffs = reactant_coeffs.reshape(len(reactant_coeffs), 1)
        reaction_coeffs = np.atleast_1d(reaction_coeffs).astype(float)
    except Exception as err:
        raise err
    
    try:
        assert len(concentrations.shape) == 1
        assert len(reactant_coeffs.shape) == 2
        assert len(reaction_coeffs.shape) == 1
    except Exception as err:
        raise err

    n_reactants = len(concentrations)
    n_reactions = len(reaction_coeffs)
    if reactant_coeffs.shape != (n_reactants, n_reactions):
        raise IndexError("Dimensions of reactant coefficients array and "
                         "numbers of reactants and reactions do not agree.")

    progress_rates = np.zeros(len(reaction_coeffs))
    for j, (rate_constant, jth_reaction_reactant_coeffs) in enumerate(zip(
            reaction_coeffs, np.rollaxis(reactant_coeffs,-1))):
        progress_rates[j] = rate_constant
        for (concentration, reactant_coeff) in zip(
                concentrations, jth_reaction_reactant_coeffs):
            progress_rates[j] *= concentration ** reactant_coeff

    return progress_rates


def calculate_reaction_rate(concentrations, reactant_coeffs, product_coeffs,
                            reaction_coeffs):
    """For a system of reactions, determine reaction rate(s).

    Parameters
    ----------
    concentrations : array_like, 1D
        An ``N``-dimensional vector of the concentrations of the
        reacting species.
    reactant_coeffs : array_like
        An ``NxM`` matrix representing stoichiometric coefficients for
        the reactants in the system of reactions, where ``N``is the
        number of reacting species and ``M`` is the number of
        elementary reactions. Each column in `reactant_coeffs`
        corresponds to an elementary reaction in the system of
        reactions.
    product_coeffs : array_like
        An ``NxM`` matrix representing stoichiometric coefficients for
        the products in the system of reactions, where ``N``is the
        number of reacting species and ``M`` is the number of
        elementary reactions. Each column in `reactant_coeffs`
        corresponds to an elementary reaction in the system of
        reactions.
    reaction_coeffs : array_like, 1D
        An ``M``-dimensinal vector of the reaction rate coefficients,
        each corresponding to one of the elementary reactions in the
        system of reactions.

    Returns
    -------
    reaction_rates : array_like, 1D
        An ``N``-dimensional vector of reaction rates, one for each of
        the reactants in the system of reactions.

    Raises
    ------
    ValueError
        Non-numeric values in function arguments.
    TypeError
        Function arguments cannot be converted to numpy arrays.
    AssertionError
        Incompatible dimensions between arguments
    IndexError
        Incompatible dimensions between arguments

    Examples
    --------
    >>> reactant_coeffs = np.rollaxis(np.array([[1,2,0], [0,0,2]]), -1)
    >>> product_coeffs = np.rollaxis(np.array([[0,0,1], [1,2,0]]), -1)
    >>> calculate_reaction_rate([1., 2., 1.], reactant_coeffs,
    ... product_coeffs, [10,10])
    array([-30., -60.,  20.])

    """
    # It's easier to convert to np arrays first.
    # Will also catch if data is not numeric
    try:
        concentrations = np.atleast_1d(concentrations).astype(float)
        reactant_coeffs = np.atleast_1d(reactant_coeffs).astype(int)
        if reactant_coeffs.ndim == 1:
            reactant_coeffs = reactant_coeffs.reshape(len(reactant_coeffs), 1)
        product_coeffs = np.atleast_1d(product_coeffs).astype(int)
        if product_coeffs.ndim == 1:
            product_coeffs = product_coeffs.reshape(len(product_coeffs), 1)
        reaction_coeffs = np.atleast_1d(reaction_coeffs).astype(float)
    except Exception as err:
        raise err

    try:
        progress_rates = calculate_progress_rate(concentrations,
                                                 reactant_coeffs,
                                                 reaction_coeffs)
    except Exception as err:
        raise err

    # `progress_rate` should have checked that dimensions of
    #   `reactant_coeffs` agrees with `concentrations` and
    #   `reaction_coeffs`
    if reactant_coeffs.shape != product_coeffs.shape:
        raise IndexError("Dimensions of reactant and product stoichiometric "
                         "coefficient arrays do not agree.")
                                  
    net_coeffs = product_coeffs - reactant_coeffs
    reaction_rates = np.zeros(len(concentrations))
    for i, ith_species_net_coeffs in enumerate(net_coeffs):
        for net_coeff, progress_rate in zip(ith_species_net_coeffs,
                                            progress_rates):
            reaction_rates[i] += net_coeff * progress_rate
    
    return reaction_rates


import xml.etree.ElementTree as ET
import numpy as np
import numbers

IDEAL_GAS_CONSTANT = 8.314


def constant_rate(k):
    """Return a constant reaction rate coefficient.

    In zeroth-order reactions, k = constant.

    Parameters
    ----------
    k
        Input rate coefficient

    Returns
    -------
    k
        The input rate coefficient unchanged.

    Notes
    -----
    Although it would be sensible if input is numeric, no exceptions
    will be raised if this is not the case.

    Examples
    --------
    >>> constant_rate(2.3)
    2.3
    """
    # TODO: no negative values
    return k


def arrhenius_rate(A, E, T, b=0, R=IDEAL_GAS_CONSTANT):
    """Return a reaction rate coefficient according to the Arrhenius equation.

    The Arrhenius equation relates the rate constant, k, of a chemical
    reaction to parameters A (pre-exponential factor), E (activation
    energy), T (absolute temperature), and b (exponential indicating
    temperature-dependence of pre-exponential factor)::

        k = A T^b exp[ -E / (RT) ]

    When ``b = 0``, the above formula corresponds to the canonical
    Arrhenius equation.  A nonzero value for b gives the modified
    Arrhenius equation.

    Parameters
    ----------
    A : scalar
        Pre-exponential factor. ``A > 0``
    E : array_like
        Activation energy. If `E` and `T` are both arrays, then their
        dimensions must match.
    T : array_like
        Absolute temperature. ``T > 0``. If `E` and `T` are both
        arrays, then their dimensions must match.
    b : scalar, optional
        Exponent representing temperature-dependence of
        pre-exponential factor. A value of 0 corresponds to the
        canonical Arrhenius equation.

    Returns
    -------
    k : scalar
        Reaction rate coefficient calculated according to function
        arguments.

    Raises
    ------
    TypeError : not a number
        Arguments must be numeric types.
    ValueError : Require > 0.
        `A` and `T` must be positive

    Examples
    --------
    >>> arrhenius_rate(10**7, 10**3, 100, b=0.5)
    30035490.889639616

    >>> arrhenius_rate(10, 100, 300)
    9.6070007471344585

    >>> arrhenius_rate(10, 100, 300, 0.5)
    166.39813402389049
    """
    try:
        assert isinstance(A, numbers.Number)
        assert isinstance(b, numbers.Number)
        assert isinstance(R, numbers.Number)
    except AssertionError:
        raise TypeError("Arguments A, b, and R must be numbers.")
    if not A > 0:
        raise ValueError("Argument `A` (%f) must be > 0." % A)
    if not np.all(T > 0):
        raise ValueError("Temperature must be positive.")

    E_arr = np.array(E)
    T_arr = np.array(T)

    if E_arr.shape and T_arr.shape:
        if E_arr.shape != T_arr.shape:
            raise TypeError(
                "Arguments E and T are arrays of different dimensions")

    rate_coeff = A * T ** b
    rate_coeff *= np.exp( -E_arr / R / T_arr)
    return rate_coeff


class ElementaryReaction():
    pass


class ReactionSystem():
    """

    """
    def __init__(elementary_reactions):
        pass

    def __repr__():
        pass

    def calculate_reaction_coefficients(temperature):
        pass

    def calculate_progress_rate():
        pass

    def calculate_reaction_rate():
        pass

    def get_species():
        pass


class XMLReader():
    """
    # XMLReader will eventually create a ReactionSystem

    """
    def __init__(self, xml_file):
        self.xml_file = xml_file
        
    

        
    def build_reaction_system():
        # return ReactionSystem()
        pass
    

    def __repr__():
        pass



if __name__ == "__main__":
    reader = XMLReader(xml_file)
    reaction_system = reader.build_reaction_system()
    print ("rates: ", reaction_system.calculate_reaction_rate())

