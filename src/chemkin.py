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

