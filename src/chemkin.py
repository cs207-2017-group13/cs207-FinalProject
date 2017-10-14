
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
    def __init__(xml_file):
        pass

    def build_reaction_system():
        # return ReactionSystem()
        pass
    

    def __repr__():
        pass



if __name__ == "__main__":
    reader = XMLReader(xml_file)
    reaction_system = reader.build_reaction_system()
    print ("rates: ", reaction_system.calculate_reaction_rate())
