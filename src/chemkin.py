
def progress_rate(X, V , k):
    """
    The function returns the progress rate of a one elementary reaction
    based on the formula ω = k*∏x^(ν′)
        
    INPUTS
    ======
    X : a list of floats, positive, required,
        concentration of species
    V : a list of floats, non negative, required,
        stoichiometric coefficients of the reactants
    k:  float, positive,required,
        reaction rate coefficient
      
  
        
    RETURNS
    =======
    progress_rate : float, unless k and Xi are non positive,
                    and Vi is negative, in which case it raises 
                    Value Error  
    NOTES
    =====
    PRE:
        - The inputs X,V and k are required
        - X,V are list of type float
        - k is float
        - Xi and k should be positive
        - Vi should be non negative


    POST:
        - raises value error when Vi is negative
        - raises value error when Xi and k are non positive
        - returns a float
    
    EXAMPLES
    ========
    >>> progress_rate([1.0,2.0,3.0],[2.0,1.0,0.0],10)
    20.0
    """
    import numpy as np
    
    product = 1
    
    if k <= 0:
        raise ValueError("The reaction rate cannot be 0 or negative.") # k should be positive
        
    for X_i, V_i in zip(X,V):
        
        if X_i <= 0:
            raise ValueError("The specie concentration cannot be 0 or negative.") # Xi should be positive
            
        if V_i < 0:
            raise ValueError("The reactant coefficient cannot be less than 0")   # The reactant coeff can't be less than zero
            
        product = product  * np.power(X_i,V_i)
        
    progress_rate = k * product
    
    return progress_rate  
   
    
def progress_rate_m(V_prime , V_double_prime, X, k_list):
    """
    The function returns the progress rate of a reaction in a system of multiple elementary reactions
        
    INPUTS
    ======
    X :              a list of floats, positive, required,
                     Concentration of species
    V_prime :        a list of floats, non negative, required,
                     Stoichiometric coefficients of reactants
    V_double_prime : a list of floats, non negative, required,
                     stoichiometric coefficients of products        
    k_list:          a list of float, positive, required,
                     reaction rate coefficient for each elementary reaction
      
  
        
    RETURNS
    =======
    progress_rates : a list of floats, unless k and Xi are is non positive, 
                    and Vi is negative, in which case it raises Value Error  
    NOTES
    =====
    PRE:
        - The inputs X,V_prime, V_double_prime and k are required
        - X,V_prime, V_double_prime are lists of type float
        - k is float
        - Xi and k should be positive
        - Vi should be non negative

    POST:
        - raises value error when Vi is negative
        - raises value error when Xi and k are non positive
        - returns a list of float
    
    EXAMPLES
    ========
    >>> progress_rate_m([[1.0,2.0,0.0], [2.0, 0.0,2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]],[1.0, 2.0, 1.0],10)
    [40.0, 10.0]
    """
    import numpy as np
    progress_rates = []
    
    for k, reaction_coeff in zip(k_list,V_prime):
        product = 1
        if k <= 0:
            raise ValueError("The reaction rate cannot be 0 or negative.") # k should be positive
            
        for X_i, V_i in zip(X,reaction_coeff):
            if X_i <= 0:
                raise ValueError("The specie concentration cannot be 0 or negative.") # Xi should be positive
            
            if V_i < 0:
                raise ValueError("The reactant coefficient can't be negative.") # Vi should not be negative
                
            product = product  * np.power(X_i,V_i)
        progress_rate = k * product       
        progress_rates.append(progress_rate) 
    return progress_rates
    
    
def reaction_rate_m(V_prime , V_double_prime, X, progress_rates):
    """
    The function returns the reaction rate of species in a system of multiple elementary reactions
        
    INPUTS
    ======
    X :              a list of floats, positive, required,
                     Concentration of species
    V_prime :        a list of floats, non negative, required,
                     Stoichiometric coefficients of reactants
    V_double_prime : a list of floats, non negative, required,
                     stoichiometric coefficients of products        
    progress_rates:  a list of floats, positive, required,
                     progress rate of each reaction
      
  
        
    RETURNS
    =======
    reaction_rate : list of float, unless Xi is non positive,
                    and Vi is negative, in which case it raises Value Error  
    NOTES
    =====
    PRE:
        - The inputs X,V_prime, V_double_prime and k are required
        - X,V_prime, V_double_prime are lists of type float
        - k is float
        - Xi and k should be positive
        - Vi should be non negative


    POST:
        - raises value error when Vi is negative
        - raises value error when Xi and k are non positive
        - returns a list of floats
    
    EXAMPLES
    ========
    >>> reaction_rate_m([[1.0,2.0,0.0], [2.0, 0.0,2.0]],[[0.0, 0.0, 2.0], [0.0, 1.0, 1.0]],[1.0, 2.0, 1.0],10)
    [-60.0, -70.0, 70.0]
    """
    import numpy as np    
    reaction_rate = 0;
    
    for X_i in X:
        if (X_i < 0):
            raise ValueError("The specie concentration can't be negative") 

    for V_i_row, V_double_i_row in zip(V_prime,V_double_prime):
        if any(V_i < 0 for V_i in V_i_row):
            raise ValueError("The stochiometric  reactant coefficient can't be negative")
        if any(V_double_i < 0 for V_double_i in V_double_i_row ):
            raise ValueError("The stochiometric  product coefficient can't be negative")

    V = np.subtract(V_double_prime,V_prime) 
    for row,progress_rate in zip(V,progress_rates):
        V_mult = np.multiply(row, progress_rate)
        reaction_rate = reaction_rate + V_mult
    return reaction_rate.tolist()
    
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
