import sys
import numpy as np
def progress_rate(v, x, k):
    """Returns the progress rate for a reaction.
    
    INPUTS
    =======
    v: array or list
       Stoichiometric coefficients
    x: array or list 
       Concentration of species
    k: float, must be positive
       Reaction rate coefficient
    
    RETURNS
    ========
    w: float
       Progress rate for a reaction
    
    EXAMPLES
    =========
    >>> progress_rate([2.0, 1.0, 1.0], [1.0, 2.0, 3.0], 10)
    20.0
    """
    if k <= 0:
        raise ValueError("Reaction rate coefficient must be positive.")
    if len(v)!=len(x) or len(x)==0 or len(v)==0:
        raise  ValueError("Invalid length of input parameters.")
    w = k
    for i in range(len(x)-1):
        w *= (x[i]**v[i])
    return w

def progress_rate_multiple(v1, x, k):
    """Returns the progress rate for a system of reactions.
    
    INPUTS
    =======
    v1: matrix
       Stoichiometric coefficients of reactants
    x: array
       Concentration of species
    k: array, all elements must be positive
       Reaction rate coefficients
    
    RETURNS
    ========
    w: array
       Progress rate for a system of reactions
    
    EXAMPLES
    =========
    >>> progress_rate_multiple([[1.0, 2.0], [2.0, 0.0], [0.0, 2.0]], [1.0, 2.0, 1.0], 10)
    [40.0, 10.0]
    """
    if np.any(k) <= 0:
        raise ValueError("Reaction rate coefficient must be positive.")
    if type(k) is not list:
        k = [k]*len(v1[0])
    if len(v1)!=len(x) or len(v1)==0 or len(x)==0 or len(v1[0])!=len(k):
        raise  ValueError("Invalid input parameters.")
    w = [0.0]*len(v1[0])
    for i in range(len(v1[0])):
        wi = k[i]
        for j in range(len(x)):
            wi *= (x[j] ** v1[j][i])
        w[i] = wi
    return w

def reaction_rate(v1, v2, x, k, progress_rate_func):
    """Returns the reaction rate of a system of irreversible reactions.
    
    INPUTS
    =======
    v1: matrix, same size as v2
       Stoichiometric coefficients of reactants
    v2: matrix, same size as v1
       Stoichiometric coefficients of products
    x: array
       Concentration of species
    k: array, all elements must be positive
       Reaction rate coefficient
    rogress_rate_func: function
       Function to calculate progress rates
    
    RETURNS
    ========
    f: matrix
       Reaction rate of a system of irreversible reactions
    
    EXAMPLES
    =========
    >>> reaction_rate([[1.0, 0.0], [2.0, 0.0], [0.0, 2.0]], [[0.0, 1.0], [0.0, 2.0], [1.0, 0.0]], [1.0, 2.0, 1.0], 10, progress_rate_multiple)
    matrix([[-30.],
            [-60.],
            [ 20.]])
    """
    if len(v2)==0 or len(v1)!=len(v2) or len(v1[0])!=len(v2[0]):
        raise  ValueError("Invalid input parameters.")
    w = progress_rate_func(v1, x, k)
    v1_1 = np.asmatrix(v1)
    v2_1 = np.asmatrix(v2)
    v = v2_1 - v1_1
    w1 = np.transpose(np.asmatrix(w))
    f = np.dot(v, w1)
    return f
