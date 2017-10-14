
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
    
