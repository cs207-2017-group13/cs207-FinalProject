import xml.etree.ElementTree as ET
import numpy as np
import numbers

class ElementaryReaction():
    
    def __init__(self,reaction_properties):
        self.reaction_properties = reaction_properties
        self.rate_type =  self.reaction_properties['rate_type']
        self.rate_params = self.reaction_properties['rate_params']
        self.reactants = self.reaction_properties['reactants']
        self.products = self.reaction_properties['products']
        
    def get_reactants(self):
        return self.reactants
    
    def get_products(self):
        return self.products
      
    def calculate_rate_coefficients(self, T):
        self. T = T
        if self.reaction_type:
            if self.reaction_type.lower() == 'constant':
                self.k = self.reaction_params['k']
                reaction_coeff = self._constant_rate(self)
            else:
                self.A = self.reaction_params['A']
                self.b = self.reaction_params['b']
                self.E = self.reaction_params['E']
                reaction_coeff = self._k_arrhenius(self)              
        else:
            raise ValueError('No value for `type of rate passed`. Pass a value to get the reaction coeff')
        

    def _constant_rate(self, k = 1.0):
    """Return a constant reaction rate coefficient.

    In zeroth-order reactions, k = constant.

    INPUTS:
    =======
    k: float, default value = 1.0
       Constant reaction rate coefficient

    RETURNS:
    ========
    k: float
       Constant reaction rate coefficient

    Notes
    -----
    Although it would be sensible if input is numeric, no exceptions
    will be raised if this is not the case.

    EXAMPLES:
    =========
    >>> k_const(5.0)
    5.0
    """
    if self.k < 0:
        raise ValueError("Negative reaction rate cannot be negative. Check the value. ")
        
    return self.k


    def _k_arrhenius(self):
    """Return a reaction rate coefficient according to the Arrhenius equation.
    
    The Arrhenius equation relates the rate constant, k, of a chemical
    reaction to parameters A (pre-exponential factor), E (activation
    energy), T (absolute temperature), and b (exponential indicating
    temperature-dependence of pre-exponential factor)::

        k = A T^b exp[ -E / (RT) ]

    When ``b = 0``, the above formula corresponds to the Arrhenius equation.  
    A nonzero value for b gives the modified Arrhenius equation.
    
    INPUTS:
    =======
    A: float
       Arrhenius prefactor, Must be positive
    b: float,
       Modified Arrhenius parameter, default value = 0.0
    E: float
       Activation energy
    T: float
       Temperature, Must be positive
    R: float, default value = 8.314
       Ideal gas constant
    
    RETURNS:
    ========
    k: float
       Modified Arrhenius reaction rate coefficient
    
    EXAMPLES:
    =========
    >>> k_mod_arr(2.0, -0.5, 3.0, 100.0)
    0.19927962618542916
    """
    if self.A < 0.0:
        raise ValueError("`A` must be positive.")

    if self.T < 0.0:
        raise ValueError("`T` in Kelvin scale  must be positive." )
   
    k_arrhenius = (self.A * (np.power(self.T,self.b))) * np.exp(-self.E / (R * self.T))
    return k_arrhenius





class ReactionSystem():
    """

    """
    def __init__(elementary_reactions):
        pass

    def __repr__():
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
