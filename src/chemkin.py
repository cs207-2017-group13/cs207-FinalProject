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
    def __init__(self, elementary_reactions, species):
        self.elementary_reactions = elementary_reactions
        self.species = species
        self.reactant_coefficients = []
        self.product_coefficients = []
        pass

    def __repr__():
        pass

    def get_reaction_coefficients(temperature):
        coefficients = [er.get_k(temperature) for er
                        in self.elementary_reactions]
        return coefficients

    def calculate_progress_rate():
        pass

    def calculate_reaction_rate():
        pass

    def get_species():
        pass

    def build_reactant_coefficient_matrix():
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i,reaction in enumerate(self.elementary_reactions):
            dict_ = reaction.get_reactant_coefficients()
            for j,species in self.species:
                mat[i,j] = dict_.get(species, 0)
            


class XMLReader():
    """
    # XMLReader will eventually create a ReactionSystem

    """
    def __init__(self, xml_file):
        self.xml_file = xml_file
        xml_tree = ET.parse(xml_file)
        self.root = xml_tree.getroot()

    def _parse_reaction(self, reaction_elt):
        """Collect individual reaction properties in a dictionary."""
        properties = reaction_elt.attrib.copy()
        properties["equation"] = reaction_elt.find("equation").text

        rate_coeff = reaction_elt.find("rateCoeff")
        rate_coeff_child = rate_coeff.getchildren()[0]
        properties["rate_type"] = rate_coeff_child.tag
        properties["rate_params"] = {}
        for child in rate_coeff_child.getchildren():
            properties["rate_params"][child.tag] = float(child.text)

        reactants = reaction_elt.find("reactants").text.split()
        properties["reactants"] = {}
        for reactant in reactants:
            species, coefficient = reactant.split(':')
            properties['reactants'][species] = coefficient

        products = reaction_elt.find("products").text.split()
        properties["products"] = {}
        for reactant in products:
            species, coefficient = reactant.split(':')
            properties['products'][species] = float(coefficient)

        return properties

    def _get_species(self):
        """Return the species involved in the reaction."""
        try:
            phase_elt = self.root.find('phase')
            species_array_elt = phase_elt.find('speciesArray')
            species_text = species_array_elt.text
        except AttributeError:
            raise LookupError("Element root>phase>speciesArray")
        species = species_text.split()
        return species

    def get_reaction_systems(self):
        """Parse all groups of reactions in xml file.
        """
        species = self._get_species()
        reaction_systems = []
        for reaction_data in self.root.findAll('reactionData'):
            elementary_reactions = []
            for reaction in reaction_data.findAll('reaction'):
                reaction_properties = self._parse_reaction(reaction)
                if reaction_properties['type'] != "Elementary":
                    raise NotImplementedError
                elementary_reaction = ElementaryReaction(reaction_properties)
                elementary_reactions.append(elementary_reaction)
            reaction_system = ReactionSystem(elementary_reactions, species)
            reaction_systems.append(reaction_system)

        return reaction_systems    

    def __repr__(self):
        return "XMLReader(%s)" % self.xml_file



if __name__ == "__main__":
    reader = XMLReader(xml_file)
    reaction_system = reader.build_reaction_system()
    print ("rates: ", reaction_system.calculate_reaction_rate(T=234))
    print ("rates: ", reaction_system.calculate_reaction_rate(T=400))

