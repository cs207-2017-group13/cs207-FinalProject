

import sys
import numbers

import xml.etree.ElementTree as ET
import numpy as np




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
      
    def calculate_rate_coefficient(self, T):
        if self.rate_type:
            if self.rate_type.lower() == 'constant':
                self.k = self.rate_params['k']
                rate_coeff = self._constant_rate()
            else:
                self.A = self.rate_params['A']
                self.b = self.rate_params['b']
                self.E = self.rate_params['E']
                rate_coeff = self._k_arrhenius(T)
        else:
            raise ValueError('No value for `type of rate passed`. Pass a value to get the reaction coeff')
        return rate_coeff
        

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


    def _k_arrhenius(self, temperature, R=8.314):
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

        if temperature < 0.0:
            raise ValueError("`T` in Kelvin scale  must be positive." )

        k_arrhenius = (self.A * temperature ** self.b) * np.exp(-self.E / (R * temperature))
        return k_arrhenius


class ReactionSystem():
    """Class for a system of reactions

    Take a list of dictionaries from the ElementaryReaction class. Build stoichiometric coefficient matrices for the 
    reactants and products, and calculate the corresponding progress rates and reaction rates.
    """
    def __init__(self, elementary_reactions, species):
        """Constructor

        INPUTS:
        =======
        elementary_reactions: a list of dictionaries 
                              each element in the list is a dictionary of the information for an elementary reaction
        species:              a list
                              species in the reaction system
        """
        self.elementary_reactions = elementary_reactions
        self.species = species
        self.reactant_coefficients = self.build_reactant_coefficient_matrix()
        self.product_coefficients = self.build_product_coefficient_matrix()


    def __repr__(self):
        """Returns a string containing basic information for the reaction system

        RETURNS:
        ========
        str: string
             Containing information on species and stoichiometric coefficients
        """
        info = "Species: {} \nStoichiometric coefficients of reactants: {} \nStoichiometric coefficients of products: {}".format(self.species, 
            self.reactant_coefficients, self.product_coefficients)
        return info


    def __len__(self):
        """Returns the number of species in the reaction system

        RETURNS:
        ========
        species_len: int
                     the number of species in the reaction system

        EXAMPLES:
        =========
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> len(reaction_system[0])
        5
        """
        return len(self.species)


    def calculate_progress_rate(self, concs, temperature):
        """Returns the progress rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        concs:    numpy array of floats 
                  concentration of species
        temperature: numpy array of floats
                     temperatures of the elementary reactions

        RETURNS:
        ========
        progress: numpy array of floats
                  size: num_reactions
                  progress rate of each reaction

        EXAMPLES:
        =========
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].calculate_progress_rate(concs, 300)
        [0.00024002941214766843, 39.005602653448953]
        """
        k = self.get_rate_coefficients(temperature)
        if type(k) is float:
            k = [k]*len(self.reactant_coefficients[0])
        if len(self.reactant_coefficients) != len(concs) or len(self.reactant_coefficients)==0 or len(self.reactant_coefficients[0])!=len(k):
            raise  ValueError("Invalid input parameters.")
        progress = k # Initialize progress rates with reaction rate coefficients
        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))
            for idx, xi in enumerate(concs):
                nu_ij = self.reactant_coefficients[idx,jdx]
                if xi  < 0.0:
                    raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij))
                
                progress[jdx] *= xi**nu_ij
        return progress   


    def calculate_reaction_rate(self, concs, temperature):
        """Returns the reaction rate of a system of irreversible, elementary reactions
        
        INPUTS:
        =======
        concs:    numpy array of floats 
                  concentration of species
        temperature: numpy array of floats
                     temperatures of the elementary reactions

        RETURNS:
        ========
        f: numpy array of floats
           size: num_species
           reaction rate of each specie

        EXAMPLES:
        =========
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].calculate_reaction_rate(concs, 300)
        array([  3.90053626e+01,  -3.90053626e+01,   3.90058427e+01,
                -3.90056027e+01,  -2.40029412e-04])
        """
        if self.reactant_coefficients.shape != self.product_coefficients.shape:
            raise  ValueError("Invalid input parameters.")
        nu = self.product_coefficients - self.reactant_coefficients
        rj = self.calculate_progress_rate(concs, temperature)
        return np.dot(nu, rj)


    def get_rate_coefficients(self, temperature):
        """Calculate reaction rate coefficients
        
        INPUTS:
        =======
        temperature: numpy array of floats
                     temperatures of the elementary reactions

        RETURNS:
        ========
        f: numpy array of floats
           reaction rate ooefficients

        EXAMPLES:
        =========
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].get_rate_coefficients(300)
        [0.00024002941214766843, 6.5009337755748255]
        """
        coefficients = [er.calculate_rate_coefficient(temperature) for er
                        in self.elementary_reactions]
        return coefficients


    def build_reactant_coefficient_matrix(self):
        """Build a reactant coefficients matrix for the reaction system

        RETURNS:
        ========
        f: numpy array of floats
           reactant coefficients matrix

        EXAMPLES:
        =========
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].build_reactant_coefficient_matrix()
        array([[ 1.,  0.],
               [ 0.,  1.],
               [ 0.,  0.],
               [ 0.,  1.],
               [ 1.,  0.]])
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_react = reaction.get_reactants()
            for j,species in enumerate(self.species):
                mat[j,i] = dict_react.get(species, 0)
        return mat


    def build_product_coefficient_matrix(self):
        """Build a product coefficients matrix for the reaction system

        RETURNS:
        ========
        f: numpy array of floats
           product coefficients matrix

        EXAMPLES:
        =========
        >>> reader = XMLReader("rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].build_product_coefficient_matrix()
        array([[ 0.,  1.],
               [ 1.,  0.],
               [ 1.,  1.],
               [ 0.,  0.],
               [ 0.,  0.]])
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_prod = reaction.get_products()
            for j,species in enumerate(self.species):
                mat[j,i] = dict_prod.get(species, 0)
        return mat

            


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
        for reaction_data in self.root.findall('reactionData'):
            elementary_reactions = []
            for reaction in reaction_data.findall('reaction'):
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



# if __name__ == "__main__":
    # reader = XMLReader(xml_file)
    # reaction_system = reader.build_reaction_system()
    # print ("rates: ", reaction_system.calculate_reaction_rate(T=234))
    # print ("rates: ", reaction_system.calculate_reaction_rate(T=400))


