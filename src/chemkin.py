import sys
import xml.etree.ElementTree as ET
import numpy as np
import numbers

IDEAL_GAS_CONSTANT = 8.314


class ReactionSystem():
    """

    """
    def __init__(self, elementary_reactions, species):
        self.elementary_reactions = elementary_reactions
        self.species = species
        self.reactant_coefficients = self.build_reactant_coefficient_matrix()
        self.product_coefficients = self.build_product_coefficient_matrix()


    def __repr__():
        print("Species: ", self.species, "\n",
              "Stoichiometric coefficients of reactants: ", self.reactant_coefficients, "\n",
              "Stoichiometric coefficients of products: ", self.product_coefficients, "\n")


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
        """
        k = self.calculate_reaction_coefficients(temperature)
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
        """
        if self.reactant_coefficients.shape != self.product_coefficients.shape:
            raise  ValueError("Invalid input parameters.")
        nu = self.product_coefficients - self.reactant_coefficients
        rj = self.calculate_progress_rate(concs, temperature)
        return np.dot(nu, rj)

    def calculate_reaction_coefficients(temperature):
        """Calculate reaction rate coefficients
        
        INPUTS:
        =======
        temperature: numpy array of floats
                     temperatures of the elementary reactions

        RETURNS:
        ========
        f: numpy array of floats
           reaction rate ooefficients
        """
        coefficients = [er.get_k(temperature) for er
                        in self.elementary_reactions]
        return coefficients

    def build_reactant_coefficient_matrix():
        """Build a reactant coefficients matrix for the reaction system

        RETURNS:
        ========
        f: numpy array of floats
           reactant coefficients matrix
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_react = reaction.get_reactant_coefficients()
            for j,species in self.species:
                mat[j,i] = dict_react.get(species, 0)

        return mat

    def build_product_coefficient_matrix():
        """Build a product coefficients matrix for the reaction system

        RETURNS:
        ========
        f: numpy array of floats
           product coefficients matrix
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_prod = reaction.get_product_coefficients()
            for j,species in self.species:
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
        properties = {}
        prop = reaction_elt.attrib

        return attributes

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
        """
        """
        species = self._get_species()
        reaction_systems = []
        for reaction_data in self.root.findAll('reactionData'):
            elementary_reactions = []
            for reaction in reaction_data.findAll('reaction'):
                reaction_properties = self._parse_reaction(reaction)
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


