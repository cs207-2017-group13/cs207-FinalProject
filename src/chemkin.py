import xml.etree.ElementTree as ET
import numpy as np
import numbers



class ElementaryReaction():
    pass


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

