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

