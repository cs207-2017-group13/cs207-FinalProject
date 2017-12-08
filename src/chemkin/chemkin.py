import xml.etree.ElementTree as ET
import numpy as np

import chemkin.thermodynamics as thermodynamics
import chemkin.simulator as simulator
# import thermodynamics
# import simulator


class ElementaryReaction():
    """Class representing a single elementary reaction.

    Instantiate using a dictionary of reaction properties (can be
    generated by `XMLReader`). This class calculates the
    temperature-dependent rate coefficient for each elementary
    reaction.

    Instances of this class are contained by instances of
    `ReactionSystem`. Methods in this class serve to facilitate
    calculations by `ReactionSystem`.

    Parameters
    ----------
    reaction_properties : dict
        A dictonary of properties for the reaction, as parsed by
        `XMLReader`.

    Methods
    -------
    get_info()
        Obtain description of elementary reaction.
    get_reactants()
        Return reactants and stochiometric coefficients.
    get_products()
        Return products and stochiometric coefficients.
    calculate_reaction_order()
        Determine the order of the reaction.
    calculate_rate_coefficient(T)
        Calculate the rate coefficient according to the reaction type.

    Examples
    --------
    >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
    ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
    ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
    ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
    ...                     'reversible': 'no', 'type' : ' Elementary'})
    """
    def __init__(self, reaction_properties):
        assert isinstance(reaction_properties, dict)
        try:
            self.reaction_properties = reaction_properties
            self.rate_type = self.reaction_properties['rate_type']
            self.rate_params = self.reaction_properties['rate_params']
            self.reactants = self.reaction_properties['reactants']
            self.products = self.reaction_properties['products']
            self.reversible = True if self.reaction_properties[
                'reversible'] == 'yes' else False
        except KeyError as err:
            print("Key {} is not present in `reaction_properties`.".format(
                str(err)))

    def __repr__(self):
        """Returns a string containing class name and number of
        species in the elementary reaction.

        Returns
        -------
        info : str
        """
        return "<%s: %d species>" % (
            self.__class__.__name__,
            len(self.reactants) + len(self.products))

    def get_info(self):
        """Returns a string containing basic information about the
        Elementary reaction.

        Returns
        -------
        info : str
        """
        info = "Reactants: {} \nProducts: {} \nRate Params: {} "
        "\nRate Type: {} \nReversible: {}".format(
            self.reactants, self.products, self.rate_params,
            self.rate_type, self.reversible)
        return info

    def get_reactants(self):
        """Returns a dictionary with reactants as keys and the
        corresponding stoichiometric coefficients as values.

        Returns
        -------
        reactants : dict

        Examples
        --------
        >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
        ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
        ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
        ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
        ...                     'reversible': 'no', 'type' : ' Elementary'})
        >>> elementary_reaction.get_reactants()
        {'H': '1', 'O2': '1'}
        """
        return self.reactants

    def get_products(self):
        """Returns a dictionary with products as keys and the
        corresponding stoichiometric coefficients as values.

        Returns
        -------
        products : dict

        Examples
        --------
        >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
        ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
        ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
        ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
        ...                     'reversible': 'no', 'type' : ' Elementary'})
        >>> elementary_reaction.get_products()
        {'O': '1', 'OH': '1'}
        """
        return self.products

    def calculate_reaction_order(self):
        """Determine reaction order from number of reactants.

        Returns
        -------
        length : tuple
                If reaction is reversible, returns 2-tuple, for forward and
                backward reaction order. Otherwise returns a 1-tuple, for the
                forward reaction order.
        """
        if self.reversible:
            return (len(self.get_reactants()), len(self.get_products()))
        else:
            return (len(self.get_reactants()),)


    def calculate_rate_coefficient(self, T=None):
        """Calculates and returns a rate coefficient based on the type
        of rate coefficient required.

        Parameters
        ----------
        T : float
            Temperature in Kelvin scale. Must be positive. Optional
            (irrelevant) for constant rate reactions.

        Returns
        -------
        rate_coeff : float
            Rate coefficient based on the type of elementary reaction.

        Examples
        --------
        >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
        ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
        ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
        ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
        ...                     'reversible': 'no', 'type' : ' Elementary'})
        >>> elementary_reaction.calculate_rate_coefficient(1000)
        5.2102032610668552
        """
        if self.rate_type:
            if self.rate_type.lower() == 'constant':
                if 'k' in self.rate_params:
                    return self._constant_rate(self.rate_params['k'])
                else:
                    return self._constant_rate()
            elif 'Arrhenius' in self.rate_type:
                assert T is not None
                A = self.rate_params['A']
                E = self.rate_params['E']
                b = self.rate_params.get('b', 0.)
                return self._k_arrhenius(A, E, T, b)
            else:
                raise NotImplementedError(
                    'Rate type other than Constant, Arrhenius, and Modified '
                    'Arrhenius is not implemented.')
        else:
            raise ValueError('No value for `rate_type` passed. '
                             'Pass a value to get the reaction coeff.')

    def _constant_rate(self, k=1.0):
        """Return a constant reaction rate coefficient.

        In zeroth-order reactions, k = constant.

        Parameters
        ----------
        k : float, optional
            Constant reaction rate coefficient

        Returns
        -------
        k : float
            Constant reaction rate coefficient

        Notes
        -----
        Although it would be sensible if input is numeric, no exceptions
        will be raised if this is not the case.

        Examples
        --------
        >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
        ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
        ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
        ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
        ...                     'reversible': 'no', 'type' : ' Elementary'})
        >>> elementary_reaction._constant_rate(5.0)
        5.0
        """
        if k < 0:
            raise ValueError("Constant reaction rate cannot be negative. "
                             "Check the value. ")
        return k

    def _k_arrhenius(self, A, E, T, b=0.0, R=8.314):
        """Calculate a reaction rate coefficient according to the
        Arrhenius equation.

        The Arrhenius equation relates the rate constant, k, of a
        chemical reaction to parameters A (pre-exponential factor), E
        (activation energy), T (absolute temperature), and b
        (exponential indicating temperature-dependence of
        pre-exponential factor)::

            k = A T^b exp[ -E / (RT) ]

        When ``b = 0``, the above formula corresponds to the Arrhenius
        equation. A nonzero value for b gives the modified Arrhenius
        equation.

        Parameters
        ----------
        A : float
            Arrhenius prefactor, Must be positive
        E : float
            Activation energy
        T : float
            Temperature, Must be positive
        b : float, optional
            Modified Arrhenius parameter, default value=0.0,
            corresponding to regular Arrhenius.
        R : float, optional
            Ideal gas constant

        Returns
        -------
        k : float
            Arrhenius reaction rate coefficient

        Examples
        --------
        >>> elementary_reaction = ElementaryReaction({'equation' : 'H + O2  [=] OH + O' ,
        ...                     'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
        ...                     'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
        ...                     'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
        ...                     'reversible': 'no', 'type' : ' Elementary'})
        >>> elementary_reaction._k_arrhenius(3520000.0, 71400.0, 1000.0, -0.7)
        5.2102032610668552
        """
        if A < 0.0:
            raise ValueError("Activation Energy `A` must be positive.")

        if T < 0.0:
            raise ValueError(
                "Temperature `T` (in Kelvin scale) cannot be negative.")

        k_arrhenius = (A * T ** b) * np.exp(-E / (R * T))
        return k_arrhenius


class ReactionSystem():
    """Class representing a system of reactions.

    Takes a list of `ElementaryReaction` instances and a list of
    chemical species names. Builds stoichiometric coefficient matrices
    for the reactants and products and calculates the corresponding
    progress rates and reaction rates.

    Parameters
    ----------
    elementary_reactions : list
        A list of `ElementaryReaction` instances that compose the
        system of reactions.
    species : list
        A list of strings identifying species in reaction system.

    Methods
    -------
    calculate_progress_rate(concs, temperature)
    calculate_reaction_rate(concs, temperature)
    get_rate_coefficients()
        Return rate coefficients (a.k.a. rate constants) for reactions
        in the forward direction.
    get_backward_rate_coefficients()
        Return rate coefficients for reactions in the reverse
        direction. (0 for irreversible reactions)
    build_reactant_coefficient_matrix()
    build_product_coefficient_matrix()
    check_reversible()
    """
    def __init__(self, elementary_reactions, species):
        self.elementary_reactions = elementary_reactions
        self.species = species
        self.reactant_coefficients = self.build_reactant_coefficient_matrix()
        self.product_coefficients = self.build_product_coefficient_matrix()
        self.reversible = self.check_reversible()
        if any(self.reversible):
            nu = self.product_coefficients - self.reactant_coefficients
            rxnset = thermodynamics.Rxnset(self.species, nu)
            self.thermochem = thermodynamics.Thermochem(rxnset)
        else:
            self.thermochem = None

    def __repr__(self):
        return "<%s: %d reactions>" % (self.__class__.__name__,
                                       len(self.elementary_reactions))

    def get_info(self):
        """Returns a string containing basic information for the reaction system

        Returns
        -------
        info : str
              Containing information on species and stoichiometric coefficients
        """
        info = '''Species: {} \nStoichiometric coefficients of reactants: {} \n
            Stoichiometric coefficients of products: {}'''.format(
                self.species,
                self.reactant_coefficients, self.product_coefficients)
        return info

    def __len__(self):
        """Returns the number of species in the reaction system.

        Returns
        -------
        species_len : int
            The number of species in the reaction system

        Examples
        --------
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> len(reaction_system[0])
        5
        """
        return len(self.species)

    def calculate_progress_rate(self, concs, temperature):
        """Returns the progress rate of a system of elementary reactions.

        If the elementary reaction is reversible, subtract the
        backward progress rate from the forward progress rate to
        obtain the net progress rate.

        Parameters
        ----------
        concs : np.ndarray
            Concentration of species
        temperature : float
            Temperature of the elementary reactions

        Returns
        -------
        progress : np.ndarray
            Progress rate of each reaction. (size = number of reactions)

        Examples
        --------
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].calculate_progress_rate(concs, 300)
        [0.00024002941214766843, 39.005602653448953]
        """
        assert len(concs) == len(self.species)
        k = self.get_rate_coefficients(temperature)
        if type(k) is float:
            k = [k]*len(self.reactant_coefficients[0])

        # backward reaction rate coefficient
        kb = self.get_backward_rate_coefficients(temperature)

        # Initialize progress rates with reaction rate coefficients
        progress = k
        for jdx, rj in enumerate(progress):
            for idx, xi in enumerate(concs):
                nu_ij = self.reactant_coefficients[idx, jdx]
                # if xi < 0.0:
                #     raise ValueError(
                #         "x{0} = {1:18.16e}: Negative concentrations are "
                #         "prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError(
                        "nu_{0}{1} = {2}:  Negative stoichiometric "
                        "coefficients are prohibited!".format(idx, jdx, nu_ij))

                progress[jdx] *= xi**nu_ij

            # substract backward progress rate if reversible
            if self.reversible[jdx]:
                for idx, xi in enumerate(concs):
                    nuij2 = self.product_coefficients[idx, jdx]
                    if nuij2 < 0:
                        raise ValueError(
                            "nu_{0}{1} = {2}: Negative stoichiometric "
                            "coefficients are prohibited!".format(
                                idx, jdx, nuij2))
                    kb[jdx] *= xi**nuij2

                progress[jdx] -= kb[jdx]

        return progress

    def calculate_reaction_rate(self, concs, temperature):
        """Returns the reaction rate of a system of irreversible,
        elementary reactions.

        Parameters
        ----------
        concs : np.ndarray
            Concentration of species
        temperature : float
            Temperature of the elementary reactions

        Returns
        -------
        f : numpy array of floats
            Reaction rate (change in concentration) of each
            species. (size: number of species)

        Examples
        --------
        >>> concs = [1., 2., 1., 3., 1.]
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].calculate_reaction_rate(concs, 300)
        array([  3.90053626e+01,  -3.90053626e+01,   3.90058427e+01,
                -3.90056027e+01,  -2.40029412e-04])
        """
        assert len(concs) == len(self.species)
        nu = self.product_coefficients - self.reactant_coefficients
        rj = self.calculate_progress_rate(concs, temperature)
        return np.dot(nu, rj)

    def get_rate_coefficients(self, temperature):
        """Calculate reaction rate coefficients.

        These rate coefficients are for reactions in the forward
        direction.

        Parameters
        ----------
        temperature : array_like
            Temperatures

        Returns
        -------
        coefficients : np.ndarray
           reaction rate ooefficients

        Examples
        --------
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].get_rate_coefficients(300)
        [0.00024002941214766843, 6.5009337755748255]
        """
        coefficients = [er.calculate_rate_coefficient(temperature) for er
                        in self.elementary_reactions]
        return coefficients

    def get_backward_rate_coefficients(self, temperature):
        """Calculate the backward rate coefficients.

        Parameters
        ----------
        temperature : array_like
            Temperatures

        Returns
        -------
        kb : np.ndarray
            Backward rate cofficients

        Examples
        --------
        >>> reader = XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].get_backward_rate_coefficients(800)
        [  3.92875283e+16   1.71276093e+12   9.99213983e+08   1.39259837e+15
           3.26602470e-02   4.70481489e+02   1.11384696e-01   2.22076430e-05
           2.50630855e-07   1.23837445e+08   4.39701018e+07]
        """
        k = self.get_rate_coefficients(temperature)
        if self.thermochem:
            kb = self.thermochem.backward_coeffs(k, temperature)
        else:
            kb = np.zeros(len(self.elementary_reactions))
        return kb

    def build_reactant_coefficient_matrix(self):
        """Build a reactant coefficients matrix for the reaction system.

        Returns
        -------
        mat : np.ndarray
           reactant stoichiometric coefficients

        Examples
        --------
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_systems = reader.get_reaction_systems()
        >>> reaction_systems[0].build_reactant_coefficient_matrix()
        array([[ 1.,  0.],
               [ 0.,  1.],
               [ 0.,  0.],
               [ 0.,  1.],
               [ 1.,  0.]])
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_react = reaction.get_reactants()
            for j, species in enumerate(self.species):
                mat[j, i] = dict_react.get(species, 0)
        return mat

    def build_product_coefficient_matrix(self):
        """Build a product coefficients matrix for the reaction system.

        Returns
        -------
        mat : np.ndarray
           product stoichiometric coefficients

        Examples
        --------
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_systems = reader.get_reaction_systems()
        >>> reaction_systems[0].build_product_coefficient_matrix()
        array([[ 0.,  1.],
               [ 1.,  0.],
               [ 1.,  1.],
               [ 0.,  0.],
               [ 0.,  0.]])
        """
        mat = np.zeros([len(self.species), len(self.elementary_reactions)])
        for i, reaction in enumerate(self.elementary_reactions):
            dict_prod = reaction.get_products()
            for j, species in enumerate(self.species):
                mat[j, i] = dict_prod.get(species, 0)
        return mat

    def check_reversible(self):
        """Check if each elementary reaction is reversible

        Returns
        -------
        reversible : list
            a list of True and False indicating if each elementary
            reaction is reversible

        Examples
        --------
        >>> reader = XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].check_reversible()
        [False, False]
        """
        reversible = [i.reversible for i in self.elementary_reactions]
        return reversible

    def setup_reaction_simulator(self, simulation_type, abundances, temperature,
                                t_span, dt=0.01, system_volume=1e-15):
        """Deterministic simulation or stochastic simulation

        Parameters
        ----------
        simulation_type : string
                          Type of simulation, deterministic or stochastic
        abundances :      array_like
                          Abundances of all species
        temperature :     array_like
                          Temperatures
        t_span :          tuple of floats
                          Time span of the reactions users want to 
                          simulate
        dt :              float
                          Size of time steps users want to simulate
        system_volume : float
                        System volume

        Returns
        -------
        kb : np.ndarray
            Backward rate cofficients

        Examples
        --------
        >>> reader = XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> reaction_system[0].get_backward_rate_coefficients(800)
        array([  3.92875283e+16,   1.71276093e+12,   9.99213983e+08,
                 1.39259837e+15,   3.26602470e-02,   4.70481489e+02,
                 1.11384696e-01,   2.22076430e-05,   2.50630855e-07,
                 1.23837445e+08,   4.39701018e+07])
        """
        choices = ['stochastic', 'deterministic']
        if simulation_type not in choices:
            raise ValueError

        if simulation_type == 'deterministic':
            determine_sim = simulator.DeterministicSimulator(
                self, abundances, temperature, t_span, dt)
            return determine_sim
        else:
            stochastic_sim = simulator.StochasticSimulator(
                self, abundances, temperature, t_span, system_volume)
            return stochastic_sim


class XMLReader():
    """Parser for chemical reaction XML files. Uses `xml.etree`.

    Parameters
    ----------
    xml_file : str
         Path to XML file.

    Attributes
    ----------
    xml_file : str
         Path to XML file.
    root : xml.etree.Element
         Top-level element in parsed tree.

    Methods
    -------
    get_reaction_system()
         Parse all groups of reactions in the XML file.

    Examples
    --------
    >>> reader = XMLReader("tests/rxns.xml")
    >>> reaction_systems = reader.get_reaction_systems()
    """
    def __init__(self, xml_file):
        self.xml_file = xml_file
        try:
            xml_tree = ET.parse(xml_file)
        except ET.ParseError:
            raise SyntaxError("Failed to parse %s" % xml_file)
        self.root = xml_tree.getroot()

    def _get_species(self):
        """Return the species involved in the reaction.

        Returns
        -------
        species : list
            A list of strings identifying the species involved in the
            reaction.

        Raises
        ------
        LookupError
            "speciesArray" element with text missing.
        """
        try:
            phase_elt = self.root.find('phase')
            species_array_elt = phase_elt.find('speciesArray')
            species_text = species_array_elt.text
        except AttributeError:
            raise LookupError("Could not find any species in "
                              "root>phase>speciesArray")
        species = species_text.split()
        return species

    def parse_reaction(self, reaction_elt):
        """Collect individual reaction properties in a dictionary.

        The first "rateCoeff" entry found is converted into two items,
        `rate_type` and `rate_params`, the former containing the
        reaction type and the latter being a dictionary of rate
        parameters.

        Parameters
        ----------
        reaction_elt : xml.etree.Element
             A "reaction" entry within a "reactionData"
             entry. Corresponds to an elementary reaction.

        Returns
        -------
        properties : dict
             Dictionary of reaction parameters.
        """
        try:
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
        except AttributeError:
            raise LookupError("Could not properly parse reaction element "
                              "in xml file %s" % self.xml_file)
        return properties

    def get_reaction_systems(self):
        """Parse all groups of reactions in xml file.

        Returns
        -------
        reaction_systems : list
             A list of `ReactionSystem` instances corresponding to
             each "reactionData" entry in the XML file.

        Raises
        ------
        NotImplementedError
            Upon encountering reaction of not of type "Elementary".
        """
        species = self._get_species()
        reaction_systems = []
        for reaction_data in self.root.findall('reactionData'):
            elementary_reactions = []
            for reaction in reaction_data.findall('reaction'):
                reaction_properties = self.parse_reaction(reaction)
                if reaction_properties['type'] != "Elementary":
                    raise NotImplementedError("Current implementation only "
                                              "supports elementary reactions")
                elementary_reaction = ElementaryReaction(reaction_properties)
                elementary_reactions.append(elementary_reaction)
            reaction_system = ReactionSystem(elementary_reactions, species)
            reaction_systems.append(reaction_system)

        return reaction_systems

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.xml_file)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
