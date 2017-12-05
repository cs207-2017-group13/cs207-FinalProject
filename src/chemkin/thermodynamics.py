"""
thermodynamics

"""
import os
import collections
import functools
import sqlite3
import numpy as np


# memoization decorator; from
# https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
class memoized():
    """Decorator. Caches a function's return value each time it is
    called. If called later with the same arguments, the cached value
    is returned (not reevaluated).
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)


# rename?
class Thermochem():
    """Methods for calculating the backward reaction rate.

    Take `Rxnset` instance and calculate backward reaction rate using
    the temperature passed into the class and the corresponding nasa
    polynomial coefficients.

    Parameters
    ==========
    rxnset : Rxnset
        Containing nasa polynomial coefficients for all species, the
        difference of product coefficients matrix and reactant
        coefficients matrix
    p0 : float
        Pressure of the reactor, default 1e5 (Pa)
    R : float
        Ideal gas constant, default 8.3144598 (J/mol K)

    Methods
    =======
    Cp_over_R(T)
        Returns specific heat of each species given by the NASA
        polynomials.
    H_over_RT(T)
        Returns the enthalpy of each species given by the NASA
        polynomials.
    S_over_R(T)
        Returns the entropy of each species given by the NASA
        polynomials.
    backward_coeffs(kf, T)
        Returns the backward reaction rate coefficient for reach
        reaction.
    """
    def __init__(self, rxnset, p0=1e5, R=8.3144598):
        assert isinstance(rxnset, Rxnset)
        self.rxnset = rxnset
        self.p0 = p0   # Pa
        self.R = R     # J / mol / K
        self.gamma = np.sum(self.rxnset.nuij, axis=0)

    def Cp_over_R(self, T):
        """Returns specific heat of each species given by the NASA
        polynomials.

        Parameters
        ==========
        T : float
            Reaction temperature

        Returns
        =======
        Cp_R : np.array
            specific heat of each species

        Examples
        ========
        >>> import chemkin.chemkin as chemkin
        >>> reader = chemkin.XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
        >>> rxnset = Rxnset(reaction_system[0].species, nu)
        >>> thermochem = Thermochem(rxnset)
        >>> thermochem.Cp_over_R(300)
        array([ 2.5       ,  2.63400196,  3.59349336,  3.46995186,  4.04072434,
                3.53457253,  4.2011126 ,  5.10656488])
        """
        a = self.rxnset.get_nasa_coefficients(T)

        Cp_R = (a[:, 0] + a[:, 1] * T + a[:, 2] * T**2.0
                + a[:, 3] * T**3.0 + a[:, 4] * T**4.0)
        return Cp_R

    def H_over_RT(self, T):
        """Returns the enthalpy of each species given by the NASA
        polynomials.

        Parameters
        ==========
        T : float
            Temperature of elementary reactions

        Returns
        =======
        H_RT : np.array
            enthalpy of each species

        Examples
        ========
        >>> import chemkin.chemkin as chemkin
        >>> reader = chemkin.XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
        >>> rxnset = Rxnset(reaction_system[0].species, nu)
        >>> thermochem = Thermochem(rxnset)
        >>> thermochem.H_over_RT(300)
        array([  8.74121997e+01,   9.99119166e+01,   1.57966367e+01,
                 2.13926454e-02,  -9.69244747e+01,   2.17928611e-02,
                 5.05805989e+00,  -5.44435841e+01])
        """
        a = self.rxnset.get_nasa_coefficients(T)

        H_RT = (a[:, 0] + a[:, 1] * T / 2.0 + a[:, 2] * T**2.0 / 3.0
                + a[:, 3] * T**3.0 / 4.0 + a[:, 4] * T**4.0 / 5.0
                + a[:, 5] / T)
        return H_RT

    def S_over_R(self, T):
        """Returns the entropy of each species given by the NASA
        polynomials.

        Parameters
        ==========
        T : float
            Temperature of elementary reactions

        Returns
        =======
        S_R : list of floats
            entropy of each species

        Examples
        ========
        >>> import chemkin.chemkin as chemkin
        >>> reader = chemkin.XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
        >>> rxnset = Rxnset(reaction_system[0].species, nu)
        >>> thermochem = Thermochem(rxnset)
        >>> thermochem.S_over_R(300)
        array([ 13.81277333,  19.38730759,  22.12090629,  15.73868267,
                22.73578462,  24.69552926,  27.58089871,  28.23850226])
        """
        a = self.rxnset.get_nasa_coefficients(T)

        S_R = (a[:, 0] * np.log(T) + a[:, 1] * T + a[:, 2] * T**2.0 / 2.0
               + a[:, 3] * T**3.0 / 3.0 + a[:, 4] * T**4.0 / 4.0 + a[:, 6])
        return S_R

    def backward_coeffs(self, kf, T):
        """Calculates backward rate coefficient.

        First calculates equilibrium constant from thermodynamic
        quantities. Then calculates backward rate coefficient from
        equilibrium constant and forward reaction rate coefficient.

        For irreversible reactions, the backward rate coefficients
        will be approximately 0.

        Parameters
        ==========
        kf : np.array
            Forward reaction rate coefficients
        T : float
            Temperature of elementary reactions

        Returns
        =======
        kf / kb : np.array
            Backward reaction rate coefficient for each reaction

        Examples
        ========
        >>> import chemkin.chemkin as chemkin
        >>> reader = chemkin.XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> k = reaction_system[0].get_rate_coefficients(300)
        >>> nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
        >>> rxnset = Rxnset(reaction_system[0].species, nu)
        >>> thermochem = Thermochem(rxnset)
        >>> thermochem.backward_coeffs(k,300)
        array([  4.27113939e+22,   1.97842662e+11,   1.11232410e+01,
                 4.43114607e+20,   2.25385449e-27,   1.32780748e-14,
                 7.24055199e-26,   2.46671212e-37,   8.33522066e-39,
                 4.32365740e-01,   4.48562565e-01])
        """
        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.rxnset.nuij.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.rxnset.nuij.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = (self.p0 / self.R / T)**self.gamma

        # Ke
        kb = fact * np.exp(delta_G_over_RT)

        return kf / kb


class Rxnset():
    """Read and store certain NASA polynomial coefficients from the
    SQL database containing coefficients for all species.

    Parameters
    ==========
    species : list
        Names of all species in the reaction system.
    nuij : np.ndarray
        The difference between product and reactant coefficient
        matrices. Dimensions n_species x n_reactions.

    Methods
    =======
    get_nasa_coefficients(T)
    read_nasa_coefficients()
    """
    def __init__(self, species, nuij):
        self.species = species
        self.nuij = nuij
        self.coefficients = self.read_nasa_coefficients()

    @memoized
    def get_nasa_coefficients(self, T):
        """Returns the corresponding NASA polynomial coefficients for
        all species at the given temperature

        Parameters
        ==========
        T : float
            Temperature of elementary reactions

        Returns
        =======
        coeffs : np.array (2D)
            NASA polynomial coefficients for all species at the given
            temperature

        Examples
        ========
        >>> import chemkin.chemkin as chemkin
        >>> reader = chemkin.XMLReader("tests/rxns_reversible.xml")
        >>> reaction_system = reader.get_reaction_systems()
        >>> nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
        >>> rxnset = Rxnset(reaction_system[0].species, nu)
        >>> rxnset.get_nasa_coefficients(300)[1]
        array([  3.16826710e+00,  -3.27931884e-03,   6.64306396e-06,
                -6.12806624e-09,   2.11265971e-12,   2.91222592e+04,
                 2.05193346e+00])
        """
        coeffs = np.zeros([len(self.species), 7])
        for i, (species, data) in enumerate(self.coefficients.items()):
            if T > data['low_T'] and T < data['mid_T']:
                coeffs[i, :] = data['low_coeffs']
            elif T >= data['mid_T'] and T < data['high_T']:
                coeffs[i, :] = data['high_coeffs']
            else:
                raise ValueError("Temperature not in range provided by NASA "
                                 "polynomial coefficients.")
        return coeffs

    def read_nasa_coefficients(self):
        """Return NASA polynomial coefficients for all species.

        Returns
        =======
        coeffs : collections.OrderedDictionary
            Dictionary of dictionaries of species, each dictionary of
            species contains low, mid, high temperature, and NASA
            polynomial coefficients for two temprature ranges for all
            species

        Examples
        ========
        """
        db_location = os.path.dirname(
            __file__) + '/NASA_polynomial_coefficients.sqlite'
        db = sqlite3.connect(db_location)
        cursor = db.cursor()
        coeffs = collections.OrderedDict()
        for species in self.species:
            query = '''SELECT low_temp, mid_temp, high_temp,
            low_coeffs, high_coeffs
            FROM coefficients WHERE species="{}"'''.format(species)
            low_T, mid_T, high_T, low_coeffs, high_coeffs = cursor.execute(
                query).fetchone()

            coeffs[species] = {}
            coeffs[species]['low_coeffs'] = np.array(
                [float(c) for c in low_coeffs.split()])
            coeffs[species]['high_coeffs'] = np.array(
                [float(c) for c in high_coeffs.split()])
            coeffs[species]['low_T'] = low_T
            coeffs[species]['mid_T'] = mid_T
            coeffs[species]['high_T'] = high_T
        db.close()
        return coeffs


if __name__ == "__main__":
    import doctest
    doctest.testmod()
