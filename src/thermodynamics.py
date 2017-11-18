import collections
import sqlite3
import numpy as np


class Thermochem():
    """Methods for calculating the backward reaction rate.
    Take Rxnset class object and calculate backward reaction rate
    using the temperature passed into the clsss and the 
    corresponding nasa polynomial coefficients.

    Parameters
    ==========
    rxnset: class
            containing nasa polynomial coefficients for all species,
            the difference of product coefficients matrix and 
            reactant coefficients matrix
    p0:     float
            pressure of the reactor, default 1e5
    R:      float
            ideal gas constant, default 8.3144598

    Methods
    =======
    Cp_over_R(T)

    H_over_RT(T)

    S_over_R(T)

    backward_coeffs(kf, T)
    """
    def __init__(self, rxnset, p0=1e5, R=8.3144598):
        self.rxnset = rxnset
        self.p0 = p0 # Pa
        self.R = R # J / mol / K
        self.gamma = np.sum(self.rxnset.nuij, axis=0)

    def Cp_over_R(self, T):
        """Returns specific heat of each species given by
           the NASA polynomials.
        
        INPUTS:
        =======
        T:    float
              Temperature of elementary reactions

        RETURNS:
        ========
        Cp_R: list of floats
              specific heat of each species

        EXAMPLES:
        =========
        """
        a = self.rxnset.get_nasa_coefficients(T)

        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0 
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)

        return Cp_R

    def H_over_RT(self, T):
        """Returns the enthalpy of each species given by
           the NASA polynomials.
        
        INPUTS:
        =======
        T:    float
              Temperature of elementary reactions

        RETURNS:
        ========
        H_RT: list of floats
              enthalpy of each species

        EXAMPLES:
        =========
        """
        a = self.rxnset.get_nasa_coefficients(T)

        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)

        return H_RT
               

    def S_over_R(self, T):
        """Returns the entropy of each species given by
           he NASA polynomials.
        
        INPUTS:
        =======
        T:   float
             Temperature of elementary reactions

        RETURNS:
        ========
        S_R: list of floats
             entropy of each species

        EXAMPLES:
        =========
        """
        a = self.rxnset.get_nasa_coefficients(T)

        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    def backward_coeffs(self, kf, T):
        """Returns the backward reaction rate
           coefficient for reach reaction.
        
        INPUTS:
        =======
        kf:   list of floats
              Forward reaction rate coefficients
        T:    float
              Temperature of elementary reactions

        RETURNS:
        ========
        kf / kb: list of floats
              backward reaction rate coefficient
              for reach reaction

        EXAMPLES:
        =========
        """
        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.rxnset.nuij.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.rxnset.nuij.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction 
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / T

        # Ke
        kb = fact**self.gamma * np.exp(delta_G_over_RT)

        return kf / kb



class Rxnset():
    """Read and store certain NASA polynomial coefficients from the
    SQL database containing coefficients for all species.

    Parameters
    ==========
    species: list of string
             name of all species in the reaction syste
    nuij:    matrix
             the difference of product coefficients matrix and 
             reactant coefficients matrix

    Methods
    =======
    get_nasa_coefficients(T)

    read_nasa_coefficients()
    """
    def __init__(self, species, nuij):
        self.species = species
        self.nuij = nuij
        self.coefficients = self.read_nasa_coefficients()

    def get_nasa_coefficients(self, T):
        """Returns the corresponding NASA polynomial coefficients 
           for all species at the given temperature
        
        INPUTS:
        =======
        T: float
           Temperature of elementary reactions

        RETURNS:
        ========
        coeffs: nd numpy array
                NASA polynomial coefficients for all species
                at the given temperature

        EXAMPLES:
        =========
        """
        coeffs = []
        for species, data in self.coefficients.items():
            if T > data['low_T'] and T < data['mid_T']:
                coeffs.append(data['low_coeffs'])
            elif T >= data['mid_T'] and T < data['high_T']:
                coeffs.append(data['high_coeffs'])
            else:
                raise ValueError("Temperature not in range provided by NASA "
                                 "polynomial coefficients.")
        return np.array(coeffs)

    def read_nasa_coefficients(self):
        """Return NASA polynomial coefficients for all species.

        RETURNS:
        ========
        coeffs: ordered dictionary
                dictionary of dictionaries of species, each dictionary
                of species contains low, mid, high temperature, and NASA 
                polynomial coefficients for two temprature ranges 
                for all species

        EXAMPLES:
        =========
        """
        db = sqlite3.connect('data/NASA_polynomial_coefficients.sqlite')
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
