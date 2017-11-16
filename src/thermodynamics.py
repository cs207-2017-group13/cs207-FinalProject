"""
thermodynamics.py

victor: My feeling is that chemkin.py is getting to be too long.

We should be able to separate out some things into this file.

"""
import collections
import sqlite3
import numpy as np


class Thermochem():
    """Methods for calculating the backward reaction rate.
    Cp_over_R: Returns specific heat of each species given by
               the NASA polynomials.
    H_over_RT:  Returns the enthalpy of each species given by
                the NASA polynomials.
    S_over_R: Returns the entropy of each species given by
              the NASA polynomials.
    backward_coeffs:  Returns the backward reaction rate
                      coefficient for reach reaction.
    Please see the notes in each routine for clarifications and 
    warnings.  You will need to customize these methods (and 
    likely the entire class) to suit your own code base.  
    Nevertheless, it is hoped that you will find these methods 
    to be of some use.
    """

    def __init__(self, rxnset, p0=1e5, R=8.3144598):
        self.rxnset = rxnset
        self.p0 = p0 # Pa
        self.R = R # J / mol / K
        self.gamma = np.sum(self.rxnset.nuij, axis=0)

    def Cp_over_R(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.get_nasa_coefficients(T)

        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0 
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)

        return Cp_R

    def H_over_RT(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.get_nasa_coefficients(T)

        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)

        return H_RT
               

    def S_over_R(self, T):

        # WARNING:  This line will depend on your own data structures!
        # Be careful to get the correct coefficients for the appropriate 
        # temperature range.  That is, for T <= Tmid get the low temperature 
        # range coeffs and for T > Tmid get the high temperature range coeffs.
        a = self.rxnset.get_nasa_coefficients(T)

        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])

        return S_R

    def backward_coeffs(self, kf, T):

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


# This class could get a different name or be incorporated into an
# existing class.
class Rxnset():
    def __init__(self, species, nuij):
        # Does it make sense to have T as a parameter for the class?
        self.species = species
        self.nuij = nuij
        self.coefficients = self.read_nasa_coefficients()

    def get_nasa_coefficients(self, T):
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
