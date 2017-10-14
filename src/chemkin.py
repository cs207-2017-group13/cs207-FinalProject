import sys
import xml.etree.ElementTree as ET
import numpy as np
import numbers

IDEAL_GAS_CONSTANT = 8.314


class ReactionSystem():
    """

    """
    def __init__(self, elementary_reactions):
        

        self.k = elementary_reactions.get_k()
        react_coeff = elementary_reactions.get_react()


    def __repr__():
        print("Reaction rate coefficients: ", k, "\n",
              "Species: ", self.get_species(), "\n",
              "Conentration of Species: ", self.concs, "\n",
              "stoichiometric coefficients: ", self.nu_react, "\n",
              )


    def calculate_progress_rate(self, nu_react, concs):
        """Returns the progress rate of a system of irreversible, elementary reactions

        INPUTS:
        =======
        nu_react: numpy array of floats, 
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reaction
        k:        array of floats
                  Reaction rate coefficient for the reaction
        concs:    numpy array of floats 
                  concentration of species

        RETURNS:
        ========
        omega: numpy array of floats
               size: num_reactions
               progress rate of each reaction

        EXAMPLES:
        =========
        >>> calculate_progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([2.0, 1.0, 1.0]), 10.0)
        array([ 40.,  20.])
        """
        if type(self.k) is float:
            self.k = [self.k]*len(nu_react[0])
        if len(nu_react) != len(concs) or len(nu_react)==0 or len(nu_react[0])!=len(self.k):
            raise  ValueError("Invalid input parameters.")
        progress = self.k # Initialize progress rates with reaction rate coefficients
        for jdx, rj in enumerate(progress):
            if rj < 0:
                raise ValueError("k = {0:18.16e}:  Negative reaction rate coefficients are prohibited!".format(rj))
            for idx, xi in enumerate(concs):
                nu_ij = nu_react[idx,jdx]
                if xi  < 0.0:
                    raise ValueError("x{0} = {1:18.16e}:  Negative concentrations are prohibited!".format(idx, xi))
                if nu_ij < 0:
                    raise ValueError("nu_{0}{1} = {2}:  Negative stoichiometric coefficients are prohibited!".format(idx, jdx, nu_ij))
                
                progress[jdx] *= xi**nu_ij
        return progress      

    def calculate_reaction_rate(self, nu_react, nu_prod, concs):
        """Returns the reaction rate of a system of irreversible, elementary reactions
        
        INPUTS:
        =======
        nu_react: numpy array of floats, 
                  size: num_species X num_reactions
                  stoichiometric coefficients for the reactants
        nu_prod:  numpy array of floats, 
                  size: num_species X num_reactions
                  stoichiometric coefficients for the products
        concs:    numpy array of floats 
                  concentration of species

        RETURNS:
        ========
        f: numpy array of floats
           size: num_species
           reaction rate of each specie
        
        EXAMPLES:
        =========
        >>> nu_react = np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
        >>> nu_prod = np.array([[0.0, 1.0], [0.0, 2.0], [1.0, 0.0]])
        >>> r = np.array([ 40.,  20.])
        >>> calculate_reaction_rate(nu_react, nu_prod, r)
        array([-20., -40.,   0.])
        """
        if len(nu_react)==0 or len(nu_react)!=len(nu_prod) or len(nu_react[0])!=len(nu_prod[0]):
            raise  ValueError("Invalid input parameters.")
        nu = nu_prod - nu_react
        rj = self.calculate_progress_rate(nu_react, concs)
        return np.dot(nu, rj)

    def get_species():
        pass



class ElementaryReaction():
    def calculate_reaction_coefficients(temperature):
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

