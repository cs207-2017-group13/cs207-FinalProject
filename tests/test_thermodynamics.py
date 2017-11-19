import chemkin.thermodynamics as thermodynamics
import chemkin.chemkin as chemkin
import numpy as np


def test_get_nasa_coefficients():
    reader = chemkin.XMLReader("tests/rxns_reversible.xml")
    reaction_system = reader.get_reaction_systems()[0]
    nu = (reaction_system.product_coefficients
          - reaction_system.reactant_coefficients)
    rxnset = thermodynamics.Rxnset(reaction_system.species, nu)
    coefficients = rxnset.get_nasa_coefficients(900)
    assert coefficients.shape == (8, 7)
    assert np.all(coefficients[0] == np.array(
        [2.50000000e+00, 7.05332819e-13, -1.99591964e-15, 2.30081632e-18,
         -9.27732332e-22, 2.54736599e+04, -4.46682853e-01]))


def test_temp_range():
    reader = chemkin.XMLReader("tests/rxns_reversible.xml")
    reaction_system = reader.get_reaction_systems()[0]
    nu = (reaction_system.product_coefficients
          - reaction_system.reactant_coefficients)
    try:
        rxnset = thermodynamics.Rxnset(reaction_system.species, nu)
        rxnset.get_nasa_coefficients(50)
    except ValueError as err:
        assert(type(err) == ValueError)
