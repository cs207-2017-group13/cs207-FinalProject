import thermodynamics
import chemkin

def test_temp_range():
    reader = chemkin.XMLReader("tests/rxns_reversible.xml")
    reaction_system = reader.get_reaction_systems()
    nu = reaction_system[0].product_coefficients - reaction_system[0].reactant_coefficients
    try:
        rxnset = thermodynamics.Rxnset(reaction_system[0].species, nu)
        rxnset.get_nasa_coefficients(50)
    except ValueError as err:
        assert(type(err) == ValueError)
