import chemkin.chemkin as chemkin
import numpy as np

def test_for_keys():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    try:
        elementary_reaction =  chemkin.ElementaryReaction(reaction_properties)
    except KeyError as err:
        assert (type(err) == KeyError)
        
    
def test_get_reactants():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    elementary_reaction =  chemkin.ElementaryReaction(reaction_properties)
    assert elementary_reaction.get_reactants() == {'H': '1', 'O2': '1'}


def test_get_products():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
    assert elementary_reaction.get_products() == {'O': '1' , 'OH': '1'}

def test_rate_type_exists():
        reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : '' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}

        try:
            elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
        except ValueError as err:
            assert (type(err) == ValueError)

            
def test_ratetype_valid():
        reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Ran' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}

        try:
            elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
            elementary_reaction.calculate_rate_coefficient(1000)
        except NotImplementedError as err:
            assert (type(err) == NotImplementedError)
            
def test_calculate_const_coeff():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'constant' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    
    elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
    assert elementary_reaction.calculate_rate_coefficient(1000) == 1.0

def test_calculate_arr_coeff():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
        
    elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
    assert elementary_reaction.calculate_rate_coefficient(1000) == 5.2102032610668552

def test_A():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : -3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    try:    
        elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
        elementary_reaction.calculate_rate_coefficient(1000)
    except ValueError as err:
        assert (type(err) == ValueError)

def test_E():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : -71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'no', 'type' : ' Elementary'}
    try:    
        elementary_reaction = chemkin.ElementaryReaction(reaction_properties)
        elementary_reaction.calculate_rate_coefficient(1000)
    except ValueError as err:
        assert (type(err) == ValueError)

def test_if_elementary():
    reader = chemkin.XMLReader("tests/test_data_elementary.xml")
    try:    
        reaction_system = reader.get_reaction_systems()
    except NotImplementedError as err:
        assert (type(err) == NotImplementedError)

def test_len():
    reader = chemkin.XMLReader("tests/test_data1.xml")
    reaction_system = reader.get_reaction_systems()
    assert len(reaction_system[0]) == 5

def test_progress_rate():
    reader = chemkin.XMLReader("tests/test_data1.xml")
    reaction_system = reader.get_reaction_systems()
    result = [2.5589111307566812, 1110.2037988957545, 0.0070156249238785568]
    assert reaction_system[0].calculate_progress_rate([3., 1., 2., 3., 1.], 425) == result

def test_neg_concentration():
    reader = chemkin.XMLReader("tests/test_data1.xml")
    reaction_system = reader.get_reaction_systems()
    try:
        reaction_system[0].calculate_progress_rate([3., -1., 2., -3., 1.], 425)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_neg_stoich_coeff():
    reader = chemkin.XMLReader("tests/rxns_negative_coeff.xml")
    reaction_system = reader.get_reaction_systems()
    try:
        reaction_system[0].calculate_progress_rate([1, 2, 3, 1, 2, 3, 1, 2], 800)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reversible_reaction_rate():
    reader = chemkin.XMLReader("tests/rxns_reversible.xml")
    reaction_system = reader.get_reaction_systems()
    assert np.allclose(list(reaction_system[0].calculate_reaction_rate([1, 2, 3, 1, 2, 3, 1, 2], 800)), 
        [2.3553537901231184e+17, -2.2321034840197091e+17, -2.6059674644924278e+17, 37395259554341.125, 
        12661546307481366.0, 2.3584012306405293e+17, -198181887233104.84, -69166904953671.859])
