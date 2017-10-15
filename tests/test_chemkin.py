import chemkin

def test_for_reversible():
    reaction_properties = {'equation' : 'H + O2  [=] OH + O' ,
                             'id' : 'reaction01', 'products' : {'O' : '1' , 'OH' : '1'}, 
                             'rate_params' : {'A' : 3520000.0, 'E' : 71400.0 , 'b' : -0.7 }, 
                             'rate_type' : 'Arrhenius' , 'reactants' : {'H' : '1' , 'O2' : '1'} , 
                             'reversible': 'yes', 'type' : ' Elementary'}
    try:
        elementary_reaction =  chemkin.ElementaryReaction(reaction_properties)
    except NotImplementedError as err:
        assert (type(err) == NotImplementedError)
    
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