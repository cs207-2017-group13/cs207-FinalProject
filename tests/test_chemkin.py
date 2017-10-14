import chemkin

def test_len():
	reader = chemkin.XMLReader("test_data1.xml")
	reaction_system = reader.get_reaction_systems()
	assert len(reaction_system[0]) == 5

def test_progress_rate():
	reader = chemkin.XMLReader("test_data1.xml")
	reaction_system = reader.get_reaction_systems()
	result = [2.5589111307566812, 1110.2037988957545, 0.0070156249238785568]
	assert reaction_system[0].calculate_progress_rate([3., 1., 2., 3., 1.], 425) == result

def test_neg_concentration():
	reader = chemkin.XMLReader("test_data1.xml")
	reaction_system = reader.get_reaction_systems()
	try:
		reaction_system[0].calculate_progress_rate([3., -1., 2., -3., 1.], 425)
	except ValueError as err:
		assert(type(err) == ValueError)

def test_reaction_rate():
	reader = chemkin.XMLReader("test_data1.xml")
	reaction_system = reader.get_reaction_systems()
	result = [ 1107.64488776, -1107.64488776,  1112.77674128, -1110.21081452, -2.56592676]
	assert reaction_system[0].calculate_reaction_rate([3., 1., 2., 3., 1.], 425) == result

