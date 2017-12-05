import chemkin.simulator as simulator
import chemkin.chemkin as chemkin


def test_neg_temp():
	concs = [1., 2., 1., 3., 1.]
	reader = chemkin.XMLReader("tests/rxns.xml")
	reaction_system = reader.get_reaction_systems()[0]
	try:
		det_sim = simulator.DeterministicSimulator(reaction_system, 
												   concs, -800, [0, 2], dt=1)
	except ValueError as err:
		assert (type(err) == ValueError)

def test_neg_temp():
	concs = [1., 2., 1., 3.]
	reader = chemkin.XMLReader("tests/rxns.xml")
	reaction_system = reader.get_reaction_systems()[0]
	try:
		det_sim = simulator.DeterministicSimulator(reaction_system, 
												   concs, 800, [0, 2], dt=1)
	except ValueError as err:
		assert (type(err) == ValueError)

def test_sim_backward_euler():
	concs = [1., 2., 1., 3., 1.]
	reader = chemkin.XMLReader("tests/rxns.xml")
	reaction_system = reader.get_reaction_systems()[0]
	det_sim = simulator.DeterministicSimulator(reaction_system,
											   concs, 800, [0, 2], dt=1)
	assert det_sim.simulate('backward_euler') == 1
