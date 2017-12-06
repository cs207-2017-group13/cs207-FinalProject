import chemkin.simulator as simulator
import chemkin.chemkin as chemkin
import numpy as np


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
    concs = np.array([1., 2., 1., 3., 1.])*1e-05
    reader = chemkin.XMLReader("tests/rxns.xml")
    reaction_system = reader.get_reaction_systems()[0]
    det_sim = simulator.DeterministicSimulator(reaction_system,
        concs, [0, 0.01], dt=0.01)
    t, y = det_sim.simulate('backward_euler')
    y0 = y[0]
    y1 = y[1]
    assert np.allclose(t, [0, 0.01])
    assert np.allclose(y0, np.array([1.00000000e-05, 2.00000000e-05, 1.00000000e-05,
        3.00000000e-05, 1.00000000e-05]))
    assert np.allclose(y1, np.array([1.03940908e-05, 1.96059092e-05, 1.04083239e-05,
        2.95987926e-05, 9.99288345e-06]))
