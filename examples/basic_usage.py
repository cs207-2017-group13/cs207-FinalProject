import chemkin.chemkin as chemkin
import numpy as np

# irreversible reactions
concs = np.array([1., 2., 1., 3., 1.])*1e-05
reader = chemkin.XMLReader("tests/rxns.xml")
reaction_system = reader.get_reaction_systems()[0]
react_rate = reaction_system.calculate_reaction_rate(concs, 800)

# reversible reactions
abundances = np.array([10., 10., 10., 10., 10., 10., 10., 10.])
concs = abundances/6.02e23/1e-15
reader = chemkin.XMLReader("tests/rxns_reversible.xml")
reaction_system = reader.get_reaction_systems()[0]
react_rate = reaction_system.calculate_reaction_rate(concs, 800)

# deterministic simulation
det_sim = reaction_system.setup_reaction_simulator('deterministic', concs, 400, [0, 1e-4], dt=1e-7)
det_sim.simulate(epsilon=1e-9)

# stochastic simulation
stoch_sim = reaction_system.setup_reaction_simulator('stochastic', abundances, 400, [0, 1e-4], 1e-15)
stoch_sim.simulate()

# compare the plots of deterministic simulation and stochastic simulation
det_sim.plot_simulation()
stoch_sim.plot_simulation()


