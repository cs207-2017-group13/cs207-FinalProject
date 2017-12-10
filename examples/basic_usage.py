import chemkin.chemkin as chemkin

# irreversible reactions
concs = np.array([1., 2., 1., 3., 1.])*1e-05
reader = chemkin.XMLReader("tests/rxns.xml")
reaction_system = reader.get_reaction_systems()[0]
react_rate = reaction_system.calculate_reaction_rate(concs, 800)

# deterministic simulation
det_sim = reaction_system.setup_reaction_simulator('deterministic', concs, 800, [0, 10], dt=0.01)
t, y = det_sim.simulate()
det_sim.plot_simulation()

# stochastic simulation
abundances = [10., 10., 10., 10., 10.]
stoch_sim = reaction_system.setup_reaction_simulator('stochastic', abundances, 800, [0, 100], 1e-15)
stoch_sim.simulate()
stoch_sim.plot_simulation()

# reversible reactions
concs = [1, 2, 3, 1, 2, 3, 1, 2]
reader = chemkin.XMLReader("tests/rxns_reversible.xml")
reaction_system = reader.get_reaction_systems()[0]
react_rate = reaction_system.calculate_reaction_rate(concs, 800)


