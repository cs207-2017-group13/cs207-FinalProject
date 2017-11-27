class ReactionSimulator():
    # could contain common plotting
    pass


class StochasticSimulator(ReactionSimulator):
    pass


class DeterministicSimulator(ReactionSimulator):
    def __init__(self, reaction_system, initial_concentrations, temperature):
        self.reaction_system = reaction_system
        # self.concentrations = initial_concentrations
        self.temperature = temperature
        # assert concentrations shape is number of species in reaction system
        # assert temperature is a positive number
        self.time = [0.]
        self.concentrations = [initial_concentrations]

    def simulate(self, t_final, method='rk45'):
        """Calculate concentrations between last time and `t_final`.

        """
        # initialize solver
        while t < t_final:
            self.step()
        # save data

    def step(self, dt):
        # placeholder
        rate = self.reaction_system.calculate_reaction_rate()
        previous_concentrations = self.concentrations[-1]
        new_concentrations = dt * rate + previous_concentrations

        
