class ReactionSimulator():
    def plot(self, t, y, species):
        # plot y's 
        for i in range(len(species)):
            plt.plot(t, y[i], label=species[i])
            plt.xlabel("Time")
            plt.ylabel("Concentration")
            plt.legend(loc='best')
            plt.show()
            # maybe save figure?


class StochasticSimulator(ReactionSimulator):
    pass


class DeterministicSimulator(ReactionSimulator):
    def __init__(self, reaction_system, initial_concentrations, temperature):
        self.reaction_system = reaction_system
        if temperature <=0:
            raise ValueError("Temperature must be positive.")
        self.temperature = temperature
        # assert concentrations shape is number of species in reaction system
        self.time = [0.]
        self.concentrations = [initial_concentrations]

    def simulate(self, t_final, method='rk45'):
        """Calculate concentrations between last time and `t_final`.

        """
        # initialize solver
        while t < t_final:
            self.step()
        # save data

    def diff_func(self, t, y):
        self.reaction_system.calculate_reaction_rate(y, self.temperature)

        
