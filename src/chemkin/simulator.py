import chemkin.ode_solver as ode_solver
# import ode_solver
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


class ReactionSimulator():
    def __init__(self):
        # need init function here?
        pass

    def plot(self, t, y, species):
        # plot y's
        y = np.array(self.concentrations)
        y = y.transpose()
        for i, species_name in enumerate(self.reaction_system.species):
            plt.plot(self.times, y[i], label=species_name)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(loc='best')
        plt.show()
        # maybe save figure?


class StochasticSimulator(ReactionSimulator):
    """Carries out stochastic reaction simulations.

    Inherits from base class `ReactionSimulator`.

    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 system_volume):
        self.reaction_system = reaction_system
        self.abundances = [initial_abundances]
        self.times = [0.]
        self.temperature = temperature
        self.system_volume = system_volume
        self.state_change_vectors = self.calculate_state_change_vectors()
        self.reaction_propensities = self.calculate_reaction_propensities(
            temperature)

    def calculate_state_change_vectors(self):
        """Set vectors that determine how abundances change.

        `self.state_change_vectors` is a 2D matrix of n_reactions x
        n_species. Note that this is flipped from other places in this code.

        """
        n_reactions = len(self.reaction_system.elementary_reactions)
        n_species = len(self.reaction_system)
        self.state_change_vectors = (
            self.reaction_system.product_coefficients
            - self.reaction_system.reactant_coefficients).T
        assert self.state_change_vectors.shape == (n_reactions, n_species)
        return

    def calculate_reaction_propensities(self, temperature):
        """

        Reaction propensity * dt gives probability.

        """
        # First, obtain deterministic rate constants
        rate_constants = self.reaction_system.get_rate_coefficients(
            temperature)
        backward_rate_constants = (
            self.reaction_system.get_backward_rate_coefficients())
        # Obtain order of each reaction
        pass

    def simulate_system(self, t_final, seed=None):
        np.random.seed(seed)
        while self.times[-1] < t_final:
            self.advance_simulation()
        return

    # I am imagining multiple possible methods
    def advance_simulation(self):
        r1 = np.random.rand()
        r2 = np.random.rand()
        pass


class DeterministicSimulator(ReactionSimulator):
    def __init__(self, reaction_system, initial_concentrations, temperature, t_span, dt=0.1):
        self.reaction_system = reaction_system
        if temperature <=0:
            raise ValueError("Temperature must be positive.")
        self.temperature = temperature
        if len(initial_concentrations) != len(self.reaction_system.species):
            raise ValueError("Invalid initial concentration.")
        self.times = [t_span[0]]
        self.concentrations = [initial_concentrations]
        self.t_span = t_span
        self.dt = dt
        self.ode_solver = ode_solver.ODE_solver(
            self.diff_func, initial_concentrations, t_span, self.dt)

    def simulate(self, method='bdf', epsilon = 1e-06):
        choices = ['backward_euler','rk45', 'bdf']
        if method not in choices:
            raise ValueError("Wrong method.")

        if method == 'bdf':
            self.times, self.concentrations = self.ode_solver.BDF(epsilon)
        elif method == 'rk45':
            self.times, self.concentrations = self.ode_solver.rk45(epsilon)
        else:
            self.times, self.concentrations = self.ode_solver.backward_euler(epsilon)
        return self.times, self.concentrations

    def simulate_step(self, method='backward_euler', epsilon = 1e-06):
        """Calculate concentrations between last time and `t_final`
        """
        # initialize solver
        if self.times[-1] < self.t_span[-1]:
            choices = ['backward_euler','rk45']
            if method not in choices:
                raise ValueError("Wrong method.")
            if method == 'rk45':
                message = self.ode_solver.rk45_step(epsilon)
            else:
                message = self.ode_solver.backward_euler_step(epsilon)
            self.dt, self.time, self.concentrations = self.ode_solver.dt, self.ode_solver.t, self.ode_solver.y
        else:
            raise IndexError("Time exceeds time span.")
        return message

    def diff_func(self, t, y):
        self.reaction_system.calculate_reaction_rate(y, self.temperature)

        
