#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import chemkin.ode_solver as ode_solver


AVOGADRO = 6.022e23


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

    Note that a reversible elementary reaction represents two
    reactions in stochastic simulation.

    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 system_volume):
        self.reaction_system = reaction_system
        self.abundances = [initial_abundances]
        self.times = [0.]
        self.temperature = temperature
        self.system_volume = system_volume
        self.state_change_matrix = self.calculate_state_change_vectors()
        self.reaction_propensities = self.calculate_reaction_propensities(
            temperature)

    def calculate_state_change_vectors(self):
        """Set vectors that determine how abundances change.

        `self.state_change_matrix` is a 2D matrix of n_reactions x
        n_species.

        """
        n_species = len(self.reaction_system)
        state_change_vector = np.zeros(n_species, dtype=int)
        state_change_matrix = []
        for reaction in self.reaction_system.elementary_reactions:
            for i, species in enumerate(self.reaction_system.species):
                state_change_vector[i] = (
                    reaction.get_products()[species]
                    - reaction.get_reactants()[species])
            state_change_matrix.append(state_change_vector)
            if reaction.reversible:
                state_change_matrix.append(-1*state_change_vector)
        return np.array(state_change_matrix)

    def calculate_reaction_propensities(self, temperature):
        """

        Reaction propensity * dt gives probability.

        """
        reaction_propensities = []
        # First, obtain deterministic rate constants
        rate_constants = self.reaction_system.get_rate_coefficients(
            temperature)
        backward_rate_constants = (
            self.reaction_system.get_backward_rate_coefficients())
        for forward_rate, backward_rate, reaction in zip(
                rate_constants, backward_rate_constants,
                self.reaction_system.elementary_reactions):
            reaction_order = reaction.calculate_reaction_order()
            for rate, order in zip(
                    [forward_rate, backward_rate], reaction_order):
                if order == 1:
                    reaction_propensity = rate
                elif order == 2:
                    reaction_propensity = rate / AVOGADRO / self.system_volume
                else:
                    raise NotImplementedError
                reaction_propensities.append(reaction_propensity)
        return reaction_propensities

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
    def __init__(self, reaction_system, initial_concentrations, temperature, t_span, dt=1):
        self.reaction_system = reaction_system
        if temperature <=0:
            raise ValueError("Temperature must be positive.")
        self.temperature = temperature
        if len(initial_concentrations) != len(self.reaction_system.species):
            raise ValueError("Invalid initial concentration.")
        self.time = [t_span[0]]
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
            self.time, self.concentrations = self.ode_solver.BDF(epsilon)
        elif method == 'rk45':
            self.time, self.concentrations = self.ode_solver.rk45(epsilon)
        else:
            self.time, self.concentrations = self.ode_solver.backward_euler(epsilon)

    def simulate_step(self, method='backward_euler', epsilon = 1e-06):
        """Calculate concentrations between last time and `t_final`
        """
        # initialize solver
        if self.time[-1] < self.t_span[-1]:
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

        
