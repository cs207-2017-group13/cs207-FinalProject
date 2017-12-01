#!/usr/bin/env python
import chemkin.ode_solver as ode_solver
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
            plt.plot(self.time, y[i], label=species_name)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(loc='best')
        plt.show()
        # maybe save figure?


class StochasticSimulator(ReactionSimulator):
    def __init__(self, reaction_system, initial_abundances, temperature):
        self.reaction_system = reaction_system
        self.abundances = [initial_abundances]
        self.temperature = temperature
        self.state_change_vectors = self.calculate_state_change_vectors()
        self.reaction_propensities = self.calculate_reaction_propensities(
            temperature)

    def calculate_state_change_vectors(self):
        pass

    def calculate_reaction_propensities(self, temperature):
        # First, obtain deterministic rate constants
        rate_constants = self.reaction_system.get_rate_coefficients(
            temperature)
        backward_rate_constants = self.reaction_system.get_backward_rate_coefficients()
        # Obtain order of each reaction
        pass

    def some_sort_of_simulator(self):
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

        
