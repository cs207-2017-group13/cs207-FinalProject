import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import chemkin.ode_solver as ode_solver


AVOGADRO = 6.022e23


class PropensityZeroException:
    pass


class ReactionSimulator():
    """Base class of simulators

    Methods
    -------
    plot()
    """
    def plot(self):
        """Plot abundances/concentrations of species
           over times

        Returns
        -------
        plot : plot of abundances/concentrations of species
               over times
        """
        plt.plot(self.times, self.abundances)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(self.reaction_system.species, loc='best')
        # file_name = "examples/figures/"+name+".png"
        # plt.savefig(file_name, dpi=125)
        plt.show()


class StochasticSimulator(ReactionSimulator):
    """Carries out stochastic reaction simulations.

    Inherits from base class `ReactionSimulator`.

    Note that a reversible elementary reaction represents two
    different reactions in a stochastic simulation.

    Parameters
    ----------
    reaction_system : ReactionSystem
        Containing all the methods and attributes of 
        `ReactionSystem` class
    initial_abundances : list
        Initial abundances of each species inside
        the reaction_system
    temperature : float
        Temperature when the reaction take place
    t_span : tuple
        Time span of the reactions users want to 
        simulate
    system_volume : float
        blah blahx

    Methods
    -------
    calculate_state_change_matrix()
        blah
    calculate_stochastic_constants(temperature)
        blah
    calculate_reaction_propensities()
        blah
    simulate()
        blah

    Examples
    --------

    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 t_span, system_volume):
        self.reaction_system = reaction_system
        self.abundances = [np.array(initial_abundances)]
        self.times = [t_span[0]]
        self.t_span = t_span
        self.temperature = temperature
        self.system_volume = system_volume
        self.state_change_matrix = self.calculate_state_change_matrix()
        self.stochastic_constants = self.calculate_stochastic_constants(
            temperature)
        # self.reaction_propensities = self.calculate_reaction_propensities()

    def calculate_state_change_matrix(self):
        """Set vectors that determine how abundances change.

        A 2D matrix with dimension n_reactions x n_species is
        calculated. Each row gives the change in abundance if the
        reaction corresponding to the row happens as the stochastic
        event.

        Returns
        -------
        state_change_matrix : np.ndarray
            A 2D matrix with dimension n_reactions x n_species.

        """
        n_species = len(self.reaction_system)
        state_change_matrix = []
        for reaction in self.reaction_system.elementary_reactions:
            state_change_vector = np.zeros(n_species, dtype=int)
            for i, species in enumerate(self.reaction_system.species):
                state_change_vector[i] = (
                    reaction.get_products().get(species, 0)
                    - reaction.get_reactants().get(species, 0))
            state_change_matrix.append(state_change_vector)
            if reaction.reversible:
                state_change_matrix.append(-1*state_change_vector)
        return np.array(state_change_matrix)

    def calculate_stochastic_constants(self, temperature):
        """Compute stochastic rate constants.

        stochastic constant * abundance * dt gives probability.

        """
        stochastic_constants = []
        # First, obtain deterministic rate constants
        rate_constants = self.reaction_system.get_rate_coefficients(
            temperature)
        backward_rate_constants = (
            self.reaction_system.get_backward_rate_coefficients(
                self.temperature))
        for forward_rate, backward_rate, reaction in zip(
                rate_constants, backward_rate_constants,
                self.reaction_system.elementary_reactions):
            reaction_order = reaction.calculate_reaction_order()
            for rate, order in zip(
                    [forward_rate, backward_rate], reaction_order):
                if order == 1:
                    stochastic_constant = rate
                elif order == 2:
                    stochastic_constant = rate / AVOGADRO / self.system_volume
                else:
                    raise NotImplementedError("Reaction order %d not valid for"
                                              " stochastic simulation" % order)
                stochastic_constants.append(stochastic_constant)
        return stochastic_constants

    def calculate_reaction_propensities(self):
        """
        Reaction propensities are determined by

        """
        reaction_propensities = []
        for state_change_vector, stochastic_constant in zip(
                self.state_change_matrix, self.stochastic_constants):
            propensity = stochastic_constant
            for change, species_abundance in zip(
                    state_change_vector, self.abundances[-1]):
                if change >= 0:
                    continue
                elif change == -1:
                    propensity *= species_abundance
                elif change == -2:
                    propensity *= species_abundance * (species_abundance - 1)
            reaction_propensities.append(propensity)
        return np.array(reaction_propensities)

    def simulate(self, seed=None):
        np.random.seed(seed)
        while self.times[-1] < self.t_span[-1]:
            try:
                self._advance_simulation()
            except PropensityZeroException:
                pass
        return

    # this is the simplest algorithm
    def _advance_simulation(self):
        reaction_propensities = self.calculate_reaction_propensities()
        print(reaction_propensities)
        propensity_cumsum = np.cumsum(reaction_propensities)
        propensity_sum = propensity_cumsum[-1]
        if propensity_sum == 0.:
            raise PropensityZeroException
        r1 = np.random.rand()
        r2 = np.random.rand()*propensity_sum
        self.times.append(self.times[-1] + ((1/propensity_sum) * np.log(1/r1)))
        reaction_index = np.where(propensity_cumsum > r2)[0][0]
        state_change_vector = self.state_change_matrix[reaction_index]
        self.abundances.append(self.abundances[-1] + state_change_vector)


class DeterministicSimulator(ReactionSimulator):
    """Class for deterministic simulation
    Inherits from base class `ReactionSimulator`.
    Simulate species abundances deterministically.

    Parameters
    ----------
    reaction_system :    `ReactionSystem` class object
                        Containing all the methods and attributes of
                        `ReactionSystem` class
    initial_abundances : list
                        Initial abundances of each species inside
                        the reaction_system
    temperature :        float
                        Temperature when the reaction take place
    t_span :             tuple of floats
                        Time span of the reactions users want to 
                        simulate
    dt :                 float
                        size of time steps users want to simulate

    Methods
    -------
    simulate(method='bdf', epsilon = 1e-06)
    diff_func(t, y)
    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 t_span, dt=0.01):
        self.reaction_system = reaction_system
        if temperature <=0:
            raise ValueError("Temperature must be positive.")
        self.temperature = temperature
        if len(initial_abundances) != len(self.reaction_system.species):
            raise ValueError("Invalid initial species abundances.")
        self.times = [t_span[0]]
        self.abundances = [initial_abundances]
        self.t_span = t_span
        self.dt = dt
        self.ode_integrator = ode_solver.ODE_solver(
            self.diff_func, initial_abundances, t_span, self.dt, False)

    def simulate(self, method='bdf', epsilon = 1e-06):
        """Simulate species abundances deterministically.

        We implemented three methods to solve the ordinary differential 
        equation. Backward differentiation formula and backward euler
        are good for stiff functions, and rk45 is accurate for non-stiff
        functions. BDF is the most suitable for solving ODE problems in
        chemical kinetics, so our default method is set as BDF.

        Parameters
        ----------
        method :  string
                 default as 'bdf'-- backward differentiation formula
                 name of the ODE solver
        epsilon : float
                 default as 1e-06
                 tolerance of error

        Returns
        -------
        self.times :      array
                         time of evaluations
        self.abundances : array
                         abundances of species at every time step

        Examples
        --------
        >>> import chemkin.chemkin as chemkin
        >>> concs = np.array([1., 2., 1., 3., 1.])*1e-05
        >>> reader = chemkin.XMLReader("tests/rxns.xml")
        >>> reaction_system = reader.get_reaction_systems()[0]
        >>> det_sim = DeterministicSimulator(reaction_system, concs, 800, [0, 0.01], dt=0.01)
        >>> det_sim.simulate()
        ([0, 0.01], [array([  1.00000000e-05,   2.00000000e-05,   1.00000000e-05,
                 3.00000000e-05,   1.00000000e-05]), array([  1.03778637e-05,   1.96221363e-05,   1.03927707e-05,
                 2.96146828e-05,   9.99254648e-06])])
        """
        choices = ['backward_euler','rk45', 'bdf']
        if method not in choices:
            raise ValueError("Wrong method.")

        if method == 'bdf':
            self.times, self.abundances = self.ode_integrator.BDF(epsilon)
        elif method == 'rk45':
            self.times, self.abundances = self.ode_integrator.rk45(epsilon)
        else:
            self.times, self.abundances = self.ode_integrator.backward_euler(epsilon)
        return self.times, self.abundances

    def diff_func(self, t, y):
        ''' Turn the calculate reaction rate function into the rhs of
            the ODE.

        In order for the ode solver to work, we need a function with
        t (time) and y as parameters. However, the original calculate reaction
        rate function is a function of y and temperature. We need to
        transform the original function.

        Parameters
        ----------
        t : float
          time
        y : array
           value of the variable we would like to integrate

        Returns
        -------
        function : right hand side of the ODE function
        '''
        if (y<0).any():
            y = [0 if i<0 else i for i in y]
        return self.reaction_system.calculate_reaction_rate(y, self.temperature)


if __name__ == "__main__":
    import doctest
    doctest.testmod()       
