import numpy as np
import matplotlib.pyplot as plt
import chemkin.ode_solver as ode_solver


AVOGADRO = 6.022e23


class PropensityZeroException(Exception):
    pass


class ReactionSimulator():
    """Base class of simulators.

    Methods
    -------
    save_data(path)

    """
    def save_data(self, path):
        """Save simulation data to file.

        Simulation data will be saved as space-separated text file.

        Parameters
        ----------
        path : str
            Location to save data.

        """
        data = np.concatenate([np.array([self.times]).T,
                               np.array(self.abundances)], axis=1)
        header = ",".join(["Times"] + self.reaction_system.species)
        if isinstance(self.abundances[0][0], int):
            format = ["%.8e"] + ["%d"] * len(self.reaction_system)
        else:
            format = ["%.8e"] + ["%.3e"] * len(self.reaction_system)
        np.savetxt(path, data, fmt=format, header=header)

    def _validate_arguments(self):
        """Ensure simulation parameters agree with `ReactionSystem`."""
        try:
            float(self.temperature)
        except:
            raise ValueError("Temperature not a number: %s" % self.temperature)
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive.")
        if len(self.abundances[0]) != len(self.reaction_system):
            raise ValueError("Number of provided abundances (%d) does not "
                             "match number of chemical species (%d)." % (
                                 len(self.abundances[0]),
                                 len(self.reaction_system)))
        assert self.t_span[1] > self.t_span[0]


class StochasticSimulator(ReactionSimulator):
    """Carries out stochastic reaction simulations.

    Inherits from base class `ReactionSimulator`.

    Note that a reversible elementary reaction represents two
    different reactions in a stochastic simulation.

    Parameters
    ----------
    reaction_system : ReactionSystem
        `ReactionSystem` instance with all the details regarding
        reaction to simulate.
    initial_abundances : list
        Initial abundances of each species inside the
        reaction_system. Abundances should be discrete numbers
        (integers).
    temperature : float
        Temperature at which the reaction take place.
    t_span : tuple
        Time span to simulate: (begin, end)
    system_volume : float
        The volume of the reaction system, in liters.

    Methods
    -------
    calculate_state_change_matrix()
        Determine how abundances change with each reaction event.
    calculate_stochastic_constants(temperature)
        Determine stochastic rate constants from deterministic rate
        constants.
    calculate_reaction_propensities()
        Determine the propensity of each reaction, the probability of
        the reaction to occur in the next interval [t, t+dt).
    simulate()
        Run stochastic simulation between `t_span[0]` and `t_span[1]`.

    Examples
    --------
    >>> # From `ReactionSystem`
    >>> import chemkin.chemkin
    >>> reader = chemkin.chemkin.XMLReader("tests/rxns.xml")
    >>> reaction_system = reader.get_reaction_systems()[0]
    >>> stochastic_simulator = reaction_system.setup_reaction_simulator(
    ... "stochastic", [10, 10, 10, 10, 10], 800., (0, 30000))

    >>> # From constructor
    >>> stochastic_simulator = StochasticSimulator(reaction_system,
    ... [10, 10, 10, 10, 10], 800., (0, 30000), 1e-15)

    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 t_span, system_volume):
        self.reaction_system = reaction_system
        self.abundances = [np.array(initial_abundances)]
        self.times = [t_span[0]]
        self.t_span = t_span
        self.temperature = temperature
        assert isinstance(system_volume, float)
        self.system_volume = system_volume
        self._validate_arguments()
        self.state_change_matrix = self.calculate_state_change_matrix()
        self.stochastic_constants = self.calculate_stochastic_constants()
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

    def calculate_stochastic_constants(self):
        """Compute stochastic rate constants (coefficients).

        Stochastic rate constants are calculated from deterministicc
        reaction rate constants. First-order stochastic rate constants
        are equivalent to stochastic rate constants, but second-order
        rate constants must be divided by the system volume and
        converted from moles.

        Reactions of other orders are not elementary reactions and are
        not supported.

        Returns
        -------
        stochastic_constants : list
            Stochastic constants, one for each reaction.

        """
        stochastic_constants = []
        # First, obtain deterministic rate constants
        rate_constants = self.reaction_system.get_rate_coefficients(
            self.temperature)
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
        """Determine probabilities for each reaction.

        The reaction propensity * dt gives the probability that a
        particular reaction will occur in the time interval [t, t+dt).

        See `manual/model_document.pdf` for scientific notes.

        Returns
        -------
        reaction_propensities : np.ndarray (1D)
            Reaction propensities, one for each reaction.

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
        """Run stochastic simulation.

        Reaction runs from time `t_span[0]` to `t_span[1]`. Times are
        put in `self.times`, and abundances at each time are put in
        `self.abundances`.

        Simulation will also stop when reaction propensities are all
        zero. This can happen in a system with irreversible reactions
        in which no reactants are remaining.

        Parameters
        ----------
        seed : int
            Sets numpy RNG seed. Default of `None` will set seed
            randomly.
        """
        np.random.seed(seed)
        while self.times[-1] < self.t_span[-1]:
            try:
                self._advance_simulation()
            except PropensityZeroException:
                break
        return self.times, self.abundances

    def _advance_simulation(self):
        """Run a single stochastic simulation step.

        Uses the basic Gillespie stochastic simulation algorithm.
        """
        reaction_propensities = self.calculate_reaction_propensities()
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

    def plot_simulation(self, show=True, savefig=None, title=None):
        """Plot abundances versus time.

        Parameters
        ----------
        show : bool
            Show the plot. Will block until plot is closed.
        savefig : str
            Save copy of figure to path. File extension will determine
            file type.
        title : str
            Add a title to the figure.

        Returns
        -------
        figure : matplotlib figure
            Figure object of the plot
        axes : matplotlib axes
            Axes object of the plot.
        """
        figure, axes = plt.subplots()
        axes.step(self.times, self.abundances, where="post")
        axes.set_xlabel("Time")
        axes.set_ylabel("Abundances")
        axes.legend(self.reaction_system.species, loc='center right')
        if title:
            axes.set_title(title)
        if show:
            plt.show()
        if savefig:
            figure.savefig(savefig)
        return figure, axes


class DeterministicSimulator(ReactionSimulator):
    """Class for deterministic simulation
    Inherits from base class `ReactionSimulator`.
    Simulate species abundances deterministically.

    Parameters
    ----------
    reaction_system : `ReactionSystem` class object
        `ReactionSystem` instance with all the details regarding
        reaction to simulate.
    initial_abundances : list
        Initial abundances (concentrations) of each species.
    temperature : float
        Temperature when the reaction take place.
    t_span : tuple
        Time span to simulate: (begin, end)
    dt : float
        Simulation time step size.

    Methods
    -------
    simulate(method='bdf', epsilon = 1e-06)
    diff_func(t, y)
    """
    def __init__(self, reaction_system, initial_abundances, temperature,
                 t_span, dt=0.01):
        self.reaction_system = reaction_system
        self.temperature = temperature
        self.t_span = t_span
        self.times = [self.t_span[0]]
        self.abundances = [initial_abundances]
        try:
            dt = float(dt)
        except:
            raise ValueError("Could not convert dt=%s to float" % dt)
        self.dt = dt
        self._validate_arguments()
        self.ode_integrator = ode_solver.ODE_solver(
            self.diff_func, initial_abundances, t_span, self.dt, False)

    def simulate(self, method='bdf', epsilon=1e-06):
        """Simulate species abundances deterministically.

        We implemented three methods to solve the ordinary differential 
        equation. Backward differentiation formula and backward euler
        are good for stiff functions, and rk45 is accurate for non-stiff
        functions. BDF is the most suitable for solving ODE problems in
        chemical kinetics, so our default method is set as BDF.

        Parameters
        ----------
        method : str
            Solver method. Defaults to 'bdf' -- backward
            differentiation formula.
        epsilon : float
            Tolerance of error. Default is 1e-06

        Returns
        -------
        self.times : list
            Times.
        self.abundances : list
            Abundances of species at every time step

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
            Time.
        y : array_like
            Value of the variable we would like to integrate.

        Returns
        -------
        function
            right hand side of the ODE function
        '''
        if (y < 0).any():
            y = [0 if i < 0 else i for i in y]
        return self.reaction_system.calculate_reaction_rate(
            y, self.temperature)

    def plot_simulation(self, show=True, savefig=None, title=None):
        """Plot concentrations versus time.

        Parameters
        ----------
        show : bool
            Show the plot. Will block until plot is closed.
        savefig : str
            Save copy of figure to path. File extension will determine
            file type.
        title : str
            Add a title to the figure.

        Returns
        -------
        figure : matplotlib figure
            Figure object of the plot
        axes : matplotlib axes
            Axes object of the plot.
        """
        figure, axes = plt.subplots()
        axes.plot(self.times, self.abundances)
        axes.set_xlabel("Time")
        axes.set_ylabel("Concentrations")
        axes.legend(self.reaction_system.species, loc='center right')
        if title:
            axes.set_title(title)
        if show:
            plt.show()
        if savefig:
            figure.savefig(savefig)
        return figure, axes


if __name__ == "__main__":
    import doctest
    doctest.testmod()
