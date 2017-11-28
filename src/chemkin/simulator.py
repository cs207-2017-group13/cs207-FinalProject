import chemkin.ode_solver as ode_solver
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

    def simulate(self, method='rk45', epsilon = 1e-06):
        choices = ['backward_euler','rk45']
        if method not in choices:
            raise ValueError("Wrong method.")
        if method == 'rk45':
            self.time, self.concentrations = self.ode_solver.rk45(epsilon)
        else:
            self.time, self.concentrations = self.ode_solver.backward_euler(epsilon)

    def simulate_step(self, method='rk45', epsilon = 1e-06):
        """Calculate concentrations between last time and `t_final`.
        Neccessary? refer to scipy.integrate.solve_ivp
        """
        # initialize solver
        if self.time[-1] < self.t_span[-1]:
            choices = ['backward_euler','rk45']
            if method not in choices:
                raise ValueError("Wrong method.")
            if method == 'rk45':
                message = self.ode_solver.rk45_step(epsilon)
                self.dt, self.time, self.concentrations = self.ode_solver.dt, self.ode_solver.t, self.ode_solver.y
            else:
                message = self.ode_solver.backward_euler_step(epsilon)
                self.dt, self.time, self.concentrations = self.ode_solver.dt, self.ode_solver.t, self.ode_solver.y
        else:
            raise IndexError("Time exceeds time span.")
        return message

    def diff_func(self, t, y):
        self.reaction_system.calculate_reaction_rate(y, self.temperature)

        
