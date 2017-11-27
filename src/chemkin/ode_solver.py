import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

class ODE_solver:
    def __init__(self, func, y0, t_final, dt):
        # what should be initialized?

    def backward_euler(self, stepsize, nsteps, initial_t, initial_y, diff_function):
        # initial_t: int, initial_y: array of concentration, 
        # diff_function: function to calculate reaction rate
        y = np.zeros(nsteps+1, len(initial_y))
        y[0] = initial_y
        for i in range(nsteps):
            y_prev = y[i]
            y_curr = y[i] + stepsize*diff_function(initial_t+(i+1)*stepsize, y_prev)
            j = 0
            # fixed point iteration
            while abs(y_curr-y_prev) > 1e-6:
                y_prev = y_curr
                y_curr = y[i] + stepsize*diff_function(initial_t+(i+1)*stepsize, y_prev)
                j += 1
                if j>10000:
                    raise RuntimeError("The sequence does not converge.")
                # assign the approximation to result
            y[i+1] = y_curr
        return y.transpose()

    def rk45(self):
        # runge kutta

    def plot(self, stepsize, nsteps, initial_t, y, species):
        # plot y's 
        t = np.zeros(nsteps+1)
        for i in range(nsteps+1):
            t[i] = initial_t + i*stepsize
        for i in range(len(species)):
            plt.plot(t, y[i], label=species[i])
            plt.xlabel("Time")
            plt.ylabel("Concentration")
            plt.legend(loc='best')
            plt.show()
            # maybe save figure?

