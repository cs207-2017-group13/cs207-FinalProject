from scipy.integrate import ode
from numpy.linalg import norm
import numpy as np

class ODE_solver:
    """ODE Solver class

    Implemented three ode solvers: backward euler method,
    Runge-Kutta-Fehlberg, and backward differentiation formula.

    Parameters
    ----------
    func :   function
             ordinary differential equation to be integrated
    y0 :     array of floats
             initial value of the variable y
    t_span : tuple of floats
             time span, first element indicating the initial time,
             second element indicating the stopping time
    dt :     float
             size of time step
    neg_y_allowed : bool
             default True
             indicating if negative y values are allowed to be passed
             to the function the user want to integrate

    Methods
    -------
    backward_euler()
    backward_euler_step()
    rk45()
    rk45_step()
    BDF()
    """
    def __init__(self, func, y0, t_span, dt, neg_y_allowed=True):
        self.diff_function = func
        self.t_span = t_span
        self.dt = dt
        self.t = [self.t_span[0]]
        self.neg_y_allowed =neg_y_allowed
        if (not self.neg_y_allowed):
            if (type(y0) is float) or (type(y0) is int):
                if y0<0:
                    self.y = [0.]
                else:
                    self.y = [y0]
            else:
                if (np.array(y0)<0).any():
                    self.y = [[0. if i < 0 else i for i in y0]]
                else:
                    self.y = [y0]
        else:
            self.y = [y0]
        

    def backward_euler(self, epsilon = 1e-06):
        """Solve the ODE using backward euler method.

        Parameters
        ----------
        epsilon : float
                 tolerance of error of the zero value

        Returns
        -------
        self.t :  array
                 all times evaluated
        self.y :  array
                 variable values at each time t

        Examples
        --------
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.5)
        >>> ode.backward_euler()
        ([1, 1.5, 2.0], [0.8, 2.3, 4.3])
        """
        while (self.t[-1]+self.dt) <= self.t_span[-1]:
            message = self.backward_euler_step(epsilon)
            if message == "Failure":
                print("Fixed point iteration does not converge.")
        if self.t[-1] < self.t_span[-1]:
            message = self.backward_euler_step(epsilon)
            if message == "Failure":
                print("Fixed point iteration does not converge.")
        return self.t, self.y

    def backward_euler_step(self, epsilon = 1e-06):
        """Solve the ODE one step/time forward 
           using backward euler method.
           

        Parameters
        ----------
        epsilon : float
                 tolerance of error of the zero value

        Returns
        -------
        message : "Success"
                 indicating if the integration is successful

        Examples
        --------
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.1)
        >>> ode.backward_euler_step()
        'Success'
        """
        if self.t[-1] + self.dt <= self.t_span[-1]:
            # fixed point iteration
            prev_y = self.y[-1]
            curr_y = self.y[-1]+ self.dt*np.array(self.diff_function(self.t[-1]+self.dt, prev_y))
            # count number of iterations to avoid looping endlessly
            j = 1
            while (norm(abs(curr_y-prev_y)) > epsilon) & (j<=1000):
                j+=1
                prev_y = curr_y
                curr_y = self.y[-1]+ self.dt*np.array(self.diff_function(self.t[-1]+self.dt, prev_y))
            if j > 1000:
                return "Failure"
            if (type(curr_y) is float) or (type(curr_y) is int):
                if (not self.neg_y_allowed) and (curr_y < 0):
                    curr_y = 0.
            else:
                if (not self.neg_y_allowed) and (curr_y < 0).any():
                    curr_y = [0. if i <0 else i for i in curr_y]
            self.y.append(curr_y)
            self.t.append(self.t[-1]+self.dt)
        elif self.t[-1] < self.t_span[-1]:
            prev_y = self.y[-1]
            curr_y = self.y[-1]+ self.dt*np.array(self.diff_function(self.t[-1]+self.dt, prev_y))
            j = 1
            while (norm(abs(curr_y-prev_y)) > epsilon) & (j<=1000):
                j += 1
                prev_y = curr_y
                curr_y = self.y[-1] + self.dt*np.array(self.diff_function(self.t[-1]+self.dt, prev_y))
            if j > 1000:
                return "Failure"
            if (type(curr_y) is float) or (type(curr_y) is int):
                if (not self.neg_y_allowed) and (curr_y < 0):
                    curr_y = 0.
            else:
                if (not self.neg_y_allowed) and (curr_y < 0).any():
                    curr_y = [0. if i <0 else i for i in curr_y]
            self.y.append(curr_y)
            self.t.append(self.t_span[-1])
        return "Success"
        
    def rk45(self, epsilon = 1e-06):
        """Solve the ODE using Runge-Kutta-Fehlberg method.

        Parameters
        ----------
        epsilon : float
                 tolerance of error of the zero value

        Returns
        -------
        self.t  : array
                 all times evaluated
        self.y  : array
                 variable values at each time t

        Examples
        --------
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.1)
        >>> ode.rk45()
        ([1, 1.1, 1.5, 2], [0.8, 1.01, 2.0500000000000003, 3.8000000000000003])
        """
        while (self.t[-1]+self.dt) <= self.t_span[-1]:
            message = self.rk45_step(epsilon)
        if self.t[-1] < self.t_span[-1]:
            message = self.rk45_step(epsilon)
        if message == "Failure":
            print("This ODE system should not be solved by RK45.",
                " Try backward euler with small step size.")
        return self.t, self.y

    def rk45_step(self, epsilon = 1e-06):
        """Solve the ODE one step/time forward 
            using Runge-Kutta-Fehlberg method.

        Parameters
        ----------
        epsilon : float
                 tolerance of error of the zero value

        Returns
        -------
        message : "Success"
                 indicating if the integration is successful

        Examples
        --------
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.1)
        >>> ode.rk45_step()
        'Success'
        """
        message = "Failure"
        # count number of iterations to avoid looping endlessly
        j=1
        while ((self.t[-1]+self.dt) <= self.t_span[-1]) & (j<1000):
            j += 1
            warning_print = False
            if (self.dt < 1e-10) & (warning_print == False):
                print("Warning: Step size too small.")
                warning_print = True
            k1 = self.dt*np.array(self.diff_function(self.t[-1], self.y[-1]))
            k2 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt/4, self.y[-1]+k1/2))
            k3 = self.dt*np.array(self.diff_function(self.t[-1]+3*self.dt/8, self.y[-1]+3*k1/32+9*k2/32))
            k4 = self.dt*np.array(self.diff_function(self.t[-1]+12*self.dt/13, 
                self.y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197))
            k5 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt, 
                self.y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104))
            k6 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt/2, 
                self.y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40))
            w1 = self.y[-1] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
            w2 = self.y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            if norm(abs(w2-w1))==0:
                if (type(w2) is float) or (type(w2) is int):
                    if (not self.neg_y_allowed) and (w2 < 0):
                        w2 = 0.
                else:
                    if (not self.neg_y_allowed) and (w2 < 0).any():
                        w2 = [0. if i <0 else i for i in w2]
                self.y.append(w2)
                self.t.append(self.t[-1]+self.dt)
                self.dt = 4*self.dt
                message = "Success"
                return message
            q = norm(0.84 * (epsilon*self.dt/abs(w2-w1))**0.25)
            # adjust step size
            if q>=1:
                if (type(w2) is float) or (type(w2) is int):
                    if (not self.neg_y_allowed) and (w2 < 0):
                        w2 = 0.
                else:
                    if (not self.neg_y_allowed) and (w2 < 0).any():
                        w2 = [0. if i <0 else i for i in w2]
                self.y.append(w2)
                self.t.append(self.t[-1]+self.dt)
                if q>= 4:
                    self.dt = 4*self.dt
                else:
                    self.dt = q*self.dt
                message = "Success"
                return message
            if q <= 0.1:
                self.dt = 0.1*self.dt
            else:
                self.dt = q*self.dt
        if ((self.t[-1]+self.dt) > self.t_span[-1]) & (j<1000):
            self.dt = self.t_span[-1]-self.t[-1]
            k1 = self.dt*np.array(self.diff_function(self.t[-1], self.y[-1]))
            k2 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt/4, self.y[-1]+k1/2))
            k3 = self.dt*np.array(self.diff_function(self.t[-1]+3*self.dt/8, self.y[-1]+3*k1/32+9*k2/32))
            k4 = self.dt*np.array(self.diff_function(self.t[-1]+12*self.dt/13, 
                self.y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197))
            k5 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt, 
                self.y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104))
            k6 = self.dt*np.array(self.diff_function(self.t[-1]+self.dt/2, 
                self.y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40))
            w2 = self.y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            if (type(w2) is float) or (type(w2) is int):
                if (not self.neg_y_allowed) and (w2 < 0):
                    w2 = 0.
            else:
                if (not self.neg_y_allowed) and (w2 < 0).any():
                    w2 = [0. if i <0 else i for i in w2]
            self.t.append(self.t_span[-1])
            self.y.append(w2)
            message = "Success"
        return message

    def BDF(self, epsilon = 1e-06):
        """Solve the ODE using backward differentiation formula.

        Parameters
        ----------
        epsilon : float
                 tolerance of error of the zero value

        Returns
        -------
        self.t :  array
                 all times evaluated
        self.y :  array
                 variable values at each time t

        Examples
        --------
        >>> func = lambda t, y: 2*t
        >>> obj = ODE_solver(func, 0.8, [1, 2], 0.5)
        >>> obj.BDF()
        ([1, 1.5, 2.0], [0.8, 2.0500030861821203, 3.8000031911685124])
        """
        r = ode(self.diff_function).set_integrator('vode', method='bdf',
            atol=epsilon, with_jacobian=False)
        r.set_initial_value(self.y[-1], self.t[-1])
        while r.successful() and r.t < self.t_span[1]:
            r.integrate(r.t + self.dt)
            if len(r.y)==1:
                y_val=r.y[0]
            else:
                y_val = r.y
            if (type(y_val) is float) or (type(y_val) is int):
                if (not self.neg_y_allowed) and (y_val < 0):
                    y_val = 0.
            else:
                if (not self.neg_y_allowed) and (y_val < 0).any():
                    y_val = [0. if i <0 else i for i in y_val]
            self.y.append(y_val)
            self.t.append(r.t)
        return self.t, self.y

if __name__ == "__main__":
    import doctest
    doctest.testmod()

