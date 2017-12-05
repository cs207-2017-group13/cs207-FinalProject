from scipy.optimize import newton
from scipy.integrate import ode

class ODE_solver:
    """ODE Solver class

    Implemented three ode solvers: backward euler method,
    Runge-Kutta-Fehlberg, and backward differentiation formula.

    Parameters:
    ===========
    func:   ordinary differential equation to be integrated
    y0:     initial value of the variable y
    t_span: time span, first element indicating the initial time,
            second element indicating the stopping time
    dt:     size of time step

    Methods:
    ========
    backward_euler()
    backward_euler_step()
    rk45()
    rk45_step()
    BDF()
    """
    def __init__(self, func, y0, t_span, dt):
        self.diff_function = func
        self.t_span = t_span
        self.dt = dt
        self.t = [self.t_span[0]]
        self.y = [y0]

    def backward_euler(self, epsilon = 1e-06):
        """Solve the ODE using backward euler method.

        INPUTS:
        =======
        epsilon: float
                 tolerance of error of the zero value

        RETURNS:
        ========
        self.t:  array
                 all times evaluated
        self.y:  array
                 variable values at each time t

        EXAMPLES:
        =========
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.5)
        >>> ode.backward_euler()
        ([1, 1.5, 2.0], [0.8, 2.3000000000000003, 4.300000000000001])
        """
        while (self.t[-1]+self.dt) <= self.t_span[-1]:
            self.backward_euler_step(epsilon)
        if self.t[-1] < self.t_span[-1]:
            self.backward_euler_step(epsilon)
        return self.t, self.y

    def backward_euler_step(self, epsilon = 1e-06):
        """Solve the ODE one step/time forward 
           using backward euler method.
           

        INPUTS:
        =======
        epsilon: float
                 tolerance of error of the zero value

        RETURNS:
        ========
        message: "Success"
                 indicating if the integration is successful

        EXAMPLES:
        =========
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.1)
        >>> ode.backward_euler_step()
        'Success'
        """
        if self.t[-1] + self.dt <= self.t_span[-1]:
            # newton raphson
            y_curr = newton(lambda update_y: update_y - self.y[-1] - 
                self.dt*self.diff_function(self.t[-1]+self.dt, update_y), self.y[-1],
                tol = epsilon)
            # assign the approximation to result
            self.y.append(y_curr)
            self.t.append(self.t[-1]+self.dt)
        elif self.t[-1] < self.t_span[-1]:
            y_curr = newton(lambda update_y: update_y - self.y[-1] - 
                self.dt*self.diff_function(self.t_span[-1], update_y), self.y[-1],
                tol = epsilon)
            # assign the approximation to result
            self.y.append(y_curr)
            self.t.append(self.t_span[-1])
        return "Success"
        
    def rk45(self, epsilon = 1e-06):
        """Solve the ODE using Runge-Kutta-Fehlberg method.

        INPUTS:
        =======
        epsilon: float
                 tolerance of error of the zero value

        RETURNS:
        ========
        self.t:  array
                 all times evaluated
        self.y:  array
                 variable values at each time t

        EXAMPLES:
        =========
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
                " Try backward euler.")
        return self.t, self.y

    def rk45_step(self, epsilon = 1e-06):
        """Solve the ODE one step/time forward 
            using Runge-Kutta-Fehlberg method.

        INPUTS:
        =======
        epsilon: float
                 tolerance of error of the zero value

        RETURNS:
        ========
        message: "Success"
                 indicating if the integration is successful

        EXAMPLES:
        =========
        >>> func = lambda t, y: 2*t
        >>> ode = ODE_solver(func, 0.8, [1, 2], 0.1)
        >>> ode.rk45_step()
        'Success'
        """
        message = "Failure"
        j=1
        while ((self.t[-1]+self.dt) <= self.t_span[-1]) & (j<1000):
            j += 1
            if self.dt < 1e-10:
            	print("Warning: Step size too small.")
            k1 = self.dt*self.diff_function(self.t[-1], self.y[-1])
            k2 = self.dt*self.diff_function(self.t[-1]+self.dt/4, self.y[-1]+k1/2)
            k3 = self.dt*self.diff_function(self.t[-1]+3*self.dt/8, self.y[-1]+3*k1/32+9*k2/32)
            k4 = self.dt*self.diff_function(self.t[-1]+12*self.dt/13, 
                self.y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197)
            k5 = self.dt*self.diff_function(self.t[-1]+self.dt, 
                self.y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104)
            k6 = self.dt*self.diff_function(self.t[-1]+self.dt/2, 
                self.y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40)
            w1 = self.y[-1] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
            w2 = self.y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            if abs(w2-w1)==0:
                self.y.append(w2)
                self.t.append(self.t[-1]+self.dt)
                self.dt = 4*self.dt
                message = "Success"
                return message
            q = 0.84 * (epsilon*self.dt/abs(w2-w1))**0.25
            if q>=1:
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
            k1 = self.dt*self.diff_function(self.t[-1], self.y[-1])
            k2 = self.dt*self.diff_function(self.t[-1]+self.dt/4, self.y[-1]+k1/2)
            k3 = self.dt*self.diff_function(self.t[-1]+3*self.dt/8, self.y[-1]+3*k1/32+9*k2/32)
            k4 = self.dt*self.diff_function(self.t[-1]+12*self.dt/13, 
                self.y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197)
            k5 = self.dt*self.diff_function(self.t[-1]+self.dt, 
                self.y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104)
            k6 = self.dt*self.diff_function(self.t[-1]+self.dt/2, 
                self.y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40)
            w2 = self.y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            self.t.append(self.t_span[-1])
            self.y.append(w2)
            message = "Success"
        return message

    def BDF(self, epsilon = 1e-06):
        """Solve the ODE using backward differentiation formula.

        INPUTS:
        =======
        epsilon: float
                 tolerance of error of the zero value

        RETURNS:
        ========
        self.t:  array
                 all times evaluated
        self.y:  array
                 variable values at each time t

        EXAMPLES:
        =========
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
            self.y.append(y_val)
            self.t.append(r.t)
        return self.t, self.y

func = lambda t, y: 2*t
obj = ODE_solver(func, 0.8, [1, 2], 0.5)
print(obj.BDF())
