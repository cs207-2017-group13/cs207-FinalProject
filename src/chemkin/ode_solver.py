from scipy.optimize import newton
from scipy.integrate import ode

class ODE_solver:
    def __init__(self, func, y0, t_span, dt):
        self.diff_function = func
        self.t_span = t_span
        self.dt = dt
        self.t = [self.t_span[0]]
        self.y = [y0]

    def backward_euler(self, epsilon = 1e-06):
        while (self.t[-1]+self.dt) <= self.t_span[-1]:
            self.backward_euler_step(epsilon)
        return self.t, self.y

    def backward_euler_step(self, epsilon = 1e-06):
        if self.t[-1] + self.dt <= self.t_span[-1]:
            # newton raphson
            y_curr = newton(lambda update_y: update_y - self.y[-1] - 
                self.dt*self.diff_function(self.t[-1]+self.dt, update_y), self.y[-1])
            # assign the approximation to result
            self.y.append(y_curr)
            self.t.append(self.t[-1]+self.dt)
        return "Success"
        
    def rk45(self, epsilon = 1e-06):
        while (self.t[-1]+self.dt) <= self.t_span[-1]:
            message = self.rk45_step(epsilon)
        if self.t[-1] < self.t_span[-1]:
            message = self.rk45_step(epsilon)
        if message == "Failure":
            print("This ODE system should not be solved by RK45.",
                " Try backward euler.")
        return self.t, self.y

    def rk45_step(self, epsilon = 1e-06):
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
