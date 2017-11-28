import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

class ODE_solver:
    def __init__(self, func, y0, t_final, dt):
        self.diff_function = func
        self.y0 = y0
        self.t_final = t_final
        self.dt = dt
        self.t = [0.]
        self.y = [y0]
        self.h = dt

    def backward_euler(self, epsilon = 1e-06):
        nsteps = t_final//dt
        y = np.zeros(nsteps+1, len(self.y0))
        y[0] = self.y0
        for i in range(nsteps):
            y_prev = y[i]
            y_curr = y[i] + self.dt*self.diff_function((i+1)*self.dt, y_prev)
            j = 0
            # fixed point iteration
            while abs(y_curr-y_prev) > epsilon:
                y_prev = y_curr
                y_curr = y[i] + self.dt*self.diff_function((i+1)*self.dt, y_prev)
                j += 1
                if j>10000:
                    raise RuntimeError("The sequence does not converge.")
            # assign the approximation to result
            y[i+1] = y_curr
            self.t.append((i+1)*self.dt)
        return y.transpose()

    def backward_euler_step(self, epsilon = 1e-06):
        if self.t[-1] + self.dt <= self.t_final:
            y_prev = self.y[-1]
            y_curr = self.y[-1] + self.dt*self.diff_function(self.t[-1]+self.dt, y_prev)
            j = 0
            # fixed point iteration
            while abs(y_curr-y_prev) > epsilon:
                y_prev = y_curr
                y_curr = self.y[-1] + self.dt*self.diff_function(self.t[-1]+self.dt, y_prev)
                j += 1
                if j>10000:
                    raise RuntimeError("The sequence does not converge.")
            # assign the approximation to result
            self.y.append(y_curr)
            self.t.append(self.t[-1]+self.dt)
        return "Success"

    def rk45(self, epsilon = 1e-06):
        h = self.dt
        y = [self.y0]
        while (self.t[-1]+h) <= self.t_final:
            k1 = h*self.diff_function(self.t[-1], y[-1])
            k2 = h*self.diff_function(self.t[-1]+h/4, y[-1]+k1/2)
            k3 = h*self.diff_function(self.t[-1]+3*h/8, y[-1]+3*k1/32+9*k2/32)
            k4 = h*self.diff_function(self.t[-1]+12*h/13, 
                y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197)
            k5 = h*self.diff_function(self.t[-1]+h, 
                y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104)
            k6 = h*self.diff_function(self.t[-1]+h/2, 
                y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40)
            w1 = y[-1] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
            w2 = y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            q = 0.84 * (epsilon*h/abs(w2-w1))**0.25
            if q>=1:
                y.append(w2)
                self.t.append(self.t[-1]+h)
            if q <= 0.1:
                h = 0.1*h
            elif q>= 4:
                h = 4*h
            else:
                h = q*h
        y = np.array(y)
        return y.transpose()

    def rk45_step(self, epsilon = 1e-06):
        while (self.t[-1]+self.h) <= self.t_final:
            k1 = self.h*self.diff_function(self.t[-1], self.y[-1])
            k2 = self.h*self.diff_function(self.t[-1]+self.h/4, self.y[-1]+k1/2)
            k3 = self.h*self.diff_function(self.t[-1]+3*self.h/8, self.y[-1]+3*k1/32+9*k2/32)
            k4 = self.h*self.diff_function(self.t[-1]+12*self.h/13, 
                self.y[-1]+1932*k1/2197-7200*k2/2197+7296*k3/2197)
            k5 = self.h*self.diff_function(self.t[-1]+self.h, 
                self.y[-1]+439*k1/216-8*k2+3680*k3/513-845*k4/4104)
            k6 = self.h*self.diff_function(self.t[-1]+self.h/2, 
                self.y[-1]-8*k1/27+2*k2-3544*k3/2565+1859*k4*4104-11*k5/40)
            w1 = self.y[-1] + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
            w2 = self.y[-1] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55
            q = 0.84 * (epsilon*self.h/abs(w2-w1))**0.25
            if q>=1:
                self.y.append(w2)
                self.t.append(self.t[-1]+self.h)
                if q>= 4:
                    self.h = 4*self.h
                else:
                    self.h = q*self.h
                return "Success"
            if q <= 0.1:
                self.h = 0.1*self.h
            else:
                self.h = q*self.h
            












    

