
import chemkin.ode_solver as ode_solver

def test_back_euler():
	 func = lambda t, y: 2*t
	 ode = ODE_solver(func, 0.8, [1, 2], 0.2)
	 assert ode.backward_euler()==([1, 1.2, 1.4, 1.5999999999999999, 1.7999999999999998, /
	 	1.9999999999999998, 2], [0.8, 1.28, 1.8399999999999999, 2.48, 3.2, 4.0, 4.8])
