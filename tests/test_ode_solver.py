
import chemkin.ode_solver as ode_solver
import numpy as np

def test_back_euler():
	func = lambda t, y: 2*t
	ode = ode_solver.ODE_solver(func, 0.8, [1, 2], 0.2)
	t, y = ode.backward_euler()
	assert np.allclose(t, [1, 1.2, 1.4, 1.5999999999999999, 
		1.7999999999999998, 1.9999999999999998, 2])
	assert np.allclose(y, [0.8, 1.28, 1.8399999999999999, 2.48, 3.2, 4.0, 4.88])


def test_rk45():
	func = lambda t, y: 2*t
	ode = ode_solver.ODE_solver(func, 0.8, [1, 2], 0.2)
	t, y = ode.rk45()
	assert np.allclose(t, [1, 1.2, 2.0])
	assert np.allclose(y, [0.8, 1.2400000000000002, 3.8000000000000007])