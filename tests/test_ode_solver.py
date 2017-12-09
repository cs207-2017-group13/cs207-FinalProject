
import chemkin.ode_solver as ode_solver
import numpy as np

def test_backward_euler():
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

def test_backward_euler_negy():
	func = lambda t, y: 2*y+t
	obj = ode_solver.ODE_solver(func, -1, [10, 10.02], 0.01, False)
	t, y = obj.backward_euler()
	assert np.allclose(t, [10, 10.01, 10.02])
	assert np.allclose(y, [0.0, 0.10214284080000001, 0.20647226983790093])

def test_BDF_negy():
	func = lambda t, y: 2*y+t
	obj = ode_solver.ODE_solver(func, -1, [10, 10.02], 0.01, False)
	t, y = obj.BDF()
	assert np.allclose(t, [10, 10.01, 10.02])
	assert np.allclose(y, [0.0, 0.10105761562779661, 0.20425719527441025])

def test_backward_euler_step_failure():
	func = lambda t, y: 2*y+t
	obj = ode_solver.ODE_solver(func, 1., [10., 11.], 1)
	assert obj.backward_euler_step() == "Failure"

def test_backward_euler_negylist():
	func = lambda t, y: 2*np.array(y)+np.array([t, t])
	obj = ode_solver.ODE_solver(func, [-1, 1], [10., 10.1], 0.1, False)
	t, y = obj.backward_euler()
	y0 = y[0]
	y1 = y[1]
	assert np.allclose(t, [10, 10.1])
	assert np.allclose(y0, [0.0, 1])
	assert np.allclose(y1, [1.26249987, 2.51249985])

def test_BDF_negylist():
	func = lambda t, y: 2*np.array(y)+np.array([t, t])
	obj = ode_solver.ODE_solver(func, [-1, 1], [10., 10.1], 0.1, False)
	t, y = obj.BDF()
	y0 = y[0]
	y1 = y[1]
	assert np.allclose(t, [10, 10.1])
	assert np.allclose(y0, [0.0, 1])
	assert np.allclose(y1, [1.1123736, 2.33377809])
