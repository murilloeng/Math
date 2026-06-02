import scipy.integrate
from numpy import pi, cos, sin, linspace

#data
ns = 1001
g = 9.81e+00
l = 1.00e+00
T = 1.00e+01
v0 = 0.00e+00
q0 = pi / 2

#system
def pendulum_ode(t, y):
	theta, omega = y
	dtheta_dt = omega
	domega_dt = -g / l * sin(theta)
	return [dtheta_dt, domega_dt]

#solve
t_eval = linspace(0, T, ns)
sol = scipy.integrate.solve_ivp(pendulum_ode, [0, T], [q0, v0], t_eval = t_eval, method='DOP853')

#write
with open("reference.txt", "w") as file:
	for i in range(len(sol.t)):
		file.write("%+.6e %+.6e \n" % (sol.t[i], sol.y[0][i]))