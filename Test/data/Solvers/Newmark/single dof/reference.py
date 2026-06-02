import scipy.integrate
from numpy import pi, cos, sin, sqrt, linspace

#data
np = 40
ns = 1001
k = 1.00e+00
c = 5.00e-02
m = 1.00e+00
x0 = 1.00e+00
v0 = 0.00e+00
Fr = 1.00e+00
wf = 2.00e+00

w0 = sqrt(k / m)
tf = 2 * pi * np / w0

#system
def pendulum_ode(t, y):
	x, v = y
	dx_dt = v
	dv_dt = Fr * cos(wf * t) - k * x - c * v
	return [dx_dt, dv_dt]

#solve
t_eval = linspace(0, tf, ns)
sol = scipy.integrate.solve_ivp(pendulum_ode, [0, tf], [x0, v0], t_eval = t_eval, method='DOP853')

#write
with open("reference.txt", "w") as file:
	for i in range(len(sol.t)):
		file.write("%+.6e %+.6e \n" % (sol.t[i], sol.y[0][i]))