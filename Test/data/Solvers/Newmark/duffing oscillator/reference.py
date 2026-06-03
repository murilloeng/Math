import scipy.integrate
from numpy import pi, cos, sin, sqrt, linspace

#data
np = 10
ns = 1001
m = 1.00e+00
c = 5.00e-02
k = 1.00e+00
F = 5.00e-01
w = 2.00e+00
v0 = 0.00e+00
x0 = 0.00e+00
knl = 3.00e+00
T = 2 * pi * np / sqrt(k / m)

#system
def single_pendulum(t, y):
	x, v = y
	dx_dt = v
	dv_dt = F * sin(w * t) - k * x - c * v - knl * x * x * x
	return [dx_dt, dv_dt]

#solve
t_eval = linspace(0, T, ns)
sol = scipy.integrate.solve_ivp(single_pendulum, [0, T], [x0, v0], t_eval = t_eval, method='DOP853')

#write
with open("reference.txt", "w") as file:
	for i in range(len(sol.t)):
		file.write("%+.6e %+.6e \n" % (sol.t[i], sol.y[0][i]))