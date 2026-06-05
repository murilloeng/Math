from scipy.integrate import solve_ivp
from numpy import sin, cos, pi, linspace

#data
ns = 2001
q1 = +pi / 4
q2 = -pi / 4
g = 9.81e+00
T = 1.00e+01
L1 = 1.00e+00
L2 = 1.00e+00
m1 = 1.00e+00
m2 = 1.00e+00
y0 = [q1, q2, 0, 0]

#system
def double_pendulum(t, y):
	#data
	q1 = y[0]
	q2 = y[1]
	w1 = y[2]
	w2 = y[3]
	dq = q2 - q1
	mt = m1 + m2
	d1 = L1 * (mt - m2 * cos(dq)**2)
	d2 = L2 * (mt - m2 * cos(dq)**2)
	#system
	a1 = (+m2 * L1 * w1**2 * sin(dq) * cos(dq) + g * m2 * sin(q2) * cos(dq) + m2 * L2 * w2**2 * sin(dq) - g * mt * sin(q1)) / d1
	a2 = (-m2 * L2 * w2**2 * sin(dq) * cos(dq) + g * mt * sin(q1) * cos(dq) - mt * L1 * w1**2 * sin(dq) - g * mt * sin(q2)) / d2
	#return
	return [w1, w2, a1, a2]

#solve
t_eval = linspace(0, T, ns)
sol = solve_ivp(double_pendulum, [0, T], y0, t_eval = t_eval, method = 'DOP853')

#write
with open("reference.txt", "w") as file:
	for i in range(len(sol.t)):
		t = sol.t[i]
		q1 = sol.y[0][i]
		q2 = sol.y[1][i]
		file.write("%+.6e %+.6e %+.6e \n" % (t, q1, q2))