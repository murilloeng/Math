set grid
set key above
set parametric
set samples 120
set trange [0.001 : 1]
set xrange [0.6 : 1.5]
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

l = 0.05
m = 1.0
c = 0.05
k = 1.0
knl = 3.0

K(A) = k + 0.75 * knl * A**2
q1(A) = c**2 - 2 * K(A) * m
q2(A) = K(A)**2 - l**2 / A**2
d(A) = q1(A)**2 - 4 * m**2 * q2(A)
w2(A, s) = d(A) < 0 ? 1/0 : (-q1(A) + s * sqrt(d(A))) / (2 * m**2)
w(A, s) = w2(A, s) < 0 ? 1/0 : sqrt(w2(A, s))

plot\
	'Test/data/Solvers/Harmonic/Duffing/numeric.txt' using ($4) : (sqrt($2**2 + $3**2)) with lines linecolor rgb "#0000ff" title "Numeric",\
	w(t, -1), t with points pointtype 7 pointsize 0.5 linecolor rgb "#ff0000" title "Exact",\
	w(t, +1), t with points pointtype 7 pointsize 0.5 linecolor rgb "#ff0000" notitle