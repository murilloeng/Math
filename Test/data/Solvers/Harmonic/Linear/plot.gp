set grid
set key above
set samples 100
set xrange [0.1 : 4]
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

k = 4.0
c = 0.1
m = 1.0
p = 1.0
A(w) = p / sqrt((k - m * w**2)**2 + (c * w)**2)

plot\
	'Test/data/Solvers/Harmonic/Linear/numeric.txt' using ($4) : (sqrt($2**2 + $3**2)) with lines linecolor rgb "#0000ff" title "Numeric",\
	A(x) with points pointtype 7 linecolor rgb "#ff0000" title "Exact"