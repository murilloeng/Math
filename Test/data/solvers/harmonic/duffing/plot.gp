set grid
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

plot\
	'Test/data/solvers/harmonic/duffing/numeric.dat' using ($8) : (sqrt($6**2 + $7**2)) with linespoints title "Harmonic 3"