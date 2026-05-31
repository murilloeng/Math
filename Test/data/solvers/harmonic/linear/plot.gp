set grid
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

plot\
	'Test/data/solvers/harmonic/linear/numeric.dat' using ($4) : (sqrt($2**2 + $3**2)) with lines title "Harmonic 1"