set grid
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

plot\
	'numeric.txt' using ($4) : (sqrt($2**2 + $3**2)) with lines title "Harmonic 1"