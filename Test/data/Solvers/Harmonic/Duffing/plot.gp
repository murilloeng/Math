set grid
set key above
set ylabel "Amplitude"
set xlabel "{/Symbol w}"

plot\
	'numeric.txt' using ($8) : (sqrt($6**2 + $7**2)) with lines title "Harmonic 3"