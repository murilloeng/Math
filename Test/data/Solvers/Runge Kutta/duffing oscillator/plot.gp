set grid
set key above
set xlabel "t"
set ylabel "x(t)"

plot\
	"Test/data/Solvers/Runge Kutta/duffing oscillator/numeric.txt" using 4:1 with lines linecolor rgb "#0000ff" title "Numeric",\
	"Test/data/Solvers/Runge Kutta/duffing oscillator/reference.txt" using 1:2 with points pointtype 7 pointsize 0.5 linecolor rgb "#ff0000" title "Reference"