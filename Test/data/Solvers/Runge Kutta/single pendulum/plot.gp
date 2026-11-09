set grid
set key above
set xlabel "t (s)"
set ylabel "{/Symbol q} (rad)"

plot\
	"Test/data/Solvers/Runge Kutta/single pendulum/numeric.txt" using 4:1 with lines linecolor rgb "#0000ff" title "Numeric",\
	"Test/data/Solvers/Runge Kutta/single pendulum/reference.txt" using 1:2 every 10 with points pointtype 7 linecolor rgb "#ff0000" title "Reference"