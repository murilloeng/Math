set grid
set key above
set xlabel "t (s)"
set ylabel "{/Symbol q} (rad)"

plot\
	"Test/data/Solvers/Newmark/single dof/numeric.txt" using 4:1 with lines linecolor rgb "#0000ff" title "Numeric",\
	"Test/data/Solvers/Newmark/single dof/reference.txt" using 1:2 with points pointtype 7 pointsize 0.5 linecolor rgb "#ff0000" title "Reference"