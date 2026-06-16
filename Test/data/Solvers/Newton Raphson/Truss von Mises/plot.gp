set grid
set key above
set samples 20
set xlabel "z"
set ylabel "{/Symbol l}"

plot\
	"numeric.txt" using 1:2 with lines linecolor rgb "#0000ff" title "Numeric",\
	x * (1 - x**2) with points pointtype 7 linecolor rgb "#ff0000" title "Reference"