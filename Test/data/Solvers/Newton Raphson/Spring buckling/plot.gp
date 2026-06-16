K = 1.00e+00
F = 1.00e+00
L = 1.00e+00
q = 1.00e-04

set grid
set key above
set samples 20
set xlabel "{/Symbol q}"
set ylabel "{/Symbol l}"

plot\
	"numeric.txt" using 1:2 with lines linecolor rgb "#0000ff" title "Numeric",\
	K / (F * L) * x / sin(x + q) with points pointtype 7 linecolor rgb "#ff0000" title "Reference"