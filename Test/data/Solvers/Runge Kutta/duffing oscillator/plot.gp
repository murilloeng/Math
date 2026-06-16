set grid
set key above
set xlabel "t"
set ylabel "x(t)"

plot\
	"numeric.txt" using 4:1 with lines linecolor rgb "#0000ff" title "Numeric",\
	"reference.txt" using 1:2 with points pointtype 7 pointsize 0.5 linecolor rgb "#ff0000" title "Reference"