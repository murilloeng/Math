set grid
set key above
set xlabel "t (s)"
set ylabel "{/Symbol q} (rad)"

plot\
	"numeric.txt" using 4:1 with lines linecolor rgb "#0000ff" title "Numeric",\
	"reference.txt" using 1:2 every 10 with points pointtype 7 linecolor rgb "#ff0000" title "Reference"