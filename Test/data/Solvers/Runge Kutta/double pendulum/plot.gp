set grid
set key above
set xlabel "t (s)"
set ylabel "{/Symbol q} (rad)"

plot\
	"numeric.txt" using 7:1 with lines linecolor rgb "#0000ff" title "Numeric {/Symbol q}_1",\
	"numeric.txt" using 7:4 with lines linecolor rgb "#ff0000" title "Numeric {/Symbol q}_2",\
	"reference.txt" using 1:2 every 10 with points pointtype 6 pointsize 0.5 linecolor rgb "#0000ff" title "Reference {/Symbol q}_1",\
	"reference.txt" using 1:3 every 10 with points pointtype 6 pointsize 0.5 linecolor rgb "#ff0000" title "Reference {/Symbol q}_2"