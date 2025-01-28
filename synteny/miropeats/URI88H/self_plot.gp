set terminal png tiny size 800,800
set output "self_plot.png"
set size 1,1
set grid
unset key
set border 15
set tics scale 0
set xlabel "contig000014,"
set ylabel "contig000014,"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:29970]
set yrange [1:29970]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "self_plot.fplot" title "FWD" w lp ls 1, \
 "self_plot.rplot" title "REV" w lp ls 2
