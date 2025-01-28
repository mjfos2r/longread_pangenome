set terminal png tiny size 800,800
set output "mummer_output/BG001_C7/BG001_C7.png"
set xtics rotate ( \
 "tig00000002" 1.0, \
 "tig00000003" 69108.0, \
 "tig00000001" 123968.0, \
 "tig00000004" 1053398.0, \
 "tig00000005" 1101882.0, \
 "" 1117848 \
)
set ytics ( \
 "tig00000002" 1.0, \
 "tig00000003" 69108.0, \
 "tig00000001" 123968.0, \
 "tig00000004" 1053398.0, \
 "tig00000005" 1101882.0, \
 "" 1117848 \
)
set size 1,1
set grid
unset key
set border 0
set tics scale 0
set xlabel "REF"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
set xrange [1:1117848]
set yrange [1:1117848]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/BG001_C7/BG001_C7.fplot" title "FWD" w lp ls 1, \
 "mummer_output/BG001_C7/BG001_C7.rplot" title "REV" w lp ls 2
