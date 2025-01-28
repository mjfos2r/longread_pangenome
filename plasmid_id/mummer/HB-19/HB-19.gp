set terminal png tiny size 800,800
set output "mummer_output/HB-19/HB-19.png"
set xtics rotate ( \
 "tig00000010" 1.0, \
 "tig00000009" 91761.0, \
 "tig00000002" 1008990.0, \
 "tig00000004" 1058888.0, \
 "tig00000006" 1110507.0, \
 "tig00000003" 1161398.0, \
 "tig00000007" 1233889.0, \
 "tig00000005" 1357638.0, \
 "tig00000008" 1402305.0, \
 "" 1516310 \
)
set ytics ( \
 "tig00000010" 1.0, \
 "tig00000009" 91761.0, \
 "tig00000002" 1008990.0, \
 "tig00000004" 1058888.0, \
 "tig00000006" 1110507.0, \
 "tig00000003" 1161398.0, \
 "tig00000007" 1233889.0, \
 "tig00000005" 1357638.0, \
 "tig00000008" 1402305.0, \
 "" 1516310 \
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
set xrange [1:1516310]
set yrange [1:1516310]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/HB-19/HB-19.fplot" title "FWD" w lp ls 1, \
 "mummer_output/HB-19/HB-19.rplot" title "REV" w lp ls 2
