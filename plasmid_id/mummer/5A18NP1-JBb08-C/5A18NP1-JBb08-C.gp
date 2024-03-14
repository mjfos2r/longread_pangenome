set terminal png tiny size 800,800
set output "mummer_output/5A18NP1-JBb08-C/5A18NP1-JBb08-C.png"
set xtics rotate ( \
 "tig00000024" 1.0, \
 "tig00000022" 17742.0, \
 "tig00000001" 50941.0, \
 "tig00000005" 1002200.0, \
 "tig00000026" 1054992.0, \
 "tig00000004" 1091915.0, \
 "tig00000003" 1171757.0, \
 "tig00000007" 1216125.0, \
 "tig00000025" 1264303.0, \
 "tig00000002" 1282194.0, \
 "tig00000016" 1300718.0, \
 "tig00000014" 1332749.0, \
 "tig00000012" 1381620.0, \
 "tig00000017" 1420559.0, \
 "tig00000020" 1489372.0, \
 "tig00000009" 1509028.0, \
 "tig00000019" 1565418.0, \
 "tig00000021" 1584170.0, \
 "tig00000008" 1605273.0, \
 "tig00000015" 1695282.0, \
 "tig00000013" 1776736.0, \
 "tig00000010" 1806860.0, \
 "tig00000018" 1847308.0, \
 "tig00000023" 1891780.0, \
 "" 1913082 \
)
set ytics ( \
 "tig00000024" 1.0, \
 "tig00000022" 17742.0, \
 "tig00000001" 50941.0, \
 "tig00000005" 1002200.0, \
 "tig00000026" 1054992.0, \
 "tig00000004" 1091915.0, \
 "tig00000003" 1171757.0, \
 "tig00000007" 1216125.0, \
 "tig00000025" 1264303.0, \
 "tig00000002" 1282194.0, \
 "tig00000016" 1300718.0, \
 "tig00000014" 1332749.0, \
 "tig00000012" 1381620.0, \
 "tig00000017" 1420559.0, \
 "tig00000020" 1489372.0, \
 "tig00000009" 1509028.0, \
 "tig00000019" 1565418.0, \
 "tig00000021" 1584170.0, \
 "tig00000008" 1605273.0, \
 "tig00000015" 1695282.0, \
 "tig00000013" 1776736.0, \
 "tig00000010" 1806860.0, \
 "tig00000018" 1847308.0, \
 "tig00000023" 1891780.0, \
 "" 1913082 \
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
set xrange [1:1913082]
set yrange [1:1913082]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/5A18NP1-JBb08-C/5A18NP1-JBb08-C.fplot" title "FWD" w lp ls 1, \
 "mummer_output/5A18NP1-JBb08-C/5A18NP1-JBb08-C.rplot" title "REV" w lp ls 2
