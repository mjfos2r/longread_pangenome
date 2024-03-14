set terminal png tiny size 800,800
set output "mummer_output/5A18NP1-JBb08-B/5A18NP1-JBb08-B.png"
set xtics rotate ( \
 "tig00000008" 1.0, \
 "tig00000024" 45296.0, \
 "tig00000026" 63174.0, \
 "tig00000006" 103970.0, \
 "tig00000009" 142672.0, \
 "tig00000025" 185768.0, \
 "tig00000019" 205473.0, \
 "tig00000016" 1122257.0, \
 "tig00000012" 1165414.0, \
 "tig00000013" 1242949.0, \
 "tig00000027" 1299048.0, \
 "tig00000007" 1317397.0, \
 "tig00000018" 1359916.0, \
 "tig00000028" 1370544.0, \
 "tig00000014" 1431557.0, \
 "tig00000023" 1479673.0, \
 "tig00000022" 1507925.0, \
 "tig00000015" 1546991.0, \
 "tig00000021" 1574840.0, \
 "tig00000020" 1645355.0, \
 "tig00000010" 1700762.0, \
 "tig00000017" 1727200.0, \
 "tig00000003" 1737983.0, \
 "tig00000011" 1830640.0, \
 "" 1847903 \
)
set ytics ( \
 "tig00000008" 1.0, \
 "tig00000024" 45296.0, \
 "tig00000026" 63174.0, \
 "tig00000006" 103970.0, \
 "tig00000009" 142672.0, \
 "tig00000025" 185768.0, \
 "tig00000019" 205473.0, \
 "tig00000016" 1122257.0, \
 "tig00000012" 1165414.0, \
 "tig00000013" 1242949.0, \
 "tig00000027" 1299048.0, \
 "tig00000007" 1317397.0, \
 "tig00000018" 1359916.0, \
 "tig00000028" 1370544.0, \
 "tig00000014" 1431557.0, \
 "tig00000023" 1479673.0, \
 "tig00000022" 1507925.0, \
 "tig00000015" 1546991.0, \
 "tig00000021" 1574840.0, \
 "tig00000020" 1645355.0, \
 "tig00000010" 1700762.0, \
 "tig00000017" 1727200.0, \
 "tig00000003" 1737983.0, \
 "tig00000011" 1830640.0, \
 "" 1847903 \
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
set xrange [1:1847903]
set yrange [1:1847903]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/5A18NP1-JBb08-B/5A18NP1-JBb08-B.fplot" title "FWD" w lp ls 1, \
 "mummer_output/5A18NP1-JBb08-B/5A18NP1-JBb08-B.rplot" title "REV" w lp ls 2
