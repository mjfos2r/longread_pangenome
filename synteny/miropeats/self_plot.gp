set terminal png tiny size 800,800
set output "self_plot.png"
set xtics rotate ( \
 "contig000001" 1, \
 "contig000002" 906672, \
 "contig000003" 966290, \
 "contig000004" 1021243, \
 "contig000005" 1075160, \
 "contig000006" 1128751, \
 "contig000007" 1159701, \
 "contig000008" 1190542, \
 "contig000009" 1221204, \
 "contig000010" 1251848, \
 "contig000011" 1282483, \
 "contig000012" 1312921, \
 "contig000013" 1343220, \
 "contig000014" 1373477, \
 "contig000015" 1403446, \
 "contig000016" 1432495, \
 "contig000017" 1459016, \
 "contig000018" 1480677, \
 "contig000019" 1502011, \
 "contig000020" 1518740, \
 "contig000021" 1533622, \
 "contig000022" 1536789, \
 "contig000023" 1537675, \
 "contig000024" 1538297, \
 "contig000025" 1538784, \
 "contig000026" 1539170, \
 "contig000027" 1539525, \
 "contig000028" 1539705, \
 "contig000029" 1539874, \
 "contig000030" 1540029, \
 "contig000031" 1540179, \
 "contig000032" 1540320, \
 "contig000033" 1540448, \
 "" 1540602 \
)
set ytics ( \
 "contig000001" 1, \
 "contig000002" 906672, \
 "contig000003" 966290, \
 "contig000004" 1021243, \
 "contig000005" 1075160, \
 "contig000006" 1128751, \
 "contig000007" 1159701, \
 "contig000008" 1190542, \
 "contig000009" 1221204, \
 "contig000010" 1251848, \
 "contig000011" 1282483, \
 "contig000012" 1312921, \
 "contig000013" 1343220, \
 "contig000014" 1373477, \
 "contig000015" 1403446, \
 "contig000016" 1432495, \
 "contig000017" 1459016, \
 "contig000018" 1480677, \
 "contig000019" 1502011, \
 "contig000020" 1518740, \
 "contig000021" 1533622, \
 "contig000022" 1536789, \
 "contig000023" 1537675, \
 "contig000024" 1538297, \
 "contig000025" 1538784, \
 "contig000026" 1539170, \
 "contig000027" 1539525, \
 "contig000028" 1539705, \
 "contig000029" 1539874, \
 "contig000030" 1540029, \
 "contig000031" 1540179, \
 "contig000032" 1540320, \
 "contig000033" 1540448, \
 "" 1540602 \
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
if(GPVAL_VERSION < 5) set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:1540602]
set yrange [1:1540602]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "self_plot.fplot" title "FWD" w lp ls 1, \
 "self_plot.rplot" title "REV" w lp ls 2
