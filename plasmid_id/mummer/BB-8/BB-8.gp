set terminal png tiny size 800,800
set output "mummer_output/BB-8/BB-8.png"
set xtics rotate ( \
 "tig00000009" 1.0, \
 "tig00000016" 63111.0, \
 "tig00000006" 111705.0, \
 "tig00000002" 169461.0, \
 "tig00000014" 203424.0, \
 "tig00000007" 258724.0, \
 "tig00000013" 306479.0, \
 "tig00000005" 368059.0, \
 "tig00000003" 426961.0, \
 "tig00000001" 522468.0, \
 "tig00000011" 1445985.0, \
 "tig00000012" 1462817.0, \
 "tig00000004" 1520501.0, \
 "tig00000010" 1575851.0, \
 "tig00000015" 1633706.0, \
 "tig00000008" 1669939.0, \
 "" 1731897 \
)
set ytics ( \
 "tig00000009" 1.0, \
 "tig00000016" 63111.0, \
 "tig00000006" 111705.0, \
 "tig00000002" 169461.0, \
 "tig00000014" 203424.0, \
 "tig00000007" 258724.0, \
 "tig00000013" 306479.0, \
 "tig00000005" 368059.0, \
 "tig00000003" 426961.0, \
 "tig00000001" 522468.0, \
 "tig00000011" 1445985.0, \
 "tig00000012" 1462817.0, \
 "tig00000004" 1520501.0, \
 "tig00000010" 1575851.0, \
 "tig00000015" 1633706.0, \
 "tig00000008" 1669939.0, \
 "" 1731897 \
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
set xrange [1:1731897]
set yrange [1:1731897]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/BB-8/BB-8.fplot" title "FWD" w lp ls 1, \
 "mummer_output/BB-8/BB-8.rplot" title "REV" w lp ls 2
