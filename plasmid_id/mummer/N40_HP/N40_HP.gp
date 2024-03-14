set terminal png tiny size 800,800
set output "mummer_output/N40_HP/N40_HP.png"
set xtics rotate ( \
 "tig00000012" 1.0, \
 "tig00000011" 46635.0, \
 "tig00000014" 81122.0, \
 "tig00000009" 120961.0, \
 "tig00000016" 171765.0, \
 "tig00000018" 191800.0, \
 "tig00000013" 215255.0, \
 "tig00000020" 242402.0, \
 "tig00000001" 262851.0, \
 "tig00000007" 1163198.0, \
 "tig00000005" 1215895.0, \
 "tig00000010" 1258841.0, \
 "tig00000008" 1297378.0, \
 "tig00000002" 1337438.0, \
 "tig00000004" 1394646.0, \
 "tig00000003" 1447696.0, \
 "tig00000006" 1499379.0, \
 "tig00000015" 1559545.0, \
 "tig00000017" 1609805.0, \
 "" 1635730 \
)
set ytics ( \
 "tig00000012" 1.0, \
 "tig00000011" 46635.0, \
 "tig00000014" 81122.0, \
 "tig00000009" 120961.0, \
 "tig00000016" 171765.0, \
 "tig00000018" 191800.0, \
 "tig00000013" 215255.0, \
 "tig00000020" 242402.0, \
 "tig00000001" 262851.0, \
 "tig00000007" 1163198.0, \
 "tig00000005" 1215895.0, \
 "tig00000010" 1258841.0, \
 "tig00000008" 1297378.0, \
 "tig00000002" 1337438.0, \
 "tig00000004" 1394646.0, \
 "tig00000003" 1447696.0, \
 "tig00000006" 1499379.0, \
 "tig00000015" 1559545.0, \
 "tig00000017" 1609805.0, \
 "" 1635730 \
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
set xrange [1:1635730]
set yrange [1:1635730]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/N40_HP/N40_HP.fplot" title "FWD" w lp ls 1, \
 "mummer_output/N40_HP/N40_HP.rplot" title "REV" w lp ls 2
