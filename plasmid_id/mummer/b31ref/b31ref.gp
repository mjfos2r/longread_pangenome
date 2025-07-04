set terminal png tiny size 800,800
set output "refmums/b31ref/b31ref.png"
set xtics rotate ( \
 "NC_001318.1" 1.0, \
 "NC_000957.1" 910724.0, \
 "NC_001904.1" 915951.0, \
 "NC_001849.2" 925336.0, \
 "NC_000955.2" 942156.0, \
 "NC_001850.1" 960932.0, \
 "NC_001903.1" 985108.0, \
 "NC_001851.2" 1011605.0, \
 "NC_001852.1" 1039759.0, \
 "NC_001853.1" 1069524.0, \
 "NC_001854.1" 1098124.0, \
 "NC_000948.1" 1125446.0, \
 "NC_000949.1" 1156195.0, \
 "NC_000950.1" 1186417.0, \
 "NC_000951.1" 1216715.0, \
 "NC_000952.1" 1246552.0, \
 "NC_000953.1" 1277351.0, \
 "NC_000954.1" 1308235.0, \
 "NC_001855.1" 1338885.0, \
 "NC_001856.1" 1375733.0, \
 "NC_001857.2" 1414561.0, \
 "NC_000956.1" 1468217.0, \
 "" 1521208 \
)
set ytics ( \
 "NC_001318.1" 1.0, \
 "NC_000957.1" 910724.0, \
 "NC_001904.1" 915951.0, \
 "NC_001849.2" 925336.0, \
 "NC_000955.2" 942156.0, \
 "NC_001850.1" 960932.0, \
 "NC_001903.1" 985108.0, \
 "NC_001851.2" 1011605.0, \
 "NC_001852.1" 1039759.0, \
 "NC_001853.1" 1069524.0, \
 "NC_001854.1" 1098124.0, \
 "NC_000948.1" 1125446.0, \
 "NC_000949.1" 1156195.0, \
 "NC_000950.1" 1186417.0, \
 "NC_000951.1" 1216715.0, \
 "NC_000952.1" 1246552.0, \
 "NC_000953.1" 1277351.0, \
 "NC_000954.1" 1308235.0, \
 "NC_001855.1" 1338885.0, \
 "NC_001856.1" 1375733.0, \
 "NC_001857.2" 1414561.0, \
 "NC_000956.1" 1468217.0, \
 "" 1521208 \
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
set xrange [1:1521208]
set yrange [1:1521208]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "refmums/b31ref/b31ref.fplot" title "FWD" w lp ls 1, \
 "refmums/b31ref/b31ref.rplot" title "REV" w lp ls 2
