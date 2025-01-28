set terminal png tiny size 800,800
set output "mummer_output/5A18NP1-JBb08-E/5A18NP1-JBb08-E.png"
set xtics rotate ( \
 "tig00000014" 1.0, \
 "tig00000004" 62461.0, \
 "tig00000017" 140220.0, \
 "tig00000013" 197476.0, \
 "tig00000019" 243398.0, \
 "tig00000012" 316948.0, \
 "tig00000002" 361009.0, \
 "tig00000010" 452307.0, \
 "tig00000015" 500428.0, \
 "tig00000011" 534168.0, \
 "tig00000006" 578484.0, \
 "tig00000009" 629410.0, \
 "tig00000005" 656771.0, \
 "tig00000018" 703204.0, \
 "tig00000016" 718456.0, \
 "tig00000021" 778584.0, \
 "tig00000007" 850342.0, \
 "tig00000020" 910456.0, \
 "tig00000003" 1843855.0, \
 "tig00000008" 1862566.0, \
 "" 1893840 \
)
set ytics ( \
 "tig00000014" 1.0, \
 "tig00000004" 62461.0, \
 "tig00000017" 140220.0, \
 "tig00000013" 197476.0, \
 "tig00000019" 243398.0, \
 "tig00000012" 316948.0, \
 "tig00000002" 361009.0, \
 "tig00000010" 452307.0, \
 "tig00000015" 500428.0, \
 "tig00000011" 534168.0, \
 "tig00000006" 578484.0, \
 "tig00000009" 629410.0, \
 "tig00000005" 656771.0, \
 "tig00000018" 703204.0, \
 "tig00000016" 718456.0, \
 "tig00000021" 778584.0, \
 "tig00000007" 850342.0, \
 "tig00000020" 910456.0, \
 "tig00000003" 1843855.0, \
 "tig00000008" 1862566.0, \
 "" 1893840 \
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
set xrange [1:1893840]
set yrange [1:1893840]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/5A18NP1-JBb08-E/5A18NP1-JBb08-E.fplot" title "FWD" w lp ls 1, \
 "mummer_output/5A18NP1-JBb08-E/5A18NP1-JBb08-E.rplot" title "REV" w lp ls 2
