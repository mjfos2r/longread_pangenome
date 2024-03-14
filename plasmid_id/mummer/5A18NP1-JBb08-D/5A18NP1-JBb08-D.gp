set terminal png tiny size 800,800
set output "mummer_output/5A18NP1-JBb08-D/5A18NP1-JBb08-D.png"
set xtics rotate ( \
 "tig00000009" 1.0, \
 "tig00000013" 17409.0, \
 "tig00000012" 46736.0, \
 "tig00000019" 77189.0, \
 "tig00000021" 121131.0, \
 "tig00000014" 152499.0, \
 "tig00000008" 246084.0, \
 "tig00000018" 285785.0, \
 "tig00000004" 325767.0, \
 "tig00000017" 370220.0, \
 "tig00000006" 411116.0, \
 "tig00000011" 463482.0, \
 "tig00000007" 503640.0, \
 "tig00000016" 569636.0, \
 "tig00000002" 578967.0, \
 "tig00000020" 686350.0, \
 "tig00000010" 730558.0, \
 "tig00000015" 774631.0, \
 "tig00000001" 826323.0, \
 "" 1787237 \
)
set ytics ( \
 "tig00000009" 1.0, \
 "tig00000013" 17409.0, \
 "tig00000012" 46736.0, \
 "tig00000019" 77189.0, \
 "tig00000021" 121131.0, \
 "tig00000014" 152499.0, \
 "tig00000008" 246084.0, \
 "tig00000018" 285785.0, \
 "tig00000004" 325767.0, \
 "tig00000017" 370220.0, \
 "tig00000006" 411116.0, \
 "tig00000011" 463482.0, \
 "tig00000007" 503640.0, \
 "tig00000016" 569636.0, \
 "tig00000002" 578967.0, \
 "tig00000020" 686350.0, \
 "tig00000010" 730558.0, \
 "tig00000015" 774631.0, \
 "tig00000001" 826323.0, \
 "" 1787237 \
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
set xrange [1:1787237]
set yrange [1:1787237]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/5A18NP1-JBb08-D/5A18NP1-JBb08-D.fplot" title "FWD" w lp ls 1, \
 "mummer_output/5A18NP1-JBb08-D/5A18NP1-JBb08-D.rplot" title "REV" w lp ls 2
