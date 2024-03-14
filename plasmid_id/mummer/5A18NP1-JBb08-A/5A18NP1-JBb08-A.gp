set terminal png tiny size 800,800
set output "mummer_output/5A18NP1-JBb08-A/5A18NP1-JBb08-A.png"
set xtics rotate ( \
 "tig00000018" 1.0, \
 "tig00000014" 21131.0, \
 "tig00000009" 61606.0, \
 "tig00000015" 89620.0, \
 "tig00000010" 133826.0, \
 "tig00000016" 174522.0, \
 "tig00000003" 220812.0, \
 "tig00000021" 315726.0, \
 "tig00000002" 355223.0, \
 "tig00000020" 373765.0, \
 "tig00000011" 427786.0, \
 "tig00000022" 512644.0, \
 "tig00000007" 530016.0, \
 "tig00000008" 575705.0, \
 "tig00000019" 613513.0, \
 "tig00000012" 708245.0, \
 "tig00000017" 750847.0, \
 "tig00000001" 788062.0, \
 "tig00000013" 1742458.0, \
 "" 1811037 \
)
set ytics ( \
 "tig00000018" 1.0, \
 "tig00000014" 21131.0, \
 "tig00000009" 61606.0, \
 "tig00000015" 89620.0, \
 "tig00000010" 133826.0, \
 "tig00000016" 174522.0, \
 "tig00000003" 220812.0, \
 "tig00000021" 315726.0, \
 "tig00000002" 355223.0, \
 "tig00000020" 373765.0, \
 "tig00000011" 427786.0, \
 "tig00000022" 512644.0, \
 "tig00000007" 530016.0, \
 "tig00000008" 575705.0, \
 "tig00000019" 613513.0, \
 "tig00000012" 708245.0, \
 "tig00000017" 750847.0, \
 "tig00000001" 788062.0, \
 "tig00000013" 1742458.0, \
 "" 1811037 \
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
set xrange [1:1811037]
set yrange [1:1811037]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/5A18NP1-JBb08-A/5A18NP1-JBb08-A.fplot" title "FWD" w lp ls 1, \
 "mummer_output/5A18NP1-JBb08-A/5A18NP1-JBb08-A.rplot" title "REV" w lp ls 2
