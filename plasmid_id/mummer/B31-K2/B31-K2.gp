set terminal png tiny size 800,800
set output "mummer_output/B31-K2/B31-K2.png"
set xtics rotate ( \
 "tig00000019" 1.0, \
 "tig00000011" 19608.0, \
 "tig00000025" 64202.0, \
 "tig00000004" 85678.0, \
 "tig00000015" 169761.0, \
 "tig00000023" 195847.0, \
 "tig00000002" 230784.0, \
 "tig00000016" 272590.0, \
 "tig00000022" 355632.0, \
 "tig00000005" 410857.0, \
 "tig00000017" 484782.0, \
 "tig00000008" 516945.0, \
 "tig00000013" 602310.0, \
 "tig00000012" 648695.0, \
 "tig00000026" 679657.0, \
 "tig00000014" 708810.0, \
 "tig00000006" 739982.0, \
 "tig00000020" 784610.0, \
 "tig00000010" 804129.0, \
 "tig00000009" 848998.0, \
 "tig00000021" 895360.0, \
 "tig00000024" 1818759.0, \
 "tig00000018" 1856068.0, \
 "" 1904926 \
)
set ytics ( \
 "tig00000019" 1.0, \
 "tig00000011" 19608.0, \
 "tig00000025" 64202.0, \
 "tig00000004" 85678.0, \
 "tig00000015" 169761.0, \
 "tig00000023" 195847.0, \
 "tig00000002" 230784.0, \
 "tig00000016" 272590.0, \
 "tig00000022" 355632.0, \
 "tig00000005" 410857.0, \
 "tig00000017" 484782.0, \
 "tig00000008" 516945.0, \
 "tig00000013" 602310.0, \
 "tig00000012" 648695.0, \
 "tig00000026" 679657.0, \
 "tig00000014" 708810.0, \
 "tig00000006" 739982.0, \
 "tig00000020" 784610.0, \
 "tig00000010" 804129.0, \
 "tig00000009" 848998.0, \
 "tig00000021" 895360.0, \
 "tig00000024" 1818759.0, \
 "tig00000018" 1856068.0, \
 "" 1904926 \
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
set xrange [1:1904926]
set yrange [1:1904926]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "mummer_output/B31-K2/B31-K2.fplot" title "FWD" w lp ls 1, \
 "mummer_output/B31-K2/B31-K2.rplot" title "REV" w lp ls 2
