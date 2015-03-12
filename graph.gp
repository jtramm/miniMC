set logscale y
set format y "%.2E"
set logscale x
set xlabel "Energy [eV]" font 'Arial,25'
set ylabel "Flux" font 'Arial,25'
set title "H/U = 10 - Equal Space Lethargy Flux vs. Energy" font 'Arial,25'
plot "data.dat" u 1:2 w lines
