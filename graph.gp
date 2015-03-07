set logscale y
set logscale x
set xlabel "Energy [eV]" font 'Arial,25'
set ylabel "Flux" font 'Arial,25'
set title "Equal Space Lethargy Flux vs. Energy" font 'Arial,30'
plot "data.dat" u 1:2 w lines
