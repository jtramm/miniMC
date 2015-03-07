set logscale y
set logscale x
set xlabel "Lethargy" font 'Arial,25'
set ylabel "Flux" font 'Arial,25'
plot "data.dat" u 1:2 w lines
