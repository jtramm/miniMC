set style data lines
set title "Variance In Slab"
set xlabel "X-Coordinate (cm)"
set ylabel "Variance (cm^-1)"
plot [0:6] "variance.dat" using 1:2
