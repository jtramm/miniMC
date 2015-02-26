set style data lines
set title "Variance of Collision Reaction Rate In Slab - 1000 bins - 100,000,000 Particles"
set xlabel "X-Coordinate (cm)"
set ylabel "Variance (cm^-1)"
plot [0:6] "variance.dat" using 1:2 lt rgb "blue" 
