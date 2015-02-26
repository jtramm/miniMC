set style data lines
set title "Collision Reaction Rate In Slab - 1000 bins - 100,000,000 Particles"
set xlabel "X-Coordinate (cm)"
set ylabel "Collision Reaction Rate (cm^-1)"
plot [0:6] "reaction_rate.dat" using 1:2
