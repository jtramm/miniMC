set style data lines
set title "Scattering Reaction Rate In Slab"
set xlabel "X-Coordinate (cm)"
set ylabel "Reaction Rate (cm^-1)"
plot [0:6] "reaction_rate.dat" using 1:2
