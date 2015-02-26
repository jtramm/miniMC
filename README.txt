1D homogenous Monte Carlo for 22.211 PSet 1
John Tramm

You can compile and run using the included makefile

>$ make
>$ make run

Requires an OpenMP supporting C compiler (i.e., anything but Apple clang).
This is for threading purposes.

Graphs for collision rate and variance can be generated as:

>$ make graph
