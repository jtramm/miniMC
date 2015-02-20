all:
	icc -std=gnu99 -O3 -xhost -ansi-alias -no-prec-div main.c -o miniMC -lm -openmp
clean:
	rm -f miniMC
