all:
	gcc -O5 main.c -o miniMC -lm -fopenmp
clean:
	rm -f miniMC
