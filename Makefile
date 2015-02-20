all:
	gcc -std=gnu99 -O5 -g main.c -o miniMC -lm -fopenmp
clean:
	rm -f miniMC
