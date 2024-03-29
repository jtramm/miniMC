COMPILER    = intel

program = minimc

source = \
main.c \
resonances.c \
FNF.c \
Faddeeva.c

obj = $(source:.c=.o)

# Standard Flags
CFLAGS := -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  LDFLAGS += 
  CFLAGS += -Ofast -ffast-math -fopenmp -Wall -flto
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  LDFLAGS += 
  CFLAGS += -O3 -xhost -ansi-alias -no-prec-div -DINTEL -openmp -ipo
endif

$(program): $(obj) minimc_header.h
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj) data.dat
run:
	./$(program)
graph:
	gnuplot graph.gp
edit:
	vim -p $(source) minimc_header.h
