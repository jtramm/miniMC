COMPILER    = gnu

# Standard Flags
CFLAGS := -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  LDFLAGS += -fopenmp
  CFLAGS += -Ofast -ffast-math -ftree-vectorize -msse2
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  LDFLAGS += -openmp
  CFLAGS += -O3 -xhost -ansi-alias -no-prec-div -DINTEL
endif

all:
	$(CC) $(CFLAGS) main.c -o miniMC $(LDFLAGS)

clean:
	rm -f miniMC
