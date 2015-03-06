#ifndef __HEADER
#define __HEADER
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<omp.h>

typedef struct{
	double Eo;
	double Tn;
	double Tg;
} Resonance;

typedef struct{
	double E;
	double sigma_g;
	double sigma_n;
	double sigma_t;
} XS;

typedef struct{
	long np;                 // Number of Particles
	int HtoU;                          // H1 to U-238
	double Eo;                          // Lethargy Reference E
	double kill;                       // Kill Energy
	double source_E;                     // Source Energy
	double low_u;       // Lethargy Low End
	double delta; // Lethargy Domain
	double temp;                         // Temperature
	Resonance * R;                // Resonances
	int nr;                             // Number of Resonances
} Variables;

Resonance * res_read(int * n_resonances);
double complex FNF( double complex Z );
XS calculate_XS( double E, double temp, Resonance * R, int nr );
void res_out( XS * xs, int gp );
void graph_driver(void);
#endif
