#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>

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

Resonance * res_read(int * n_resonances);
double complex FNF( double complex Z );
XS calculate_XS( double E, double temp, Resonance * R, int nr );
void res_out( XS * xs, int gp );
void graph_driver(void);
