#ifndef __HEADER
#define __HEADER
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<omp.h>
#include"Faddeeva.h"

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
	double delta_u; // Lethargy Domain
	double temp;                         // Temperature
	Resonance * R;                // Resonances
	int nr;                             // Number of Resonances
} Input;

Resonance * res_read(int * n_resonances);
double complex FNF( double complex Z );
XS FNF_calculate_XS( double E, double temp, Resonance * R, int nr );
XS MIT_calculate_XS( double E, double temp, Resonance * R, int nr );
void res_out( XS * xs, int gp );
void graph_driver(void);
int find_u_bin(double E, Input input);
int find_E_bin(double E, Input input);
double find_u_bin_u(int i, Input input);
double find_u_bin_E( int i, Input input );
double find_E_bin_E(int i, Input input);
void print_flux(Input input, double * flux);
void tally( double E, Input input, double * flux, double Sigma_t );
void abs_tally( double E, Input input, double * flux, double * Ra, double Sigma_t );
void calculate_sigma_a(double E_low, double E_high,
		double * flux, double * Ra, Input input );
double run_slowing_down_problem(long HtoU, long np, XS (*Get_XS)( double, double, Resonance *, int) );
#endif
