#include "minimc_header.h"

XS FNF_calculate_XS( double E, double temp, Resonance * R, int nr )
{
	double sigma_pot = 11.2934;
	double k = 8.6173324e-5;
	double A = 238.05078826;
	XS xs = {0};
	xs.E = E;

	for( int j = 0; j < nr; j++ )
	{
		double r = 2603911.0 / R[j].Eo * (A+1) / A;
		double q = 2.0 * sqrt(r * sigma_pot);
		double T = R[j].Tn + R[j].Tg;
		double x = 2.0 * (E - R[j].Eo) / T;
		double xi = T * sqrt(A / (4.0 * k * temp * R[j].Eo));
		double complex faddeeva_in = x + I;
		faddeeva_in *= xi;
		double complex faddeeva_out = xi * FNF( faddeeva_in);
		double psi = sqrt(M_PI) * creal(faddeeva_out); 
		double chi = sqrt(M_PI) * cimag(faddeeva_out);
		xs.sigma_g += R[j].Tn * R[j].Tg / (T*T) * sqrt(R[j].Eo / E) * r * psi;
		xs.sigma_n += R[j].Tn * R[j].Tn / (T*T) * ( r * psi + q * T/R[j].Tn * chi ); 
	}
	xs.sigma_n += sigma_pot;
	xs.sigma_t = xs.sigma_g + xs.sigma_n;

	return xs;
}

XS MIT_calculate_XS( double E, double temp, Resonance * R, int nr )
{
	double sigma_pot = 11.2934;
	double k = 8.6173324e-5;
	double A = 238.05078826;
	XS xs = {0};
	xs.E = E;

	for( int j = 0; j < nr; j++ )
	{
		double r = 2603911.0 / R[j].Eo * (A+1) / A;
		double q = 2.0 * sqrt(r * sigma_pot);
		double T = R[j].Tn + R[j].Tg;
		double x = 2.0 * (E - R[j].Eo) / T;
		double xi = T * sqrt(A / (4.0 * k * temp * R[j].Eo));
		double complex faddeeva_in = x + I;
		faddeeva_in *= xi;
		double complex faddeeva_out = xi * Faddeeva_w( faddeeva_in,0.0);
		double psi = sqrt(M_PI) * creal(faddeeva_out); 
		double chi = sqrt(M_PI) * cimag(faddeeva_out);
		xs.sigma_g += R[j].Tn * R[j].Tg / (T*T) * sqrt(R[j].Eo / E) * r * psi;
		xs.sigma_n += R[j].Tn * R[j].Tn / (T*T) * ( r * psi + q * T/R[j].Tn * chi ); 
	}
	xs.sigma_n += sigma_pot;
	xs.sigma_t = xs.sigma_g + xs.sigma_n;

	return xs;
}

void res_out( XS * xs, int gp )
{
	FILE * fp = fopen("data.dat", "w");
	for( int i = 0; i < gp; i++ )
	{
		fprintf(fp, "%e\t%e\t%e\t%e\n",
				xs[i].E,
				xs[i].sigma_g,
				xs[i].sigma_n,
				xs[i].sigma_t);
	}
	fclose(fp);
}

Resonance * res_read(int * n_resonances)
{
	FILE * fp = fopen("resonances.dat", "r");
	if( fp == NULL )
	{
		printf("ERROR! No resonances file found!\n");
		exit(1);
	}

	int lines = 0;
	char c;
	while ( (c = getc(fp)) != EOF)
	{
		if (c == '\n')
			lines++;
	}
	rewind(fp);

	Resonance * R = (Resonance *) malloc( lines * sizeof(Resonance));

	for( int i = 0; i < lines; i++ )
	{
		double dummy;
		double * dum = &dummy;
		fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				&R[i].Eo,
				dum,
				&R[i].Tn,
				&R[i].Tg,
				dum,
				dum);
	}

	fclose(fp);
	*n_resonances = lines;

	return R;
}
