#include"minimc_header.h"

#define NBINS 1000


int main(int argc, char * argv[])
{
	double FNF, MIT, SNF, SIMD;
	XS (*Get_XS)( double, double, Resonance *, int);
	printf("Threads = %d\n", omp_get_max_threads());

	long np = 100000;
	printf("NP = %ld\n", np);
	printf("SNF Slowing Down==========================================\n");
	Get_XS = &Simple_FNF_calculate_XS;
	SNF = run_slowing_down_problem(100000, np, Get_XS);
	printf("FNF Slowing Down==========================================\n");
	Get_XS = &FNF_calculate_XS;
	FNF = run_slowing_down_problem(100000, np, Get_XS);
	printf("SIMD Slowing Down==========================================\n");
	Get_XS = &SIMD_FNF_calculate_XS;
	SIMD = run_slowing_down_problem(100000, np, Get_XS);
	printf("MIT Slowing Down==========================================\n");
	Get_XS = &MIT_calculate_XS;
	MIT = run_slowing_down_problem(100000, np, Get_XS);

	printf("SUMMARY===================================================\n");
	printf("SNF  Particles/sec: %.2lf\n", np / SNF);
	printf("FNF  Particles/sec: %.2lf\n", np / FNF);
	printf("MIT  Particles/sec: %.2lf\n", np / MIT);
	printf("SIMD Particles/sec: %.2lf\n", np / SIMD);
	printf("FNF  Advantage: %.2lf%%\n", ((MIT-FNF) / MIT) * 100.0);
	printf("SNF  Advantage: %.2lf%%\n", ((MIT-SNF) / MIT) * 100.0);
	printf("SIMD Advantage: %.2lf%%\n", ((MIT-SIMD) / MIT) * 100.0);
	return 0;
}

double run_slowing_down_problem(long HtoU, long np, XS (*Get_XS)( double, double, Resonance *, int) )
{
	printf("Hydrogen to Uranium Ratio (H/U): %ld\n", HtoU);

	XS h1 = {0.0, 0.0, 20.0, 20.0};
	double flux[NBINS] = {0};
	double Ra[NBINS] = {0};

	Input input;
	input.np = np;                              // Number of Particles
	input.HtoU = HtoU;                          // H1 to U-238
	input.source_E = 2000000;                   // Source Energy
	input.Eo = input.source_E;                  // Lethargy Reference E
	input.kill = 0.6;                           // Kill Energy
	input.low_u = log(input.Eo/input.source_E); // Lethargy Low End
	input.delta_u = log(input.Eo/input.kill) - input.low_u; // Lethargy Domain
	input.temp = 300;                           // Temperature
	input.R = res_read(&input.nr);              // Resonances
	printf("resonances: %d\n", input.nr);
	//input.nr = 3;                             // Number of Resonances

	int escapes = 0;
	double start = omp_get_wtime();
	#pragma omp parallel default (none) \
	shared(input, h1, flux, Ra, Get_XS) \
	reduction(+:escapes)
	{
		unsigned seed = (omp_get_thread_num()+1)*time(NULL) +1337;
		double local_flux[NBINS] = {0};
		double local_Ra[NBINS] = {0};

		// Loop over particles
		#pragma omp for schedule(dynamic, 10)
		for( int i = 0; i < input.np; i++ )
		{
			double E = input.source_E;

			// Particle Life Loop
			while(1)
			{
				// Kill Particle once it clears Resonances
				if( E < input.kill )
				{
					escapes++;
					break;
				}

				// Get XS's
				XS u238 = (*Get_XS)( E, input.temp, input.R, input.nr );

				// Calculate total XS
				double Sigma_t = u238.sigma_t + input.HtoU * h1.sigma_t;

				// Sample Collision Nucleus
				double r = (double) rand_r(&seed) / RAND_MAX;
				r *= Sigma_t;

				if( r <= u238.sigma_t ) // U-238 collision
				{
					// Sample Collison Type
					r = (double) rand_r(&seed) / RAND_MAX;
					r *= u238.sigma_t;
					if( r <= u238.sigma_g ) // Absorption
					{
						abs_tally(E, input, local_flux, local_Ra, Sigma_t);
						// Neutron Death
						break;
					}
					else // Elastic Scatter
					{
						tally(E, input, local_flux, Sigma_t);
						// Scatter
						double alpha = 0.9833336251;
						double lowE = alpha*E;
						double delta = E - lowE;
						r = (double) rand_r(&seed) / RAND_MAX;
						E = r*delta + lowE;

					}
				}
				else // H1 Scatter
				{
					tally(E, input, local_flux, Sigma_t);
					// Sample New Energy
					r = (double) rand_r(&seed) / RAND_MAX;
					E = r*E;
				}
			} // Neutron Life Loop
		}// Particles Loop

		// Accumulate Flux
		#pragma omp critical
		{
			for( int i = 0; i < NBINS; i++ )
			{
				flux[i] += local_flux[i]/input.np;
				Ra[i] += local_Ra[i] / input.np;
			}
		}
	} // Exit Parallel Region
	double end = omp_get_wtime();

	// Calculate Group XS's
	calculate_sigma_a(6,10,flux,Ra,input);
	calculate_sigma_a(10,25,flux,Ra,input);
	calculate_sigma_a(25,50,flux,Ra,input);

	// Print Flux to File
	print_flux(input, flux);
	// Compute Resonance Escape Probability
	printf("Resonance Escape Probability = %.2lf%%\n",(double)escapes/input.np * 100.0);
	printf("Particles: %.2e\n", (double) input.np);
	printf("Particles per second = %.2e\n", input.np / (end-start));
	free(input.R);

	return end-start;
}

void calculate_sigma_a(double E_low, double E_high,
		double * flux, double * Ra, Input input )
{
	long bin_low, bin_high;
	bin_high = find_u_bin(E_low, input);
	bin_low = find_u_bin(E_high, input);

	double sigma_a = 0;

	double total_Ra = 0;
	double total_flux = 0;
	for( int i = bin_low; i < bin_high; i++ )
	{
		total_Ra += Ra[i];
		total_flux += flux[i];
	}
	sigma_a = total_Ra / total_flux;

	printf("Group: (%2.0lf-%2.0lf) [eV]    Absorption XS: %8.3lf [b]\n",
			E_low, E_high, sigma_a);	
}

int find_u_bin(double E, Input input)
{
	double u = log(input.Eo/E);
	double val = (u-input.low_u) / (input.delta_u);
	int bin = val * NBINS;
	return bin;
}

int find_E_bin(double E, Input input)
{
	return (E / input.source_E) * NBINS;
}

double find_u_bin_u(int i, Input input)
{
	double del = (input.delta_u)/NBINS;
	double bin = input.low_u + i*del;
	return bin;
}

double find_u_bin_E( int i, Input input )
{
	double u = find_u_bin_u(i, input);
	double E = input.Eo / exp(u);
	return E;
}

double find_E_bin_E(int i, Input input)
{
	return (input.source_E / NBINS) * i;
}

void print_flux(Input input, double * flux)
{
	FILE *fp = fopen("data.dat", "w");

	for(int i = 0; i < NBINS; i++ )
		fprintf(fp, "%e\t%e\n", find_u_bin_E(i,input),flux[i]);
	fclose(fp);
}

void tally( double E, Input input, double * flux, double Sigma_t )
{
	// Flux Tally
	int bin = find_u_bin(E,input);
	flux[bin] += 1.0/Sigma_t;
}

void abs_tally( double E, Input input, double * flux, double * Ra, double Sigma_t )
{
	int bin = find_u_bin(E,input);
	// Flux Tally
	flux[bin] += 1.0/Sigma_t;
	// Group XS Tally
	Ra[bin] += 1.0;
}
