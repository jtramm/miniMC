#include"minimc_header.h"

#define NBINS 1000

int main(int argc, char * argv[])
{
	if( argc == 2 )
		omp_set_num_threads(atoi(argv[1]));

	XS h1 = {0.0, 0.0, 20.0, 20.0};

	double flux[NBINS] = {0};
	double group[3] = {0};

	Input input;
	input.np = 10000000;                      // Number of Particles
	input.HtoU = 10;                          // H1 to U-238
	input.Eo = 1000;                          // Lethargy Reference E
	input.kill = 0.025;                       // Kill Energy
	input.source_E = 250;                     // Source Energy
	input.low_u = log(input.Eo/input.source_E);       // Lethargy Low End
	input.delta = log(input.Eo/input.kill) - input.low_u; // Lethargy Domain
	input.temp = 300;                         // Temperature
	input.R = res_read(&input.nr);                // Resonances
	input.nr = 3;                             // Number of Resonances

	int escapes = 0;
	double start = omp_get_wtime();
	#pragma omp parallel default (none) \
	shared(input, h1, flux, group) \
	reduction(+:escapes)
	{
		unsigned seed = (omp_get_thread_num()+1)*time(NULL) +1337;
		double local_flux[NBINS] = {0};
		double local_group[3] = {0};

		#pragma omp for schedule(dynamic, 10)
		for( int i = 0; i < input.np; i++ )
		{
			double E = input.source_E;

			while(1)
			{
				// Kill Particle once it clears Resonances
				if( E < input.kill )
				{
					escapes++;
					break;
				}

				// Get XS's
				XS u238 = calculate_XS( E, input.temp, input.R, input.nr );

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
						tally(E, input, local_flux, local_group, Sigma_t);
						// Neutron Death
						break;
					}
					else // Elastic Scatter
					{
						tally(E, input, local_flux, local_group, Sigma_t);
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
					tally(E, input, local_flux, local_group, Sigma_t);
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
				flux[i] += local_flux[i]/input.np;
			for( int i = 0; i < 3; i++ )
				group[i] += local_group[i]/input.np;
		}
	} // Exit Parallel Region

	double end = omp_get_wtime();
	print_flux(input, flux);
	printf("Threads: %d\n", omp_get_max_threads());
	printf("Particles: %.2e\n", (double) input.np);
	printf("Time = %.2lf sec\n", end-start);
	printf("Particles per second = %.2e\n", input.np / (end-start));
	// Compute Resonance Escape Probability
	printf("Resonance Escape Probability = %.2lf%%\n",(double)escapes/input.np * 100.0);
	for( int i = 0; i < 3; i++ )
		printf("Group %d XS: %lf\n", i, (double) group[i]);
	return 0;
}

int find_u_bin(double E, Input input)
{
	double u = log(input.Eo/E);
	double val = (u-input.low_u) / (input.delta);
	int bin = val * NBINS;
	return bin;
}

double find_u_bin_u(int i, Input input)
{
	double del = (input.delta)/NBINS;
	double bin = input.low_u + i*del;
	return bin;
}

void print_flux(Input input, double * flux)
{
	FILE *fp = fopen("data.dat", "w");

	for(int i = 0; i < NBINS; i++ )
		fprintf(fp, "%e\t%e\n", find_u_bin_u(i,input),flux[i]);
	fclose(fp);
}

void tally( double E, Input input, double * flux, double * group, double Sigma_t )
{
	// Flux Tally
	int bin = find_u_bin(E,input);
	flux[bin] += 1.0/Sigma_t;

	// Group XS Tally
	if( E > 6.0 && E <= 10.0 )
		group[0] += 1.0/Sigma_t;
	else if ( E > 10.0 && E <= 25.0 )
		group[1] += 1.0/Sigma_t;
	else if ( E > 25.0 && E <= 50.0 )
		group[2] += 1.0/Sigma_t;
}
