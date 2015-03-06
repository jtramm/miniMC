#include"minimc_header.h"

#define NBINS 1000;

int main(int argc, char * argv[])
{
	if( argc == 2 )
		omp_set_num_threads(atoi(argv[1]));

	double particle_source_E = 250;
	int np = 1000;
	double temp = 300;
	double Eo = 1000000; // Lethargy Reference
	int nr;
	Resonance * R = res_read(&nr);
	nr = 30;
	double HtoU = 10;
	XS h1 = {0.0, 0.0, 20.0, 20.0};
	double start = omp_get_wtime();

	int escapes = 0;
	#pragma omp parallel default (none) \
	shared(np, temp, nr, R, HtoU, h1, partile_source_E) \
	reduction(+:escapes)
	{
		unsigned seed = (omp_get_thread_num()+1)*time(NULL) +1337;

		double u_grid[NBINS];
		double flux[NBINS];
		for( int i = 0; i < NBINS; i++ )
		{
			double bot = log( Eo / particle_source_E ); 
			double del = bot / NBINS;
			u_grid[i] = bot + i*del;
		}

		#pragma omp for schedule(dynamic, 10)
		for( int i = 0; i < np; i++ )
		{
			// Particle Born @ 250 eV
			double E = particle_source_E;

			while(1)
			{
				// Kill Particle once it clears Resonances
				if( E < 0.6 )
				{
					escapes++;
					break;
				}

				// Travel some distance (don't care)

				// Get XS's
				XS u238 = calculate_XS( E, temp, R, nr );

				// Calculate total XS
				double Sigma_t = u238.sigma_t + HtoU * h1.sigma_t;

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
						// Tally

						// Neutron Death
						break;
					}
					else // Elastic Scatter
					{
						// Scatter
						double alpha = 0.9833336251;
						double lowE = alpha*E;
						double delta = E - lowE;
						r = (double) rand_r(&seed) / RAND_MAX;
						E = r*delta + lowE;

						// Flux Tally
					}
				}
				else // H1 Scatter
				{
					// Sample New Energy
					r = (double) rand_r(&seed) / RAND_MAX;
					E = r*E;

					// Flux Tally
				}

			} // Neutron Life Loop

		}// Particles Loop
	} // Exit Parallel Region

	double end = omp_get_wtime();
	printf("Threads: %d\n", omp_get_max_threads());
	printf("Time = %.2lf sec\n", end-start);
	printf("Particles per second = %.2e\n", np / (end-start));
	// Compute Resonance Escape Probability
	printf("Resonance Escape Probability = %.2lf%%\n",(double)escapes/np * 100.0);

	return 0;
}
