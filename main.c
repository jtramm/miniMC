#include"minimc_header.h"

#define NBINS 1000

int find_u_bin(double E, double Eo, double kill, double source_E)
{
	double low = log(Eo/source_E);
	double high = log(Eo/kill);

	double u = log(Eo/E);
	double val = (u-low) / (high-low);

	int bin = val * NBINS;
	return bin;
}

double find_u_bin_u(int i, double Eo, double kill, double source_E)
{
	double low = log(Eo/source_E);
	double high = log(Eo/kill);
	double del = (high-low)/NBINS;
	double bin = low+ i*del;
	return bin;
}

void print_flux(double Eo, double kill, double source_E, double * flux)
{
	FILE *fp = fopen("data.dat", "w");

	for(int i = 0; i < NBINS; i++ )
		fprintf(fp, "%e\t%e\n", find_u_bin_u(i,Eo,kill,source_E),flux[i]);
	fclose(fp);
}


int main(int argc, char * argv[])
{
	if( argc == 2 )
		omp_set_num_threads(atoi(argv[1]));

	double source_E = 250;
	int np = 100000;
	double temp = 300;
	double Eo = 1000; // Lethargy Reference
	int nr;
	Resonance * R = res_read(&nr);
	nr = 3;
	double HtoU = 10;
	XS h1 = {0.0, 0.0, 20.0, 20.0};
	double start = omp_get_wtime();
	double kill = 0.025;

	double flux[NBINS] = {0};

	int escapes = 0;
	#pragma omp parallel default (none) \
	shared(np, temp, nr, R, HtoU, h1, source_E, kill, Eo, flux) \
	reduction(+:escapes)
	{
		unsigned seed = (omp_get_thread_num()+1)*time(NULL) +1337;
		double local_flux[NBINS] = {0};

		#pragma omp for schedule(dynamic, 10)
		for( int i = 0; i < np; i++ )
		{
			double E = source_E;

			while(1)
			{
				// Kill Particle once it clears Resonances
				if( E < kill )
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
						int bin = find_u_bin(E,Eo,kill,source_E);
						local_flux[bin] += 1.0/Sigma_t;

						// Neutron Death
						break;
					}
					else // Elastic Scatter
					{
						// Flux Tally
						// Tally
						int bin = find_u_bin(E,Eo,kill,source_E);
						local_flux[bin] += 1.0/Sigma_t;
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
					// Flux Tally
					int bin = find_u_bin(E,Eo,kill,source_E);
					local_flux[bin] += 1.0/Sigma_t;
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
				flux[i] += local_flux[i];
		}
	} // Exit Parallel Region

	double end = omp_get_wtime();
	print_flux(Eo, kill, source_E, flux);
	printf("Threads: %d\n", omp_get_max_threads());
	printf("Particles: %.2e\n", (double) np);
	printf("Time = %.2lf sec\n", end-start);
	printf("Particles per second = %.2e\n", np / (end-start));
	// Compute Resonance Escape Probability
	printf("Resonance Escape Probability = %.2lf%%\n",(double)escapes/np * 100.0);

	return 0;
}
