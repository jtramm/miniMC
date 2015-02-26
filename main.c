#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#define NBINS 1000

float rn(unsigned long * seed);
void plot( int nbins, float * mean, float * variance );

int main(void)
{
	long n_particles           = 10000000;
	int global_tally[NBINS]    = {0};
	int global_variance[NBINS] = {0};
	long n_collisions          = 0;
	double start               = omp_get_wtime();
	omp_set_num_threads(omp_get_num_procs());
	printf("Simulating %e particles in %d bins...\n", (float) n_particles,
			NBINS);
	printf("Beginning simulation on %d threads...\n", omp_get_num_procs());

	#pragma omp parallel default(none) \
	shared(n_particles, global_tally, global_variance) \
	reduction(+:n_collisions) 
	{
		// stack is faster for small # of bins
		int local_tally[NBINS]    = {0};
		int local_variance[NBINS] = {0};
		int particle_tally[NBINS] = {0};

		// Intialize Seed
		unsigned long seed = 1337 * time(NULL) + omp_get_thread_num();

		// Loop over particles
		#pragma omp for schedule(dynamic, 100)
		for( long i = 0; i < n_particles; i++ )
		{
			// Particle State
			float x, direction;
			int region;

			// Birth Particle Location
			x = (float) rn(&seed) * 2.f;
			if( x < 2.f )
				region = 1;
			else
				region = 2;

			// Birth Particle Direction
			direction = rn(&seed) * 2.f - 1.f; 

			// Particle Tracking Loop
			while(1)
			{
				// Set XS's
				float sigma_t, sigma_a;
				if( region == 1 )
				{
					sigma_t = 1.f; sigma_a = 0.5f;
				}
				else
				{
					sigma_t = 1.5f; sigma_a = 1.2f;
				}

				// Sample Flight Distance and Project onto X-axis
				float dist = -logf( rn(&seed) ) / sigma_t;
				dist = dist * direction;

				// Move Particle
				if( region == 1 )
				{
					x += dist;
					if( x < 0 ) // Particle Escapes Left
						break;
					if( x > 2.f ) // Particle Travels 1 -> 2
					{
						x = 2.f;
						region = 2;
						continue;
					}
				}
				else
				{
					x += dist;
					if( x > 6.f ) // Particle Escapes Left
						break;
					if( x < 2.f ) // Particle Travels 2 -> 1
					{
						x = 2.f;
						region = 1;
						continue;
					}
				}

				// Sample Interaction
				float interaction = rn(&seed) * sigma_t;
				if( interaction < sigma_a )
					break;
				else
				{
					particle_tally[(int) (x*(NBINS/6.f))]++;
					direction = rn(&seed) * 2.f - 1.f; 
					n_collisions++;
				}
			}

			// Accumulate Tallies for Neutron history into Local running tally
			#pragma simd
			for( int j = 0; j < NBINS; j++ )
			{
				int tally = particle_tally[j];
				particle_tally[j] = 0;
				local_tally[j] += tally;
				local_variance[j] += tally * tally;
			}

		}

		// Global Accumulation of Tallies
		#pragma omp critical
		{
			for( int i = 0; i < NBINS; i++ )
			{
				global_tally[i] += local_tally[i];
				global_variance[i] += local_variance[i];
			}
		}
	}

	float mean[NBINS] = {0};
	float variance[NBINS] = {0};
	float phi[NBINS] = {0};

	// Compute statistics
	for( int i = 0; i < NBINS; i++ )
	{
		mean[i] = (float) global_tally[i] / n_particles;
		variance[i] = sqrtf( 1.f / (n_particles-1) * 
				( (float) global_variance[i] /
				  n_particles - mean[i] * mean[i] )); 
	}

	// Compute Phi
	for( int i = 0; i < NBINS; i++ )
	{
		float V = 6.f / NBINS;
		if( V * i < 2.f )
			phi[i] = mean[i] * 1.f / V;
		else
			phi[i] = mean[i] * 1.f / ( V * 1.5 );
	}

	double end = omp_get_wtime();
	printf("Neutrons:   %g\n", (float) n_particles);
	printf("Collisions: %ld\n", n_collisions); 
	printf("Runtime:    %.3lf seconds\n", end-start);
	printf("Neutrons/s: %g\n", n_particles/(end-start));
	plot( NBINS, phi, variance );

	return 0;
}

// Park & Miller LCG from
// Numerical Recipes Vol. 2
float rn(unsigned long * seed)
{
    float ret;
    unsigned long n1;
    unsigned long a = 16807;
    unsigned long m = 2147483647;
    n1 = ( a * (*seed) ) % m;
    *seed = n1;
    ret = (float) n1 / m;
    return ret;
}

void plot( int nbins, float * mean, float * variance )
{
	printf("writing data to files...\n");
	FILE * fp_rr = fopen("reaction_rate.dat", "w");
	FILE * fp_var = fopen("variance.dat", "w");
	fprintf(fp_rr, "x\tmean\n");
	fprintf(fp_var, "x\tvariance\n");
	for( int i = 0; i < nbins; i++ )
	{
		fprintf(fp_rr, "%e\t%e\n",
				6.f/nbins * i,
				mean[i]);
		fprintf(fp_var, "%e\t%e\n",
				6.f/nbins * i,
				variance[i]);
	}
	fclose(fp_rr);
	fclose(fp_var);
}
