#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

// Direction: 0 = left, 1 = right
// Region: 1 = left, 2 = right

#define NBINS 1000

float rn(unsigned long * seed);
void plot( int nbins, float * mean, float * variance );

int main(void)
{
	long n_particles = 100000000;
	int global_tally[NBINS] = {0};
	int global_variance[NBINS] = {0};

	omp_set_num_threads(omp_get_num_procs());

	double start = omp_get_wtime();

	long n_collisions = 0;

	#pragma omp parallel default(none) \
	shared(n_particles, global_tally, global_variance) \
	reduction(+:n_collisions) 
	{
		// stack is faster for small # of bins
		int local_tally[NBINS] = {0};
		int local_variance[NBINS] = {0};
		int particle_tally[NBINS] = {0};

		/*
		#ifdef INTEL
		int * local_tally = (int *) _mm_malloc( NBINS * sizeof(int), 64 );
		memset(local_tally, 0, NBINS * sizeof(int));
		int * local_variance = (int *) _mm_malloc( NBINS * sizeof(int), 64 );
		memset(local_variance, 0, NBINS * sizeof(int));
		int * particle_tally = (int *) _mm_malloc( NBINS * sizeof(int), 64 );
		memset(particle_tally, 0, NBINS * sizeof(int));
		#else
		//int * local_tally = (int *) calloc( NBINS, sizeof(int) );
		//int * local_variance = (int *) calloc( NBINS, sizeof(int) );
		//int * particle_tally = (int *) calloc( NBINS, sizeof(int) );
		#endif
		*/

		unsigned long seed = 1337 + time(NULL) + omp_get_thread_num();

		#pragma omp for schedule(dynamic, 100)
		// Loop over particles
		for( long i = 0; i < n_particles; i++ )
		{
			// Particle State
			float x;
			float direction;
			int region;

			// Birth Particle Location
			x = (float) rn(&seed) * 2.f;
			if( x < 2.f )
				region = 1;
			else
				region = 2;

			// Birth Particle Direction
			direction = rn(&seed) * 2.f - 1.f; 

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

	// Compute statistics
	for( int i = 0; i < NBINS; i++ )
	{
		mean[i] = (float) global_tally[i] / n_particles;
		variance[i] = sqrtf( 1.f / (n_particles-1) * ( (float) global_variance[i] / n_particles - mean[i] * mean[i] )  ); 
	}


	double end = omp_get_wtime();

	printf("Neutrons:   %g\n", (float) n_particles);
	printf("Runtime:    %.3lf seconds\n", end-start);
	printf("Neutrons/s: %g\n", n_particles/(end-start));

	plot( NBINS, mean, variance );

	return 0;
}

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
	FILE * fp = fopen("data.dat", "w");
	fprintf(fp, "x\tmean\tmin\tmax\n");
	for( int i = 0; i < nbins; i++ )
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n",
				6.f/nbins * i,
				mean[i],
				mean[i] - variance[i]/2.f,
				mean[i] + variance[i]/2.f );
	}
	fclose(fp);
}
