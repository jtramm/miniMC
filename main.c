#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

// Direction: 0 = left, 1 = right
// Region: 1 = left, 2 = right

int main(void)
{
	unsigned int seed = time(NULL);

	int n_particles = 1000;

	int n_absorbed = 0;
	int n_scattered = 0;

	int global_tally[600] = {0};
	int global_variance[600] = {0};
	int local_tally[600] = {0};
	int local_variance[600] = {0};
	int particle_tally[600] = {0};

	float mean[600] = {0};
	float variance[600] = {0};

	// Loop over particles
	for( int i = 0; i < n_particles; i++ )
	{
		// Particle State
		float x;
		float direction;
		int region;
		int alive = 1;

		// Birth Particle Location
		x = ( (float) rand_r(&seed) / RAND_MAX ) * 6.f;
		if( x < 2.f )
			region = 1;
		else
			region = 2;

		while(alive)
		{
			// Sample Particle Direction
			direction = rand_r(&seed) % 2; 

			// Set XS's
			float sigma_t;
			float sigma_a;
			if( region == 1 )
			{
				sigma_t = 1.f;
				sigma_a = 0.5f;
			}
			else
			{
				sigma_t = 1.5f;
				sigma_a = 1.2f;
			}

			// Sample Flight Distance
			float dist = -logf( (float) rand_r(&seed) / RAND_MAX  ) / sigma_t;

			// Move Particle
			if( region == 1 )
			{
				if( direction = 0 )
					if( dist > x ) // Particle Escapes Left
					{
						alive = 0;
						break;
					}
					else
						if( dist > 2.f - x ) // Particle Travels 1 -> 2
						{
							x = 2.f;
							region = 2;
							break;
						}
			}
			else
			{
				if( direction = 0 )
					if( dist > x - 2.f ) // Particle Travels 2 -> 1
					{

						x = 2.f;
						region = 1;
						break;
					}
					else
						if( dist > 6.f - x ) // Particle Escapes Right
						{
							alive = 0;
							break;
						}
			}


			// Sample Interaction
			float interaction = (float) rand_r(&seed) / RAND_MAX * sigma_t;
			if( interaction < sigma_a )
			{
				n_absorbed++;
				alive = 0;
				break;
			}
			else
			{
				n_scattered++; // also need to tally x location
				particle_tally[(int) (x*100)]++;
			}
		}

		// Accumulate Tallies for Neutron history into Local running tally
		for( int j = 0; j < 600; j++ )
		{
			local_tally[j] += particle_tally[j];
			local_variance[j] += particle_tally[j] * particle_tally[j];
			particle_tally[j] = 0;
		}

	}

	// Global Accumulation of Tallies
	for( int i = 0; i < 600; i++ )
	{
		global_tally[i] = local_tally[i];
		global_variance[i] = local_variance[i];
	}

	// Compute statistics
	for( int i = 0; i < 600; i++ )
	{
		mean[i] = (float) global_tally[i] / n_particles;
		variance[i] = sqrtf( 1.f / (n_particles-1) * ( (float) global_variance[i] / n_particles - mean[i] * mean[i] )  ); 
	}
	
	for( int i = 0; i < 600; i ++ )
		printf("Bin: %d\tMean: %f\tVariance: %f\n", i, mean[i], variance[i]);

	return 0;
}
