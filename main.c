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

	int n_particles = 100;

	int n_absorbed = 0;
	int n_scattered = 0;

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
				n_scattered++; // also need to tally x location
		}
	}
	printf("%d neutrons run. %d scattered\n", n_particles, n_scattered);

	return 0;
}
