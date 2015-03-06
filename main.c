#include"minimc_header.h"


int main(void)
{
	int np = 1000;
	double temp = 300;
	int nr;
	Resonance * R = res_read(&nr);
	double HtoU = 100000;
	XS h1 = {0.0, 0.0, 20.0, 20.0};

	unsigned seed = time(NULL) +1337;

	int escapes = 0;

	for( int i = 0; i < np; i++ )
	{
		// Particle Born @ 250 eV
		double E = 250;

		while(1)
		{
			// Kill Particle once it clears Resonances
			if( E < 0.025 )
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

			if( r <= sigma_t.u238 ) // U-238 collision
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
}

return 0;
}
