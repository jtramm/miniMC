// 22.211 PSET 1 - Problem 2 - Distribution Verification Code

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int main(void)
{
	double val = 0;
	long iterations = 10000000;
	unsigned int seed = 1337 * (time(NULL) + 13);
	for( int i = 0; i < iterations; i ++ )
	{
		double r = (double) rand_r(&seed) / RAND_MAX; 
		double part = sqrt(16 * r*r - 16 * r + 5 ) + 4 *r - 2;
		val += pow(part, 1.0/3.0) - pow(part,-1.0/3.0);
	}

	double average = val / iterations;

	printf("Average = %lf\n", average);

	return 0;
}
