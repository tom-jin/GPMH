//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
gsl_rng *r;

void rkernel(double *source, double *destination)
{
  *destination = *source + gsl_ran_gaussian(r, 1);
	*(destination+1) = *(source+1) + gsl_ran_gaussian(r, 4);
	return;
}

double dkernel(double *source, double *destination)
{
	return gsl_ran_gaussian_pdf(*source - *destination, 1) * gsl_ran_gaussian_pdf(*(source+1) - *(destination+1), 1);
}

double dtarget(double *value)
{
	if(*(value+1)< -100)
		return 0;
	if(*(value+1)> 100)
		return 0;
	return 0.1 * gsl_ran_gaussian_pdf(*value-12, 1) + 0.2 * gsl_ran_gaussian_pdf(*value-8, 1) + 0.3 * gsl_ran_gaussian_pdf(*value-4, 1) + 0.4 * gsl_ran_gaussian_pdf(*value, 1);// * gsl_ran_flat_pdf(*(value+1), -100, 100);
}

void cGPMH(double *samples, void *target_dummy, void *rkernel_dummy, void *dkernel_dummy, double *init, int *num_samples, int *N)
{
	int n = 0;
  int dim = 2;
	double *proposals, *acceptance;
	
	// Seed random number generator
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_mt19937);

	// malloc arrays
	proposals = (double*)malloc((*N+1) * dim * sizeof(double));
	acceptance = (double*)malloc((*N+1) * sizeof(double));

	while(n < *num_samples) {
		// MCMC Update
    
		if (n == 0)
			memcpy(proposals+(dim*(*N)), init, dim * sizeof(double));
		else
			gsl_ran_sample(r, proposals+(dim*(*N)), 1, samples+dim*(n-(*N)), (*N), dim * sizeof(double));

		for(int i=0; i < (*N); i++) {
			rkernel(proposals+(dim*(*N)), proposals+(dim*i));
		}

		// Calculate acceptance probability
		for(int i=0; i < (*N)+1; i++) {
			acceptance[i] = dtarget(proposals+(dim*i));
			for (int j=0; j < (*N)+1; j++) {
				if (i==j) continue;
				acceptance[i] *= dkernel(proposals+(dim*i), proposals+(dim*j));
			}
		}

		gsl_ran_discrete_t *g = gsl_ran_discrete_preproc((*N)+1, acceptance);

		// Sample
		for(int i=0; i < (*N); i++) {
			int a = gsl_ran_discrete(r, g);
			samples[dim*(n+i)] = proposals[dim*a];
			memcpy(samples+dim*(n+i), proposals+(dim*a), dim * sizeof(double));
		}
		gsl_ran_discrete_free(g);

		//Update counter.
		n += (*N);
	}
  samples[0] = 1234;
}
