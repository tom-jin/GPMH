#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return EXIT_FAILURE;}} while(0)

__device__ double cuda_gaussian_pdf(double x, double var)
{
  return (1/(var * sqrt(2 * M_PI)))*exp(- (x * x)/(2* var * var));
}

__global__ void cuda_rkernel(double *source, double *destination, double *random, int dim, int N)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < N)
	{
		for(int i=0; i<dim; i++)
		{
			*(destination+(tid*dim)+i) = *(source+i) + *(random+(tid*dim)+i);
		}
	}
}

__device__ double cuda_dkernel(double *source, double *destination, int dim)
{
	double density = 1.0;
	for(int i=0; i<dim; i++)
	{
		density *= cuda_gaussian_pdf(*(source+i) - *(destination+i), 1);
	}
	return density;
}

__device__ double cuda_dtarget(double *value)
{
	if(*(value+1)< -10)
		return 0;
	if(*(value+1)> 10)
		return 0;
	return 0.1 * cuda_gaussian_pdf(*value-12, 1) + 0.2 * cuda_gaussian_pdf(*value-8, 1) + 0.3 * cuda_gaussian_pdf(*value-4, 1) + 0.4 * cuda_gaussian_pdf(*value, 1);
}

__global__ void cuda_acceptance(double *proposals, double *acceptance, int dim, int N)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < N)
	{
    double p[] = {*(proposals+(dim*tid)), *(proposals+(dim*tid)+1)};
		double a = cuda_dtarget(p);
		for (int j=0; j < N; j++)
		{
			if (tid==j) continue;
			a *= cuda_dkernel(proposals+(dim*tid), proposals+(dim*j), dim);
		}
		acceptance[tid] = a;
	}
}

__global__ void cuda_sample(double *proposals, double *acceptance, double *sample, double *random, int dim, int Nprop, int Nsamp)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < Nsamp)
	{
		// This needs replacing with a reduction.
		double total = 0.0;
		double r = *(random+tid);
		*(sample + (dim*tid)) = r;

		for(int i=0; i<Nprop; i++) {
			total += *(acceptance+i);
		}
		r *= total;
		*(sample + (dim*tid) + 1) = *(acceptance+tid);
		total = 0.0;
		for(int i=0; i<Nprop; i++) {
			total += *(acceptance+i);
			if(r < total) {
				*(sample + (dim*tid)) = *(proposals + (dim*i));
				*(sample + (dim*tid)+1) = *(proposals + (dim*i)+1);
				break;
			}
		}
	}
}

void cudaGPMH(double *samples, void *target_dummy, void *rkernel_dummy, void *dkernel_dummy, double *init, int *num_samples, int *N)
{
  int n = 0;
	const int dim = 2;

	double *dev_samples, *dev_proposals, *dev_acceptance, *dev_rand;
	curandGenerator_t gen;

	// Seed random number generator
  srand((unsigned int)time(0));
	curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gen, (unsigned int) time(NULL));

	// malloc arrays
	samples = (double*)malloc(*num_samples * dim * sizeof(double));
	cudaMalloc((void**)&dev_samples, *num_samples * dim * sizeof(double));
	cudaMalloc((void**)&dev_proposals, ((*N)+1)*dim*sizeof(double));
	cudaMalloc((void**)&dev_acceptance, ((*N)+1)*sizeof(double));
	cudaMalloc((void**)&dev_rand, (*N)*dim*sizeof(double));

	cudaMemcpy(dev_proposals+(dim*(*N)), &init, dim * sizeof(double), cudaMemcpyDefault);

	while(n < *num_samples) {
		// MCMC Update
		curandGenerateNormalDouble(gen, dev_rand, (*N)*dim, 0.0, 1.0);
		cuda_rkernel<<<1, (*N)>>>(dev_proposals+((*N)*dim), dev_proposals, dev_rand, dim, (*N));

		// Calculate acceptance probability
		cuda_acceptance<<<1, (*N)+1>>>(dev_proposals, dev_acceptance, dim, (*N)+1);

		// Sample
		curandGenerateUniformDouble(gen, dev_rand, (*N));
		cuda_sample<<<1, (*N)>>>(dev_proposals, dev_acceptance, dev_samples+(n*dim), dev_rand, dim, (*N)+1, (*N));

    int i = rand() % (*N);
		cudaMemcpy(dev_proposals+(dim*(*N)), dev_samples+(dim*(n+i)), dim * sizeof(double), cudaMemcpyDefault);

    //Update counter.
		n += (*N);
	}
  cudaMemcpy(samples, dev_samples, *num_samples*dim*sizeof(double), cudaMemcpyDefault);

	cudaFree(dev_samples);
	cudaFree(dev_proposals);
	cudaFree(dev_acceptance);
    cudaFree(dev_rand);

	return;
}
