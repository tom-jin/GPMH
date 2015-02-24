// C chunk called by rcGMH to compute the vector A of likelihoods, given Y
// Current limitations: only works in dim = 1, with target N(10,1) and proposal N(x,1).

#include <R.h>
#include <Rmath.h>
#include <omp.h>

void Cgpmh(double *A, double *Y,int *N){  
  
  int j,k;
  double K;
  int n=N[0];
  
  for (j=0; j<(n+1); j++){
    
    K=1;
    
    #pragma omp parallel for reduction(*:K)
    for (k=0; k<j; k++){
      K = K * dnorm(Y[k], Y[j], 1,0);  
    }      
    
    #pragma omp parallel for reduction(*:K)
    for (k=j+1; k<(n+1); k++){
      K = K * dnorm(Y[k], Y[j], 1,0);  
    }
    
    A[j] = K*dnorm(Y[j],10,1,0);
  }
  
}