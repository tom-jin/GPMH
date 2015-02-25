#' @title Generalised Metropolis-Hastings algorithm (R version) 
#'
#'  @description Runs the Generalised Metropolis-Hastings algorithm for sampling
#'  from a given target distribution. At each iteration, the algorithm draws N 
#'  new points from the proposal distribution given the current state, then it 
#'  computes the likelihoods of all N+1 points (the N new points plus the 
#'  current state). Then, these N+1 points are sampled with replacement 
#'  according to their likelihoods and the new sample of N+1 points is obtained.
#'  Finally, a state from these N+1 points is randomly chosen to generate the N 
#'  new points in the successive iteration.
#'  @param target The target distribution that the algorithm aims at sampling from.
#'  @param kernel Input function to draw samples from the proposal distribution.
#'  @param dkernel Input function to evaluate the proposal distribution pointwise.
#'  @param init.state Sets the starting value for the first point of the Markov Chain. 
#'  @param n Number of total iterations of the algorithm.
#'  @param N Number of samples to draw from the proposal distribution at each step.
#'  @return Returns a vector of \code{n*(N+1)} samples from the target distribution
#'  @example demo/DemoR.R
#'  @export rGPMH

rGPMH <- function(target, kernel, dkernel, init.state, n, N = 8) {
  d <- length(init.state)
  I <- 1 
  X <- matrix(NA, nrow = n*(N + 1), ncol = d)
  X[1, ] <- init.state
  Y <- matrix(NA, nrow = N+1, ncol = d)
  Y[1, ] <- X[1, ]
  
  for(j in 2:(N+1)) {
    Y[j, ] <- kernel(X[1, ]) # generate N new points from the proposal
  }
  
  K <- numeric(N+1)
  A <- numeric(N+1) 
  
  for (i in 1:n) { 
    for(j in 1:(N+1)) {
      for (k in 1:(N+1)) {
        K[k] <- dkernel(Y[j, ], Y[k, ])  
      }
      A[j]<- prod(K[-j]) * target(Y[j, ]) # likelihood of each of the N+1 points
    }
    
    # sample from the N+1 points to obtain the new N+1 MCMC samples 
    X[((i-1)*N+i):(i*(N + 1)), ] <- Y[sample(seq(1:(N+1)), replace=TRUE, prob=A), ]
    
    # Now we set Y for the next iteration:
    I <- sample(((i-1)*N+i):(i*(N + 1)), 1)  # draw uniformly the index I...
    Y[1, ] <- X[I, ] # ...of the X[I,] that will generate N new points
    for(j in 2:(N+1)) {
      Y[j, ] <- kernel(X[I, ]) # generate N new points
    }
  }
  
  output <- list(x = X[-1, ], n = n, call = match.call())
  class(output) <- "GPMH"
  return(output)
}
