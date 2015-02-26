#'  @title Metropolis-Hastings algorithm 
#'
#'  @description Runs the standard Metropolis-Hastings algorithm for sampling 
#'  from a given target distribution. 
#'  @param target The target distribution that the algorithm aims at sampling from.
#'  @param kernel Input function to draw samples from the proposal distribution.
#'  @param dkernel Input function to evaluate the proposal distribution pointwise.
#'  @param init.state Sets the starting value for the first point of the Markov Chain. 
#'  @param n Number of total iterations of the algorithm.
#'  @return Returns a vector of \code{n} samples from the target distribution.
#'  @example demo/DemoMH.R
#'  @export MH
MH <- function(target, kernel, dkernel, init.state, n)
{
  d <- length(init.state)
  X <- matrix(NA, nrow = n + 1, ncol = d)
  X[1, ] <- init.state
  Y <- kernel(X[1, ])
    for (i in seq(1, (n))) {
      Y <- kernel(X[i, ])
      a <- min(1, target(Y)*dkernel(Y, X[i, ])/(target(X[i, ])*dkernel(X[i, ], Y)))
      if (runif(1) < a) 
        X[i + 1, ] <- Y
      else X[i + 1, ] <- X[i, ]
    }

  output <- list(x = X[-1, ], n = n, call = match.call())
  class(output) <- "GPMH"
  return(output)
}