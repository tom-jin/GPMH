#' @export
cudaGPMH <- function(target, rkernel, dkernel, init, n, N) 
  .C("cudaGPMH", as.double(rep(0.0, 2*n)), as.integer(0), as.integer(0), as.integer(0), as.double(init), as.integer(n), as.integer(N))[[1]]

