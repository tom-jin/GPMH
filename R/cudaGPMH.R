cGPMH <- function(target, rkernel, dkernel, init, n, N) 
  .C("cGPMH", as.double(rep(0.0, n)), as.integer(0), as.integer(0), as.integer(0), as.double(init), as.integer(n), as.integer(N))[[1]]

