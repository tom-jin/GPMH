cgmh <- function(A, Y, N) {
  .C("Cgmh", as.double(A), as.double(Y), as.integer(N))[[1]]   
}