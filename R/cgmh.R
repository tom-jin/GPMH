#' Wrapper function for C code.
#' 
#' @param A Vector of returned acceptance values.
#' @param Y Vector of proposed moves.
#' @param N Number of parallel samples.
#' @return The sum of \code{x} and \code{y}.
cgmh <- function(A, Y, N) {
  .C("Cgmh", as.double(A), as.double(Y), as.integer(N))[[1]]   
}