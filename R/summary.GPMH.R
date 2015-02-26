#' @export
summary.GPMH <- function(object, ...) {
  return(describe(obj$x)[c(2:5, 8:12)])
}