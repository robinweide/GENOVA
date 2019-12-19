#' Rescale a vector.
#'
#' @param x A vector which needs to be resized.
#' @param newrange A vector of the wanted length.
#' @return A rescaled vector.
rescale <- function(x, newrange = range(x)) {
  xrange <- range(x)
  mfac <- (newrange[2] - newrange[1]) / (xrange[2] - xrange[1])
  newrange[1] + (x - xrange[1]) * mfac
}
