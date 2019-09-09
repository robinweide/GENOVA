#' A quick lookup of square region from a HiC matrix.
#'
#' Extracts a symmetric matrix. *X* and *Y* renote the indexes to extracts.
#'
#' @param data HiC-pro matrix
#' @param X The X index-range.
#' @param Y The Y index-range.
#' @return A data.table with normalised counts.
#' @import data.table
select.sub.2D <- function(data, X, Y) {
  x <- rep(X[1]:X[length(X)], X[length(X)] - X[1] + 1)
  y <- rep(Y[1]:Y[length(Y)], each = length(Y[1]:Y[length(Y)]))
  data.sub <- data[base::list(x, y)]
  data.sub <- data.sub[!is.na(data.sub$V3)]
  data.sub
}
