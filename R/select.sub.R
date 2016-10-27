#' A quick lookup of square region from a HiC matrix.
#'
#' Extracts a symmetric matrix around the diagonal from *start* to *stop*.
#'
#' @param data HiC-pro matrix
#' @param start Start of index-range.
#' @param end End of index-range.
#' @return A data.table with normalised counts.
select.sub <- function( data, start, end ){
  x <- rep(start:end, end-start+1)
  y <- rep(start:end, each=end-start+1)
  data.sub <- data[list(x,y)]
  data.sub <- data.sub[!is.na(data.sub$V3)]
  return(data.sub)
}
