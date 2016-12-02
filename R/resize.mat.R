#' Rescale a matrix.
#'
#' @param mat A vector which needs to be resized.
#' @param ndim Dimentions of the wanted matrix: can be obtained by dim(matrixOfChoice).
#' @return A rescaled vector.
resize.mat <- function(mat, ndim=dim(mat)){
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  # interpolation
  ans[ncord] <- GENOVA::fields.interp.surface(obj, loc)
  ans
}
