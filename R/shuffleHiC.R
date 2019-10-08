#' Shuffle a symmetric matrix
#'
#' ShuffleHiC shuffles a matrix from `select_subset()`` diagonal-wise. That is, every diagonal in the matrix is shuffled independently.
#' @examples
#' \dontrun{
#' ## Shuffle a region on chromosome 1
#' HiCregion <- select_subset(EXP, chrom = "chr1", start = 15e6, end = 20e6)
#' HiCregion_shuffled <- shuffleHiC(HiCregion)
#'
#' # plot observed/expected
#' image(HiCregion$z / HiCregion_shuffled$z)
#' }
#' @param MAT A Hi-C matrix: can be produced by `select_subset()`.
#' @param symmetric In some special cases you might not want to get a symmetric matrix back. Set this argument on FALSE.
#' @return A shuffled version of the input: either a list with the shuffled `z` (when using `select_subset`) or a shuffled matrix.
#' @export
#'
#'
shuffleHiC <- function(MAT, symmetric = T) {
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  # works only with symmetric matrices, where the diagonal starts at 1,1 and ends at n,n

  # check if input is from select_subset (list with c(x),c(y), mat(z)) or just a matrix
  SS <- F
  MAT_BK <- MAT
  if (is.list(MAT)) {
    SS <- T
    MAT <- MAT$z
  }

  id <- which(MAT > -1, arr.ind = T)

  v <- MAT[MAT > -1]

  sh <- tapply(v, abs(id[, 1] - id[, 2]), sample)

  id.o <- id[order(abs(id[, 1] - id[, 2])), ]

  MAT.sh <- MAT

  MAT.sh[cbind(id.o[, 1], id.o[, 2])] <- unlist(sh)

  # make symmetric again
  if (symmetric) {
    MAT.sh[lower.tri(MAT.sh)] <- t(MAT.sh)[lower.tri(MAT.sh)]
  }

  # if from ss(), put shuffled Z back in Z
  if (SS) {
    MAT_BK$z <- MAT.sh
    MAT.sh <- MAT_BK
  }

  return(MAT.sh)
}
