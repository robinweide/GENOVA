#' Match bed-like entries to Hi-C bin indices
#'
#' @inheritParams PESCAn
#' @param mode A \code{character} of length 1 indicating what position of the
#'   \code{bed} argument to match with the indices. Possible values:
#'   \code{"center"}, \code{"start"} or \code{"end"}.
#'
#' @return An \code{integer} vector of length \code{nrow(bed)} and parallel to
#'   \code{bed} with indices to the Hi-C matrix.
#'
#' @details Out of bounds values are matched to nearest bin.
#' @export
bed2idx <- function(ABS, bed, mode = c("centre", "start", "end")) {
  if (!inherits(bed, "data.frame")) {
    bed <- as.data.frame(bed)
  }

  # American/British spelling
  mode <- gsub("center", "centre", mode)
  mode <- match.arg(mode, c("centre", "start", "end"))

  # Reformat bed depending on mode
  bed <- cbind.data.frame(
    V1 = bed[, 1],
    V2 = switch(mode,
      "centre" = (bed[, 2] + bed[, 3]) / 2,
      "start" = bed[, 2],
      "end" = bed[, 3]
    )
  )

  # Assign entries to shared chromosomes
  chroms <- intersect(ABS[, 1], bed[, 1])
  bed_group <- match(bed[, 1], chroms)
  ABS_group <- match(ABS[, 1], chroms)

  # Split by chromosome
  bed_chrom <- split(bed[, 2], bed_group)
  ABS_chrom <- split(ABS[, c(2, 4)], ABS_group)

  # Match bed entry to idx
  out <- mapply(function(i, j) {
    j[pmax(findInterval(i, j[, 1]), 1), 2]
  }, i = bed_chrom, j = ABS_chrom)
  unsplit(out, bed_group)
}

#' Upper triangle sparse symmetric triplet matrix to dense matrix
#'
#' Convenience function for returning square, dense symmetric matrices from the
#' upper triangle of triplet format symmetric matrices.
#'
#' @param x A \code{numeric} vector containing the data for matrix elements.
#' @param i A \code{integer} vector parallel to \code{x} containing row
#'   positions.
#' @param j A \code{integer} vector parallel to \code{x} containing column
#'   positions.
#' @param dim A \code{integer} of length 1 specifying the dimensions of the
#'   output matrix.
#' @param offset A \code{integer} of length 1 noting a potential offset in the
#'   \code{i} and \code{j} arguments.
#'
#'   Particularly useful when lookup up square regions around the diagonal in
#'   Hi-C data.
#'
#' @note For reasons of speed, no checks are performed wether the input is
#'   compatible with sensible output.
#'
#' @return A \code{dim * dim} sized \code{matrix}.
#'
#' @keywords internal
dt_matrix <- function(x, i, j, dim, offset) {
  m <- .Internal(matrix(0, dim, dim, FALSE, NULL, FALSE, FALSE))
  i <- .Internal(matrix(c(i, j, j, i) - offset,
                        2 * length(i), 2, FALSE, NULL, FALSE, FALSE))
  m[i] <- c(x, x)
  m
}
