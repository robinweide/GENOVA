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


# taken from ggplot
try_require <- function(package, fun, source = NULL) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }

  if(source == 'BIOC'){
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from Bioconductor and try again.", call. = FALSE)
  } else   if(source == 'github'){
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from github and try again.", call. = FALSE)
  } else {
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install and try again.", call. = FALSE)
  }

}
