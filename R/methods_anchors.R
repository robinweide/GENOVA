# Strategies --------------------------------------------------------------

#' Anchor positions for PE-SCAn
#'
#' Takes all pairwise interactions of locations on the genome and maps these to
#' indices of the Hi-C matrix.
#'
#' @inheritParams PESCAn
#' @param ABS The indices slot of a GENOVA experiment.
#' @param RES The resolution slot of a GENOVA experiment.
#' @param mode A \code{character} vector of length 1 indicating which
#'   interaction to take. Possible values: \code{"cis"}, \code{"trans"} or
#'   \code{"both"}.
#'
#' @return A \code{matrix} with two columns.
#'
#' @details The resulting matrix contains ordered indices to the Hi-C matrix
#'   slot in the GENOVA experiment.
#'
#'   The \code{mode} argument determines what pairwise interactions are
#'   reported. \code{"cis"} returns pairwise interactions within a chromosome;
#'   \code{"trans"} gives these between chromosomes and \code{"bot"}.
#'
#' @seealso \code{\link[GENOVA]{PESCAn}} for context.
#'
#'   \code{\link[GENOVA]{bed2idx}} for general conversion of BED-like
#'   \code{data.frame}s to Hi-C indices.
#' @export
anchors_PESCAn <- function(ABS, RES, bed,
                           dist_thres = c(5e6L, Inf),
                           min_compare = 10L,
                           mode = c("cis", "trans", "both")) {
  mode <- match.arg(mode)

  # Match and clean indices
  newbed <- cbind(bed, idx = bed2idx(ABS, bed, mode = "centre"))
  newbed <- newbed[!is.na(newbed$idx), ]
  newbed <- newbed[order(newbed$idx), ]

  if (mode == "cis") {

    # Generate cis combinations
    idx <- split(newbed$idx, newbed[, 1])
    idx <- idx[lengths(idx) >= max(min_compare, 2)]
    if (length(idx) == 0) {
      stop("Too few comparable positions per chromosome.", call. = FALSE)
    }
    idx <- lapply(idx, function(x) t(combn(x, 2)))
    idx <- do.call(rbind, idx)
    is_cis <- TRUE
  } else {

    # Generate all combinations
    idx <- t(combn(newbed$idx, 2))
    is_cis <- newbed[match(idx, newbed$idx), 1]
    is_cis <- matrix(is_cis, ncol = 2)
    is_cis <- is_cis[, 1] == is_cis[, 2]

    if (mode == "trans") {

      # Exclude cis when trans
      idx <- idx[!is_cis, ]
      is_cis <- FALSE
    }
  }

  # Filtering distances only makes sense in cis
  if (mode != "trans") {
    # Covert dists to resolution space
    dist_thres <- dist_thres / RES
    dist_thres <- sort(dist_thres)

    # Exclude cis combinations based on min/max_dist
    cis <- idx[is_cis, ]
    dist <- abs(cis[, 1] - cis[, 2])
    cis <- cis[dist >= dist_thres[1] & dist <= dist_thres[2], ]

    # Recombine and order
    idx <- rbind(idx[!is_cis, ], cis)
    idx <- idx[order(idx[, 1], idx[, 2]), ]
  }

  idx
}

#' Anchor positions for APA
#'
#' Transforms a BEDPE formatted \code{data.frame} to indices of the Hi-C matrix.
#'
#' @inheritParams APA
#' @inheritParams anchors_PESCAn
#'
#' @return A \code{matrix} with two columns.
#'
#' @details The resulting matrix contains ordered indices to the Hi-C matrix
#'   slot in the GENOVA experiment.
#'
#'   The \code{mode} argument determines what entries are
#'   reported. \code{"cis"} returns entries within chromosomes;
#'   \code{"trans"} gives these between chromosomes and \code{"bot"}.
#'
#' @seealso \code{\link[GENOVA]{APA}} for context.
#'
#'   \code{\link[GENOVA]{bed2idx}} for general conversion of BED-like
#'   \code{data.frame}s to Hi-C indices.
#' @export
anchors_APA <- function(ABS, RES, bedpe,
                        dist_thres = c(0, Inf),
                        mode = c("both", "cis", "trans")) {
  mode <- match.arg(mode)

  # Convert bedpe to idx matrix
  newbed <- cbind(bedpe[, 1:6],
    idx1 = bed2idx(ABS, bedpe[, 1:3], mode = "centre"),
    idx2 = bed2idx(ABS, bedpe[, 4:6], mode = "centre")
  )
  newbed <- na.exclude(newbed)
  newbed <- newbed[!duplicated(newbed[, 7:8]), ]

  # Shortcut when no additional filtering is needed
  if (identical(dist_thres, c(0, Inf)) && mode == "both") {
    # Ordering
    idx <- cbind(
      pmin(newbed$idx1, newbed$idx2),
      pmax(newbed$idx1, newbed$idx2)
    )
    idx <- idx[order(idx[, 1]), ]
    return(idx)
  }

  newbed$is_cis <- newbed[, 1] == newbed[, 4]

  if (mode == "trans") {
    newbed <- newbed[!newbed$is_cis, ]
  } else {
    # Convert distances to bins
    dist_thres <- dist_thres / RES
    dist_thres <- sort(dist_thres)

    # Filter cis on distances
    cis <- newbed[newbed$is_cis, ]
    dist <- abs(cis$idx1 - cis$idx2)
    cis <- cis[dist >= dist_thres[1] & dist <= dist_thres[2], ]
    if (mode == "cis") {
      newbed <- cis
    } else {
      newbed <- rbind(newbed[!newbed$is_cis, ], cis)
    }
  }

  # Ordering
  idx <- cbind(
    pmin(newbed$idx1, newbed$idx2),
    pmax(newbed$idx1, newbed$idx2)
  )
  idx <- idx[order(idx[, 1]), ]
  return(idx)
}

#' Anchor positions for ATA
#'
#' Transforms a BED-formatted \code{data.frame} containing TAD positions into
#' indices to the Hi-C matrix.
#'
#' @inheritParams anchors_PESCAn
#' @param bed A \code{data.frame} with 3 columns in BED format, containing TAD
#'   boundary positions for one TAD at each row.
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum TAD sizes to include.
#' @param padding A \code{numeric} of length 1 to determine the padding around
#'   TADs, expressed in TAD widths.
#'
#' @return A \code{matrix} with two columns.
#'
#' @details The resulting matrix contains ordered indices to the Hi-C matrix
#'   slot in the GENOVA experiment.
#'
#' @seealso \code{\link[GENOVA]{ATA}} for context.
#'
#'   \code{\link[GENOVA]{bed2idx}} for general conversion of BED-like
#'   \code{data.frame}s to Hi-C indices.
#' @export
anchors_ATA <- function(ABS, bed,
                        dist_thres = c(225000, Inf),
                        padding = 1) {
  if (!inherits(bed, "data.frame")) {
    bed <- as.data.frame(bed)[, 1:3]
  }

  # Setup parameters
  width <- abs((bed[, 3] - bed[, 2]))
  mid <- round((bed[, 2] + bed[, 3]) / 2)
  keep <- width  > dist_thres[1] & width < dist_thres[2]

  # Resize regions
  bed <- data.frame(bed[,1],
                    mid - width * padding,
                    mid + width * padding)[keep, ]

  # Translate to Hi-C indices
  idx <- cbind(
    bed2idx(ABS, bed, mode = "start"),
    bed2idx(ABS, bed, mode = "end")
  )
  # Sort
  idx <- cbind(
    pmin(idx[, 1], idx[, 2]),
    pmax(idx[, 1], idx[, 2])
  )
  idx <- idx[order(idx[, 1]), ]

  # Attribute to let matrix lookup methods know it is performing ATA
  attr(idx, "type") <- "TADs"
  return(idx)
}

# Manipulations ---------------------------------------------------------

#' Shift anchors
#'
#' Shifts anchors upstream by a specified amount, unless they become out of
#' bounds. In that case, shift the anchors downstream.
#'
#' @inheritParams rep_mat_lookup
#' @inheritParams anchors_PESCAn
#' @param shift An \code{integer} of length 1 indicating how many bins the
#'   anchors should be shifted.
#'
#' @return A \code{matrix} with two columns.
#'
#' @details The resulting matrix contains indices to the Hi-C matrix slot in the
#'   GENOVA experiment.
#'
#'   An index is considered out of bounds when that index plus the shift size
#'   and maximum relative position would belong to a different chromosome.
#'
#' @seealso \code{\link[GENOVA]{PESCAn}} for context.
#'
#'   \code{\link[GENOVA]{anchors_filter_oob}} for general out of bounds
#'   filtering of anchors.
#'
#' @export
anchors_shift <- function(ABS, anchors, rel_pos, shift = 1) {
  # Translate indices to chromosomes
  chrom <- ABS[match(anchors, ABS[, 4]), 1]
  shifted <- ABS[match(anchors + shift + max(rel_pos), ABS[, 4]), 1]

  # Check chromosome is same after shift
  inbounds <- matrix(chrom == shifted, ncol = 2)
  inbounds <- apply(inbounds, 1, all)

  # Shift upstream unless out of bounds
  shifted <- ifelse(inbounds, 1, -1) * shift

  anchors + c(shifted, shifted)
}

# Quality control ---------------------------------------------------------

#' Filter anchors that become out of bounds.
#'
#' Discards anchors for which a lookup will result in any lookup bins that are
#' out of bounds.
#'
#' @inheritParams anchors_shift
#'
#' @return A \code{matrix} with two columns.
#'
#' @details The resulting matrix contains indices to the Hi-C matrix slot in the
#'   GENOVA experiment.
#'
#'   An index is considered out of bounds when that index plus the maximum, or
#'   minus the minimum, relative position would belong to a different chromosome.
#'
#' @export
anchors_filter_oob <- function(ABS, anchors, rel_pos) {

  if (!is.null(attr(anchors, "type"))) {
    if (attr(anchors, "type") == "TADs") {
      left  <- ABS[match(anchors[, 1], ABS[, 4]), 1]
      right <- ABS[match(anchors[, 2], ABS[, 4]), 1]
      keep <- left == right
      anchors <- anchors[keep, ]
      attr(anchors, "type") <- "TADs"
      return(anchors)
    }
  }

  # Match idx +/- relative position to chrom
  plus <- ABS[match(anchors + max(rel_pos), ABS[, 4]), 1]
  minus <- ABS[match(anchors + min(rel_pos), ABS[, 4]), 1]

  # Check wether chromosomes have changed
  inbounds <- matrix(plus == minus, ncol = 2)
  inbounds <- apply(inbounds, 1, all)

  # Return anchors that are not out of bounds
  anchors[inbounds, , drop = FALSE]
}
