# User facing -------------------------------------------------------------

#' Repeated Hi-C matrix lookup
#'
#' Analyses like aggregrate peak analysis (APA) and paired-end spatial chromatin
#' analysis (PE-SCAn) require to look up regions of a Hi-C matrix repeatedly.
#' This function acts as a common method for the lookup and summary of such
#' repeated lookups for these analyses.
#'
#' @inheritParams PESCAn
#' @param anchors A \code{matrix} with two columns containing pre-computed
#'   anchor indices. See the \code{\link[GENOVA]{anchors}} documentation.
#' @param rel_pos An \code{integer} vector indicating relative positions in
#'   bins, indicating ranges around anchors to lookup.
#'
#' @return A \code{list} of the same length as \code{explist} wherein list
#'   elements contain the results of the repeated lookup per experiment.
#'
#' @details For each row in the \code{anchors} argument a region of the Hi-C
#'   matrix is looked up corresponding to that anchor. This data is then
#'   summarised by taking the mean of each position relative to the anchor
#'   across all the anchors.
#'
#'   Anchors are subject to a filtering step wherein anchors are discarded when
#'   they are within \code{rel_pos} range of a chromosome start or end. This
#'   ensures the anchors all report data from the same chromosome in the x- or
#'   y-direction.
#'
#'   For shifted anchors, an attempt is made to shift the anchors in the
#'   opposite direction before they are discarded.
#'
#'   When a region corresponding to a non-shifted anchor is looked up and is
#'   found to have no contacts within that region, it is discarded. Shifted
#'   regions are only discarded when the corresponding non-shifted anchor is
#'   discarded.
#'
#'   Anchors typically contain a '\code{type}' attribute which informs
#'   \code{rep_mat_lookup} how the lookup should occur. The \code{APA},
#'   \code{PESCAn} and \code{ARA} anchors look up regions of dimensions
#'   \code{length(rel_pos)} x \code{length(rel_pos)}. The \code{ARA} anchors
#'   transpose these square regions when given a '\code{-}' direction before
#'   summary occurs. The \code{ATA} anchors look up \code{anchors[, 2] -
#'   anchors[, 1]} sized square regions and resizes these to a
#'   \code{max(rel_pos)} square region through bilinear interpolation before
#'   summary.
#'
#' @section Resolution recommendation: 10kb-40kb
#'
#' @export
#'
#' @family aggregate repeated matrix lookup analyses
rep_mat_lookup <- function(explist, anchors, rel_pos, shift = 0,
                           outlier_filter = c(0, 1), raw = FALSE) {

  # Finish off anchors
  anchors <- anchors_finish(explist[[1]]$IDX, anchors, rel_pos, shift)
  anch_id <- eval(attr(anchors, "anch_id"))
  shft_id <- eval(attr(anchors, "shft_id"))

  dnames <- format(rel_pos * attr(explist[[1]], "res"), trim = TRUE)
  rawnames <- if (is.null(rownames(anchors))) {
    paste0(anchors[anch_id, 1], ";", anchors[anch_id, 2])
  } else {
    rownames(anchors)
  }
  
  group <- if ("group" %in% names(attributes(anchors))) {
    inverse.rle(attr(anchors, "group"))
  } else {
    NULL
  }

  # Set data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)

  # Decide on engine
  run_engine <- switch(attr(anchors, "type"),
                       "TADs" = "engine_resizing",
                       "ARA"  = "engine_flipping",
                       "engine_default")
  run_engine <- getFromNamespace(run_engine, "GENOVA")

  # Loop over experiments, perform matrix lookup
  results <- lapply(seq_along(explist), function(i) {

    # Lookup matrices
    master <- run_engine(explist[[i]]$MAT, anchors, rel_pos)
    arr <- master[anch_id,,, drop = FALSE]
    mat_mu <- summarise_lookup(arr, outlier_filter, 
                               group = group[anch_id])
    dimnames(mat_mu$mat) <- list(
      rev(dnames), dnames, unique(group[anch_id])
    )[seq_along(dim(mat_mu$mat))]
    if (raw) {
      dimnames(arr) <- list(rawnames[anch_id], rev(dnames), dnames)
      arr <- arr[mat_mu$keep, , , drop = FALSE]
      attr(arr, "group") <- group[anch_id][mat_mu$keep]
    } else {
      arr <- NULL
    }

    # Calculate shifted values
    if (shift > 0) {
      # shifted_arr <- run_engine(explist[[i]]$ICE, shift_anchors, rel_pos)
      shifted_arr <- master[shft_id,,, drop = FALSE]
      shifted_mu <- summarise_lookup(shifted_arr,
                                     outlier_filter,
                                     keep = mat_mu$keep,
                                     group = group[shft_id]
      )$mat
      dimnames(shifted_mu) <- dimnames(mat_mu$mat)
      if (raw) {
        shifted_arr <- shifted_arr[mat_mu$keep, , , drop = FALSE]
        dimnames(shifted_arr) <- dimnames(arr)
        attr(shifted_arr, "group") <- group[shft_id][mat_mu$keep]
      } else {
        shifted_arr <- NULL
      }
    } else {
      shifted_arr <- NULL
      shifted_mu <- NULL
    }

    # Calculate observed over expected
    obsexp <- if (!is.null(shifted_mu)) {
      mat_mu$mat / median(shifted_mu, na.rm = TRUE)
    } else {
      NULL
    }

    results <- list(
      obsexp = obsexp,
      signal = mat_mu$mat,
      signal_raw = arr,
      shifted = shifted_mu,
      shifted_raw = shifted_arr
    )

    null_results <- vapply(results, is.null, logical(1))

    return(
      results[!null_results]
    )
  })

  # Format experiment names
  names(results) <- if (is.null(names(explist))) {
    vapply(explist, function(exp) {
      attr(exp, "samplename")
    }, character(1L))
  } else {
    names(explist)
  }

  merge_res <- names(results[[1]])
  results <- lapply(merge_res, function(res) {
    resname <- res
    objects <- lapply(results, `[[`, res)
    # Return array if rawdata
    if (resname %in% c("shifted_raw", "signal_raw")) {
      return(objects)
    }
    # Simplify matrix to array
    newobject <- do.call(c, objects)
    dim(newobject) <- c(dim(objects[[1]]),
                        length(results))
    dimnames(newobject) <- c(dimnames(objects[[1]]),
                             list(names(results)))
    newobject
  })
  names(results) <- merge_res

  results
}

# Internals --------------------------------------------------------

#' Engine for repeated matrix lookups
#'
#' It looks up square parts of equal size in the Hi-C matrix of a single experiment.
#'
#' @param MAT The Hi-C matrix slot of a GENOVA experiment object.
#' @inheritParams rep_mat_lookup
#'
#' @return A threedimensional \code{array} wherein the first dimension is
#'   parallel to the rows in \code{anchors}, and the second and third dimensions
#'   are of length \code{m}.
#'
#' @seealso \code{\link[GENOVA]{rep_mat_lookup}}
#'
#' @keywords internal
engine_default <- function(MAT, anchors, rel_pos) {
  class(anchors) <- "matrix"

  # Setup indices
  grid <- expand.grid(rel_pos, rel_pos)
  idx <- expand.grid(seq_len(nrow(anchors)), seq_len(nrow(grid)))
  idx <- list(grid[idx$Var2, 1] + anchors[idx$Var1, 1],
              grid[idx$Var2, 2] + anchors[idx$Var1, 2])

  # Get data
  mats <- MAT[idx, V3]
  dim(mats) <- c(nrow(anchors), length(rel_pos)[c(1, 1)])
  mats
}

#' @keywords internal
engine_flipping <- function(MAT, anchors, rel_pos) {

  mats <- engine_default(MAT = MAT, anchors = anchors, rel_pos = rel_pos)

  # Split off reverse orientation
  dir <- inverse.rle(attr(anchors, "dir"))
  dir <- which(dir == "-")

  seq <- seq_along(rel_pos)
  mats[dir, seq, seq] <- aperm(mats[dir, rev(seq), rev(seq)], c(1,3,2))
  attr(mats, "dir") <- attr(anchors, "dir")
  mats
}

#' Engine for repeated matrix lookup and resizing
#'
#' It looks up square parts of unequal size around the diagonal of a single Hi-C
#' experiment and resizes these to a single size.
#'
#' @param ICE The Hi-C matrix slot of a GENOVA experiment object.
#' @param anchors A \code{matrix} with two columns containing pre-computed
#'   indices between which a square is looked up.
#' @param rel_pos A \code{integer} sequence ranging from \code{[0-n]} wherein
#'   \code{n} is the target size to which individual squares are resized.
#'
#'   \code{anchors} specify two points on the diagonal between which the region
#'   should be looked up.
#'
#' @return A threedimensional \code{array} wherein the first dimension is
#'   parallel to the rows in \code{anchors}, and the second and third dimensions
#'   are of length \code{m}.
#'
#' @seealso \code{\link[GENOVA]{rep_mat_lookup}}
#'
#' @keywords internal
engine_resizing <- function(MAT, anchors, rel_pos) {
  template <- matrix(0, tail(rel_pos, 1), tail(rel_pos, 1))
  coords <- as.matrix(expand.grid(rel_pos, rel_pos))
  max_pos <- max(rel_pos)

  # Array transpose to match rep_mat_lookup expectation
  base::aperm(vapply(seq_len(nrow(anchors)), function(j) {
    # Prep parameters
    i <- seq.int(anchors[j, 1], anchors[j, 2])
    len <- length(i)
    seqs <- as.double(seq.int(len))
    min <- i[1] - 1

    # Lookup matrix
    mat <- MAT[CJ(V1 = i, V2 = i),
               dt_matrix(V3, V1, V2, len, min),
               nomatch = NULL]

    # Interpolate coordinates
    xy <- (coords - 1) * (len - 1) / (max_pos - 1) + 1
    x <- xy[, 1]
    y <- xy[, 2]

    # Match to current coordinates
    dx <- x - {flx <- floor(x)}
    dy <- y - {fly <- floor(y)}

    # Mellow out the extremes
    dx[flx == len] <- 1
    dy[fly == len] <- 1
    flx[flx == len] <- len - 1
    fly[fly == len] <- len - 1

    # Interpolate values
    template[coords] <-
      mat[cbind(flx, fly)] * (1 - dx) * (1 - dy) +
      mat[cbind(flx + 1, fly)] * dx * (1 - dy) +
      mat[cbind(flx, fly + 1)] * (1 - dx) * dy +
      mat[cbind(flx + 1, fly + 1)] * dx * dy

    template
  }, template), c(3, 1, 2))
}

#' Summarise the results of a repeated matrix lookup
#'
#' This function either filters out all the array slices that are solely
#' composed of \code{NA}s or those specified. It sets all other \code{NA}s to
#' \code{0} and then takes the mean along the first dimension of a
#' threedimensional array.
#'
#' @inheritParams APA
#' @param array A threedimensional \code{array}.
#' @param keep A \code{logical} vector of length \code{dim(array)[1]} indicated
#'   which slices to keep. If \code{NULL}, independently filters out
#'   all-\code{NA} slices.
#'
#' @return A \code{list} of length 2 with a \code{matrix} giving mean values of
#'   \code{array} along the first dimension and the \code{keep} as used in the
#'   calculation of the means.
#'
#' @seealso \code{\link[GENOVA]{rep_mat_lookup}}
#'
#' @keywords internal
summarise_lookup <- function(array, outlier_filter = c(0, 1), 
                             keep = NULL, group = NULL) {
  if (is.null(keep)) {
    keep <- rowSums(array, na.rm = TRUE) != 0
  }
  
  if (is.null(group)) {
    group <- rep.int(1, dim(array)[1])
  }
  group <- group[keep]
  grouplvls <- unique(group)

  # Remove all-NA slices, set other NAs to 0
  array <- array[keep, , , drop = FALSE]
  array[is.na(array)] <- 0

  # Do outlier filtering
  if (!identical(outlier_filter, c(0, 1))) {
    outlier_filter <- sort(outlier_filter)
    thres <- quantile(array, outlier_filter)
    if (outlier_filter[2] != 1) {
      array <- pmin(array, thres[2])
    }
    if (outlier_filter[1] != 0) {
      array <- pmax(array, thres[1])
    }
  }

  # Summarise
  mat <- vapply(grouplvls, function(lvl) {
    colMeans(array[group == lvl, , , drop = FALSE])
  }, matrix(NA_real_, dim(array)[2], dim(array)[3]))
  
  # Drop last dimension if empty
  if (tail(dim(mat), 1) == 1L) {
    dim(mat) <- head(dim(mat), -1)
  }

  return(list(mat = mat, keep = keep))
}
