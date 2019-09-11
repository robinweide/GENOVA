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
#'   anchor indices.
#' @param rel_pos An \code{integer} vector indicating relative positions in
#'   bins, indicating ranges around anchors to lookup.
#'
#' @return A \code{list} of the same length as \code{explist} wherein list
#'   elements contain the results of the repeated lookup per experiment.
#' @export
#'
#' @seealso \code{\link[GENOVA]{APA}} and \code{\link[GENOVA]{PESCAn}}
rep_mat_lookup <- function(
                           explist, anchors, rel_pos, shift = 0, outlier_filter = c(0, 1), raw = FALSE) {
  anchors <- anchors_filter_oob(
    explist[[1]]$ABS,
    anchors, rel_pos
  )

  if (shift > 0) {
    shift_anchors <- anchors_shift(
      explist[[1]]$ABS,
      anchors, rel_pos, shift
    )
  }

  dnames <- format(rel_pos * explist[[1]]$RES, trim = TRUE)

  # Set data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)

  lookup_fun <- if (!is.null(attr(anchors, "type")) &
                    attr(anchors, "type") == "TADs") {
    getFromNamespace("lookup_resizer", "GENOVA")
  } else {
    getFromNamespace("matrix_lookup", "GENOVA")
  }

  # Loop over experiments, calculate PESCAn
  results <- lapply(seq_along(explist), function(i) {

    # Calculate true values
    arr <- lookup_fun(explist[[i]]$ICE, anchors, rel_pos)
    mat_mu <- summarise_lookup(arr, outlier_filter)
    dimnames(mat_mu$mat) <- list(rev(dnames), dnames)
    if (raw) {
      dimnames(arr) <- list(
        paste0(anchors[, 1], ",", anchors[, 2]),
        rev(dnames), dnames
      )
      arr <- arr[mat_mu$keep, , ]
    } else {
      arr <- NULL
    }

    # Calculate shifted values
    if (shift > 0) {
      shifted_arr <- matrix_lookup(explist[[i]]$ICE, shift_anchors, rel_pos)
      shifted_mu <- summarise_lookup(shifted_arr,
                                     outlier_filter,
                                     keep = mat_mu$keep
      )$mat
      dimnames(shifted_mu) <- dimnames(mat_mu$mat)
      if (raw) {
        shifted_arr <- shifted_arr[mat_mu$keep, , ]
        dimnames(shifted_arr) <- dimnames(arr)
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
      exp$NAME
    }, character(1L))
  } else {
    names(explist)
  }

  merge_res <- names(results[[1]])
  results <- lapply(merge_res, function(res) {
    objects <- lapply(results, `[[`, res)
    # Return array if rawdata
    if (inherits(objects[[1]], "array")) {
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
#' @param ICE The Hi-C matrix slot of a GENOVA experiment object.
#' @inheritParams rep_mat_lookup
#'
#' @return A threedimensional \code{array} wherein the first dimension is
#'   parallel to the rows in \code{anchors}, and the second and third dimensions
#'   are of length \code{m}.
#'
#' @seealso \code{\link[GENOVA]{rep_mat_lookup}}
#'
#' @keywords internal
matrix_lookup <- function(ICE, anchors, rel_pos) {
  # Setup a grid of relative positions
  grid <- expand.grid(rel_pos, rel_pos)
  idx <- seq_len(nrow(grid))

  # Loop over grid entries, get scores at all anchor positions
  mats <- vapply(idx, function(i) {
    ICE[list(
      grid[i, 1] + anchors[, 1],
      grid[i, 2] + anchors[, 2]
    )]$V3
  }, as.numeric(anchors[, 1]))

  # Cast results in array
  array(mats,
    dim = c(
      nrow(anchors),
      length(rel_pos)[c(1, 1)]
    )
  )
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
lookup_resizer <- function(ICE, anchors, rel_pos) {
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

    # Setup indices
    idx <- data.table(V1 = rep(i, each = len),
                      V2 = rep.int(i, len),
                      key = c("V1", "V2"))

    # Lookup matrix
    mat <- ICE[idx, dt_matrix(V3, V1, V2, len, min),
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
summarise_lookup <- function(
                             array, outlier_filter = c(0, 1), keep = NULL) {
  if (is.null(keep)) {
    keep <- apply(array, 1, function(slice) {
      !all(is.na(slice))
    })
  }

  # Remove all-NA slices, set other NAs to 0
  array <- array[keep, , ]
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
  mat <- apply(array, 2:3, function(x) .Internal(mean(x)))
  return(list(mat = mat, keep = keep))
}
