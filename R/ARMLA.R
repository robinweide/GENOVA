# Note: ARMLA is the abbreviation for Aggregate Repeated Matrix Lookup Analysis.

# User facing -------------------------------------------------------------

#' Aggregrate Peak Analysis
#'
#' Perform multiple matrix lookup in Hi-C matrices for a twodimensional set of
#' regions, for example loops.
#'
#' @param explist Either a single GENOVA experiment or a list of GENOVA
#'   experiments.
#' @param bedpe A \code{data.frame} with 6 columns in BEDPE format containing
#'   the regions to be anchored: chrom1/start1/end1/chrom2/start2/end2.
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum distances in basepairs between anchorpoints.
#' @param size_bin The size of the lookup regions in bins (i.e. a score of 21
#'   yields an output with 10 Hi-C bins both up- and downstream of the anchor).
#' @param size_bp Alternative parametrisation for the lookup regions, expressed
#'   in basepairs. Not used when the argument \code{size_bin} is set.
#' @param outlier_filter A \code{numeric} of length 2 between \code{[0-1]}
#'   indicating quantiles of data to be used as thresholds. Data outside these
#'   thresholds are set to the nearest threshold. Setting this to \code{c(0, 1)}
#'   performs no outlier correction.
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses these instead.
#' @param raw A \code{logical} of length 1: should the raw array underlying the
#'   summary matrices be returned in the output?
#'
#' @return An \code{APA_discovery} object with the APA results.
#'
#' @seealso The \code{\link[GENOVA]{discovery}} class.
#'   \code{\link[GENOVA]{visualise}} for visualisation of the results.
#'
#' @export
#'
#' @examples
#' # Typical usage: APA for loops
#' apa <- APA(list(WT = WT_10kb, KO = KO_10kb), bedpe = WT_loops)
#'
#' # Alternative usage with pre-calculated anchors
#' anchors <- anchors_APA(WT_10kb$ABS, WT_10kb$RES,
#'   bedpe = WT_loops
#' )
#' apa <- APA(list(WT = WT_10kb, KO = KO_10kb), anchors = anchors)
#'
#' # Visualising results
#' visualise(apa)
APA <- function(explist, bedpe,
                dist_thres = NULL,
                size_bin = 21, size_bp = NULL,
                outlier_filter = c(0, 0.995),
                anchors = NULL, raw = TRUE) {
  # Verify experiment compatability
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- attr(explist[[1]], "res")
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)

  if (is.null(dist_thres)) {
    dist_thres <- c((diff(range(rel_pos)) + 3) * res, Inf)
  }

  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_APA(
      explist[[1]]$IDX,
      res,
      bedpe,
      dist_thres
    )
  }

  results <- rep_mat_lookup(explist, anchors,
    rel_pos = rel_pos,
    shift = 0, outlier_filter = outlier_filter,
    raw = raw
  )

  structure(results, class = c("APA_discovery", "ARMLA_discovery"),
            resolution = res, package = "GENOVA")
}

#' Paired-end spatial chromatin analysis
#'
#' Performs an all-to-all Hi-C contact analysis for specified regions.
#'
#' @inheritParams APA
#' @param bed A \code{data.frame} with 3 columns in BED format, containing the
#'   regions to anchor in pairwise manner.
#' @param shift An \code{integer} of length 1 indicating how many basepairs the
#'   anchors should be shifted. Essentially performs circular permutation of
#'   \code{size} for a reasonable estimate of background. The argument is
#'   ignored when \code{shift <= 0}.
#' @param min_compare An \code{integer} of length 1 indicating the minimum
#'   number of pairwise interactions on a chromosome to consider.
#'
#' @return A \code{PESCAn_discovery} object with the PE-SCAn results.
#'
#' @seealso The \code{\link[GENOVA]{discovery}} class.
#'   \code{\link[GENOVA]{visualise}} for visualisation of the results.
#'
#' @export
#'
#' @examples
#' # Typical usage: PESCAn for super enhancers using a 1 MB
#' # circular permutation on a pair of experiments.
#' pescan <- PESCAn(
#'   explist = list(WT_40kb, KO_40kb),
#'   bed = super_enhancers,
#'   shift = 1e6
#' )
#'
#' # Alternative usage with pre-calculated anchors and no permutation
#' anchors <- anchors_PESCAn(WT_40kb$ABS, WT_40kb$RES,
#'   genes_tss,
#'   dist_thres = c(5e6, 15e6)
#' )
#' pescan <- PESCAn(
#'   explist = list(WT_40kb),
#'   anchors = anchors,
#'   shift = 0
#' )
#'
#' # Visualising PE-SCAns
#' autoplot(pescan)
PESCAn <- function(explist, bed, shift = 1e6L,
                   dist_thres = c(5e6L, Inf),
                   size_bin = NULL, size_bp = 4e5,
                   outlier_filter = c(0, 1),
                   min_compare = 10,
                   anchors = NULL, raw = FALSE) {
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- attr(explist[[1]], "res")
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)
  shift <- round(shift / res)

  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_PESCAn(
      explist[[1]]$IDX, attr(explist[[1]], "res"),
      bed, dist_thres,
      min_compare = min_compare
    )
  }

  results <- rep_mat_lookup(explist, anchors,
    rel_pos = rel_pos,
    shift = shift,
    outlier_filter = outlier_filter,
    raw = raw
  )

  structure(results, class = c("PESCAn_discovery", "ARMLA_discovery"),
            resolution = res, package = "GENOVA")
}

#' Aggregate TAD analysis
#'
#' Extracts Hi-C matrices around TADs, resizes these to a uniform size and
#' averages the results for all TADs.
#'
#' @inheritParams PESCAn
#' @param bed A \code{data.frame} with 3 columns in BED format, containing the
#'   TAD boundary positions per row.
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum sizes of TADs to include in basepairs.
#' @param size A code \code{integer} vector of length 1 noting the dimensions of
#'   the output.
#'
#' @return An \code{ATA_discovery} object with the ATA results.
#'
#' @seealso The \code{\link[GENOVA]{discovery}} class.
#'   \code{\link[GENOVA]{visualise}} for visualisation of the results.
#' @export
ATA <- function(explist, bed,
                dist_thres = c(225000, Inf),
                size = 100L,
                outlier_filter = c(0, 1),
                anchors = NULL, raw = TRUE) {
  # Verify experiment compatibility
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- attr(explist[[1]], "res")
  rel_pos <- seq.int(size)

  if (is.null(anchors)) {
    anchors <- anchors_ATA(
      explist[[1]]$IDX,
      bed,
      dist_thres
    )
  }
  pad <- if ("padding" %in% names(attributes(anchors))) {
    attr(anchors, "padding")
  } else {
    "unknown"
  }

  results <- rep_mat_lookup(explist, anchors, rel_pos = rel_pos,
                            shift = 0, outlier_filter = outlier_filter,
                            raw = raw)

  structure(results, class = c("ATA_discovery", "ARMLA_discovery"),
            package = "GENOVA", resolution = res, padding = pad)
}


#' Aggregate Region Analysis
#'
#' Extracts Hi-C matrices centered around regions and averages the results for
#' all regions.
#'
#' @inheritParams PESCAn
#' @param bed A \code{data.frame} with 3 columns in BED format, containing the
#'   regions to anchor in pairwise manner. Entries wherein the second column is
#'   larger than the third column are considered in the reverse direction.
#'
#' @return An \code{ARA_discovery} object with the ARA results.
#'
#' @details By default, \code{ARA} also calculates the results for shifted
#'   anchors and normalises the \code{"obsexp"} slot by off-diagonal bands.
#'
#'   The '\code{bed}' argument can take in oriented entries, wherein entries
#'   with a start site larger than the end site are considered to be in the
#'   reverse direction. The reverse sites are flipped during analysis, so the
#'   orientation is the same as in the forward sites.
#'
#' @seealso The \code{\link[GENOVA]{discovery}} class.
#'   \code{\link[GENOVA]{visualise}} for visualisation of the results.
#'
#' @export
#'
#' @examples
#' # Typical usage
#' ara <- ARA(list(WT_20kb, KO_20kb), ctcf_sites)
#'
#' # Alternative usage with pre-calculated anchors
#' anchors <- anchors_ARA(WT_20kb$ABS, ctcf_sites)
#' ara <- ARA(list(WT_20kb, KO_20kb))
#'
#' # Visualisation
#' visualise(ara)
ARA <- function(explist, bed, shift = 1e6,
                size_bin = 21, size_bp = NULL,
                outlier_filter = c(0, 1),
                anchors = NULL, raw = FALSE) {
  # Verify experiment compatability
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- attr(explist[[1]], "res")
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)
  shift <- round(shift / res)

  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_ARA(
      explist[[1]]$IDX,
      bed
    )
  }

  results <- rep_mat_lookup(explist, anchors,
                            rel_pos = rel_pos,
                            shift = shift,
                            outlier_filter = outlier_filter,
                            raw = raw
  )
  results$signal <- results$signal + aperm(results$signal, c(2, 1, 3))
  if ("shifted" %in% names(results)) {
    results$shifted <- results$shifted + aperm(results$shifted, c(2, 1, 3))
  }
  if (all("obsexp" %in% names(results))) {
    # Expected per band from diagonal
    obs <- results$signal
    exp <- results$shifted
    rows <- row(obs[,,1]) ; cols <- col(obs[,,1])
    i <- cols - rows
    i <- i - min(i) + 1
    nexp <- apply(exp, 3, function(x) {
      mu <- vapply(split(x, i), function(x){.Internal(mean(x))}, numeric(1))
      mu[i]
    })
    dim(nexp) <- dim(obs)
    results$obsexp <- obs / nexp
  }
  structure(results, class = c("ARA_discovery", "ARMLA_discovery"),
            resolution = res, package = "GENOVA")
}

# Internals ---------------------------------------------------------------

#' Check compatability of a list of GENOVA experiments
#'
#' Checks if the indices (ABS slot) across experiments are identical.
#'
#' @inheritParams APA
#'
#' @return A \code{list} of GENOVA experiment(s).
#'
#' @keywords internal
check_compat_exp <- function(explist) {
  if (!is.list(explist)) {
    stop("Supply either a GENOVA experiment or list of GENOVA experiments",
      call. = FALSE
    )
  }

  # Re-list of only one experiment was given
  if (any(c("MAT", "IDX") %in% names(explist))) {
    explist <- list(explist)
  }

  # Test equality of experiments in list
  if (length(explist) > 1) {
    equal <- vapply(seq_along(explist)[-1], function(i) {
      all.equal(explist[[1]]$IDX, explist[[i]]$IDX)
    }, logical(1))

    if (any(!equal)) {
      stop(paste(
        "Indices of experiment(s)",
        paste(which(!equal) + 1, collapse = " & "),
        "are not equal to indices of experiment 1"
      ), call. = FALSE)
    }
  }

  return(explist)
}

#' Parse size and resolution to relative positions
#'
#' @param res The RES slot of a GENOVA Hi-C experiment.
#' @inheritParams APA
#'
#' @return A \code{integer} vector of relative positions.
#' @keywords internal
parse_rel_pos <- function(res, size_bin, size_bp) {
  if (is.null(size_bin) & !(is.null(size_bp))) {
    # Translate size to relative positions
    rel_pos <- seq_len(((size_bp / res) * 2 + 1))
    # Center size around 0
    rel_pos <- floor(rel_pos - median(rel_pos))
  } else if (!is.null(size_bin)) {
    rel_pos <- seq.int(size_bin)
    rel_pos <- floor(rel_pos - median(rel_pos))
  } else {
    stop("Supply either 'size_bin' or 'size_bp'",
      call. = FALSE
    )
  }
  return(rel_pos)
}
