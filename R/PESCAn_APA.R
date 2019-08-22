# User facing -------------------------------------------------------------

#' Aggregrate Peak Analysis
#'
#' @param explist Either a single GENOVA experiment or a list of GENOVA
#'   experiments.
#' @param bedpe A \code{data.frame} with 6 columns in BEDPE format containing the
#'   regions to be anchored: chrom1/start1/end1/chrom2/start2/end2.
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
#' @param raw A \code{logical} of length 1: should the raw array underlying
#'   the summary matrices be returned in the output?
#'
#' @return A \code{list} of the same length as \code{explist} wherein list
#'   elements contain the results of the APA per experiment.
#' @export
APA <- function(explist, bedpe,
                 dist_thres = NULL,
                 size_bin = 21, size_bp = NULL,
                 outlier_filter = c(0, 1),
                 anchors = NULL, raw = TRUE
) {
  # Verify experiment compatability
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- explist[[1]]$RES
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)

  if (is.null(dist_thres)) {
    dist_thres <- c((diff(range(rel_pos)) + 3) * res)
  }

  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_APA(
      explist[[1]]$ABS,
      explist[[1]]$RES,
      bedpe,
      dist_thres
    )
  }

  results <- rep_mat_lookup(explist, anchors, rel_pos = rel_pos,
                            shift = 0, outlier_filter = outlier_filter,
                            raw = raw)

  structure(results, class = "APA_results")
}

#' Paired-end spatial chromatin analysis
#'
#' Performs an all-to-all Hi-C contact analysis for specified regions.
#'
#' @inheritParams APA2
#' @param bed A \code{data.frame} with 3 columns in BED format, containing the
#'   regions to anchor in pairwise manner.
#' @param shift An \code{integer} of length 1 indicating how many basepairs the
#'   anchors should be shifted. Essentially performs circular permutation of
#'   \code{size} for a reasonable estimate of background. The argument is
#'   ignored when \code{shift <= 0}.
#' @param min_compare An \code{integer} of length 1 indicating the minimum
#'   number of pairwise interactions on a chromosome to consider.
#'
#' @return A \code{list} of the same length as \code{explist} wherein list
#'   elements contain the results of the PE-SCAn per experiment.
#'
#' @export
#'
#' @examples
#' # Typical usage: PESCAn for super enhancers using a 1 MB
#' # circular permutation on a pair of experiments.
#' pescan <- PESCAn2(explist = list(WT_40kb, KO_40kb),
#'                   bed = super_enhancers,
#'                   shift = 1e6)
#'
#' # Alternative usage with pre-calculated anchors and no permutation
#' anchors <- anchors_PESCAn(WT_40kb$ABS, WT_40kb$RES,
#'                           genes_tss,
#'                           dist_thres = c(5e6, 15e6))
#' pescan <- PESCAn2(explist = list(WT_40kb),
#'                   anchors = anchors,
#'                   shift = 0)
PESCAn <- function(explist, bed, shift = 1e6L,
                    dist_thres = c(5e6L, Inf),
                    size_bin = NULL, size_bp = 4e5,
                    outlier_filter = c(0, 1),
                    min_compare = 10,
                    anchors = NULL, raw = FALSE
) {
  explist <- check_compat_exp(explist)

  # Initialise parameters
  res <- explist[[1]]$RES
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)
  shift <- round(shift / res)

  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_PESCAn(
      explist[[1]]$ABS, explist[[1]]$RES,
      bed, dist_thres, min_compare = min_compare
    )
  }

  results <- rep_mat_lookup(explist, anchors, rel_pos = rel_pos,
                            shift = shift,
                            outlier_filter = outlier_filter,
                            raw = raw)

  structure(results, class = "PESCAn_results")
}

# Internals ---------------------------------------------------------------

#' Check compatability of a list of GENOVA experiments
#'
#' Checks if the indices (ABS slot) across experiments are identical.
#'
#' @inheritParams APA2
#'
#' @return A \code{list} of GENOVA experiment(s).
#'
#' @keywords internal
check_compat_exp <- function(explist) {
  if (!is.list(explist)) {
    stop("Supply either a GENOVA experiment or list of GENOVA experiments",
         call. = FALSE)
  }

  # Re-list of only one experiment was given
  if (any(c("ICE","ABS","RES") %in% names(explist))) {
    explist <- list(explist)
  }

  # Test equality of experiments in list
  if (length(explist) > 1) {
    equal <- vapply(seq_along(explist)[-1], function(i){
      all.equal(explist[[1]]$ABS, explist[[i]]$ABS)
    }, logical(1))

    if(any(!equal)){
      stop(paste(
        "Indices of experiment(s)",
        paste(which(!equal) + 1, collapse = " & "),
        "are not equal to indices of experiment 1"
      ), call. = FALSE
      )
    }
  }

  return(explist)
}

#' Parse size and resolution to relative positions
#'
#' @param res The RES slot of a GENOVA Hi-C experiment.
#' @inheritParams APA2
#'
#' @return A \code{integer} vector of relative positions.
#' @keywords internal
parse_rel_pos <- function(res, size_bin, size_bp) {
  if (is.null(size_bin) & !(is.null(size_bp))) {
    # Translate size to relative positions
    rel_pos  <- seq_len(((size_bp / res)*2 + 1))
    # Center size around 0
    rel_pos  <- floor(rel_pos - median(rel_pos))
  } else if (!is.null(size_bin)) {

    rel_pos <- seq.int(size_bin)
    rel_pos <- floor(rel_pos - median(rel_pos))

  } else {
    stop("Supply either 'size_bin' or 'size_bp'",
         call. = FALSE)
  }
  return(rel_pos)
}
