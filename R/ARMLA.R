# Note: ARMLA is the abbreviation for Aggregate Repeated Matrix Lookup Analysis.

# User facing -------------------------------------------------------------

#' Aggregrate Peak Analysis
#'
#' Perform multiple matrix lookup in Hi-C matrices for a twodimensional set of
#' locations, for example loops.
#'
#' @param explist Either a single GENOVA \code{contacts} object or a list of
#'   GENOVA \code{contacts} objects.
#' @param bedpe A BEDPE-formatted \code{data.frame} with the following 6
#'   columns: \enumerate{ \item A \code{character} giving the
#'   chromosome names of the first coordinate. \item An \code{integer} giving
#'   the start positions of the first coordinate. \item An \code{integer} giving
#'   the end positions of the first coordinate. \item A \code{character} giving
#'   the chromosome names of the second coordinate. \item An \code{integer}
#'   giving the start positions of the second coordinate. \item An
#'   \code{integer} giving the end positions of the second coordinate. }
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum distances in basepairs between anchorpoints.
#' @param size_bin,size_bp The size of the lookup regions in bins (i.e. a score
#'   of 21 yields an output with 10 Hi-C bins both up- and downstream of the
#'   anchor). When \code{NULL} (default), it is internally set to \code{21} when
#'   the \code{size_bp} is also \code{NULL}. \code{size_bp} is an alternative
#'   parametrisation for the lookup regions, expressed in basepairs.
#'   \code{size_bp} is not used when the argument \code{size_bin} is set.
#' @param outlier_filter A \code{numeric} of length 2 between \code{[0-1]}
#'   indicating quantiles of data to be used as thresholds. Data outside these
#'   thresholds are set to the nearest threshold. Setting this to \code{c(0, 1)}
#'   performs no outlier correction.
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses this argument instead. See 
#'   \code{\link[GENOVA:anchors]{anchors_APA()}}.
#' @param raw A \code{logical} of length 1: should the raw array underlying the
#'   summary matrices be returned in the output? Should be \code{TRUE} if the
#'   intention is to use the \code{quantify} function.
#'
#' @return An \code{APA_discovery} object containing the following slots:
#'   \describe{ \item{signal}{An \code{array} with the dimensions
#'   \code{size_bin} x \code{size_bin} x \code{length(explist)} containing mean
#'   contact values for bins surrounding the anchors.} \item{signal_raw}{A
#'   \code{list} with \code{length(explist)} elements for each contacts object,
#'   wherein an element is an n x \code{size_bin} x \code{size_bin} array with
#'   contact values for each anchor. 'n' is the number of non-empty valid
#'   anchors.}}
#'
#' @details For each row in the '\code{bedpe}' or '\code{anchors}' argument, an
#'   \code{size_bin} x \code{size_bin} region centered on that location is
#'   retrieved. This data is then summarised by taking the mean for every
#'   element in these matrices across all locations.
#'
#'   The '\code{bedpe}' argument is converted internally to an '\code{anchors}'
#'   object.
#'
#' @seealso The \code{\link[GENOVA]{rep_mat_lookup}} function that performs the
#'   lookup and summary for the \code{APA} function and others. \cr The
#'   \code{\link[GENOVA]{discovery}} class for a general description of
#'   \code{discovery} classes. \cr The \code{\link[GENOVA]{visualise}} function
#'   for visualisation of the results. \cr The \code{\link[GENOVA]{quantify}}
#'   function for quantification of loop strenghts. \cr The
#'   \code{\link[GENOVA]{anchors}} documentation for more information about
#'   anchors.
#'   
#' @family aggregate repeated matrix lookup analyses
#'
#' @section Resolution recommendation: 10kb-20kb
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Typical usage: APA for loops
#' apa <- APA(list(WT = WT_10kb, KO = KO_10kb), bedpe = WT_loops)
#'
#' # Alternative usage with pre-calculated anchors
#' anchors <- anchors_APA(WT_10kb$ABS, resolution(WT_10kb),
#'   bedpe = WT_loops
#' )
#' apa <- APA(list(WT = WT_10kb, KO = KO_10kb), anchors = anchors)
#'
#' # Visualising results
#' visualise(apa)
#' }
APA <- function(explist, bedpe,
                dist_thres = NULL,
                size_bin = NULL, size_bp = NULL,
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

  structure(results, class = c("APA_discovery", "ARMLA_discovery", "discovery"),
            resolution = res, package = "GENOVA")
}

#' Paired-end spatial chromatin analysis
#'
#' Performs an all-to-all Hi-C contact analysis for specified regions.
#'
#' @inheritParams APA
#' @param bed A BED-formatted \code{data.frame} with the following 3 columns:
#'   \enumerate{ \item A \code{character} giving the chromosome names. \item An
#'   \code{integer} with start positions. \item An \code{integer} with end
#'   positions. }
#' @param shift An \code{integer} of length 1 indicating how many basepairs the
#'   anchors should be shifted. Essentially performs circular permutation of
#'   \code{size} for a reasonable estimate of background. The argument is
#'   ignored when \code{shift <= 0}.
#' @param min_compare An \code{integer} of length 1 indicating the minimum
#'   number of pairwise interactions on a chromosome to consider.
#' @param size_bin,size_bp The size of the lookup regions in bins (i.e. a score
#'   of 21 yields an output with 10 Hi-C bins both up- and downstream of the
#'   anchor). When \code{NULL} (default), it is internally set to \code{21} when
#'   the \code{size_bp} is also \code{NULL}. \code{size_bp} is an alternative
#'   parametrisation for the lookup regions, expressed in basepairs.
#'   \code{size_bp} is not used when the argument \code{size_bin} is set.
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses this argument instead. See
#'   \code{\link[GENOVA:anchors]{anchors_PESCAn()}}.
#'
#' @return An \code{PESCAn_discovery} object containing the following slots:
#'   \describe{ \item{obsexp}{An \code{array} with the dimensions
#'   \code{size_bin} x \code{size_bin} x \code{length(explist)} containing fold
#'   change values for the signal over the median shifted values.}
#'   \item{signal}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{length(explist)} containing mean contact values for
#'   bins surrounding the anchors.} \item{signal_raw}{A \code{list} with
#'   \code{length(explist)} elements for each contacts object, wherein an
#'   element is an n x \code{size_bin} x \code{size_bin} array with contact
#'   values for each anchor. 'n' is the number of non-empty valid anchors.}
#'   \item{shifted}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{length(explist)} containing mean contact values for
#'   bins that are \code{shift} basepairs away from the anchors.}
#'   \item{shifted_raw}{A \code{list} with \code{length(explist)} elements for
#'   each contacts object, wherein an alemeent is an n x \code{size_bin} x
#'   \code{size_bin} array with contact values for each shifted anchors. 'n' is 
#'   the number of non-empty valid (unshifted) anchors.}}
#'
#' @section Resolution recommendation: 20kb-40kb
#'
#' @seealso The \code{\link[GENOVA]{rep_mat_lookup}} function that performs the
#'   lookup and summary for the \code{PESCAn} function and others. \cr The
#'   \code{\link[GENOVA]{discovery}} class for a general description of
#'   \code{discovery} classes. \cr The \code{\link[GENOVA]{visualise}} function
#'   for visualisation of the results. \cr The \code{\link[GENOVA]{quantify}}
#'   function for quantification of interaction strenghts. \cr The
#'   \code{\link[GENOVA]{anchors}} documentation for more information about
#'   anchors.
#'   
#' @family aggregate repeated matrix lookup analyses
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Typical usage: PESCAn for super enhancers using a 1 MB
#' # circular permutation on a pair of experiments.
#' pescan <- PESCAn(
#'   explist = list(WT_40kb, KO_40kb),
#'   bed = super_enhancers,
#'   shift = 1e6
#' )
#'
#' # Alternative usage with pre-calculated anchors and no permutation
#' anchors <- anchors_PESCAn(WT_40kb$IDX, resolution(WT_40kb),
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
#' visualise(pescan)
#' }
PESCAn <- function(explist, bed, shift = 1e6L,
                   dist_thres = c(5e6L, Inf),
                   size_bin = NULL, size_bp = NULL,
                   outlier_filter = c(0, 1),
                   min_compare = 10,
                   anchors = NULL, raw = TRUE) {
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

  structure(results, 
            class = c("PESCAn_discovery", "ARMLA_discovery", "discovery"),
            resolution = res, package = "GENOVA")
}

#' Aggregate TAD analysis
#'
#' Extracts Hi-C matrices around TADs, resizes these to a uniform size and
#' averages the results for all TADs.
#'
#' @inheritParams PESCAn
#' @param bed A \code{data.frame} with in BED format, containing the
#'   a TAD region per row with the following three columns:
#'   \enumerate{ \item A \code{character} giving the chromosome names. \item An
#'   \code{integer} with start positions. \item An \code{integer} with end
#'   positions. }
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum sizes of TADs to include in basepairs.
#' @param size A code \code{integer} vector of length 1 noting the dimensions of
#'   the output.
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses this argument instead. See
#'   \code{\link[GENOVA:anchors]{anchors_ATA()}}.
#'
#' @return An \code{ATA_discovery} object containing the following slots:
#'   \describe{ \item{signal}{An \code{array} with the dimensions
#'   \code{size} x \code{size} x \code{length(explist)} containing mean
#'   bilinearly interpolated values across TADs.} \item{signal_raw}{A
#'   \code{list} with \code{length(explist)} elements for each contacts object,
#'   wherein an element is an n x \code{size} x \code{size} array with
#'   bilinearly interpolated contact values for each TAD. 'n' is the number of 
#'   non-empty valid TADs.}}
#'
#' @section Resolution recommendation: 10kb-40kb
#'
#' @seealso The \code{\link[GENOVA]{rep_mat_lookup}} function that performs the
#'   lookup and summary for the \code{ATA} function and others. \cr The
#'   \code{\link[GENOVA]{discovery}} class for a general description of
#'   \code{discovery} classes. \cr The \code{\link[GENOVA]{visualise}} function
#'   for visualisation of the results. \cr The \code{\link[GENOVA]{quantify}}
#'   function for quantification of TAD strenghts. \cr The
#'   \code{\link[GENOVA]{anchors}} documentation for more information about
#'   anchors.
#' @family aggregate repeated matrix lookup analyses
#'   
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

  structure(results, 
            class = c("ATA_discovery", "ARMLA_discovery", "discovery"),
            package = "GENOVA", resolution = res, padding = pad)
}


#' Aggregate Region Analysis
#'
#' Extracts Hi-C matrices centered around regions and averages the results for
#' all regions. Can take orientations of regions into account.
#' 
#' @inheritParams PESCAn
#' @param bed A BED-formatted \code{data.frame} with the following 3 columns:
#'   \enumerate{ \item A \code{character} giving the chromosome names. \item An
#'   \code{integer} with start positions. \item An \code{integer} with end
#'   positions. } Note that rows wherein the second column is larger than the
#'   third column are considered to be in the reverse direction.
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses this argument instead. See
#'   \code{\link[GENOVA:anchors]{anchors_ARA()}}.
#' @param size_bin,size_bp The size of the lookup regions in bins (i.e. a score
#'   of 21 yields an output with 10 Hi-C bins both up- and downstream of the
#'   anchor). When \code{NULL} (default), it is internally set to \code{21} when
#'   the \code{size_bp} is also \code{NULL}. \code{size_bp} is an alternative
#'   parametrisation for the lookup regions, expressed in basepairs.
#'   \code{size_bp} is not used when the argument \code{size_bin} is set.
#' 
#' @return An \code{ARA_discovery} object containing the following slots:
#'   \describe{ \item{obsexp}{An \code{array} with the dimensions
#'   \code{size_bin} x \code{size_bin} x \code{length(explist)} containing fold
#'   change values for the signal over the (banded) median, shifted values.}
#'   \item{signal}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{length(explist)} containing mean contact values for
#'   bins surrounding the anchors.} \item{signal_raw}{A \code{list} with
#'   \code{length(explist)} elements for each contacts object, wherein an
#'   element is an n x \code{size_bin} x \code{size_bin} array with contact
#'   values for each anchor. 'n' is the number of non-empty valid anchors.}
#'   \item{shifted}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{length(explist)} containing mean contact values for
#'   bins that are \code{shift} basepairs away from the anchors.}
#'   \item{shifted_raw}{A \code{list} with \code{length(explist)} elements for
#'   each contacts object, wherein an alemeent is an n x \code{size_bin} x
#'   \code{size_bin} array with contact values for each shifted anchors. 'n' is 
#'   the number of non-empty valid (unshifted) anchors.}}
#'
#' @details By default, \code{ARA} also calculates the results for shifted
#'   anchors and normalises the \code{"obsexp"} slot by off-diagonal bands.
#'
#'   The '\code{bed}' argument can take in oriented entries, wherein entries
#'   with a start site larger than the end site are considered to be in the
#'   reverse direction. The reverse sites are flipped during analysis, so the
#'   orientation is the same as in the forward sites.
#'
#' @section Resolution recommendation: 10kb-40kb
#'
#' @seealso The \code{\link[GENOVA]{rep_mat_lookup}} function that performs the
#'   lookup and summary for the \code{ARA} function and others. \cr The
#'   \code{\link[GENOVA]{discovery}} class for a general description of
#'   \code{discovery} classes. \cr The \code{\link[GENOVA]{visualise}} function
#'   for visualisation of the results. \cr The \code{\link[GENOVA]{quantify}}
#'   function for quantification of ARA results. \cr The
#'   \code{\link[GENOVA]{anchors}} documentation for more information about
#'   anchors.
#' @family aggregate repeated matrix lookup analyses
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Typical usage
#' ara <- ARA(list(WT_20kb, KO_20kb), ctcf_sites)
#'
#' # Alternative usage with pre-calculated anchors
#' anchors <- anchors_ARA(WT_20kb$IDX, ctcf_sites)
#' ara <- ARA(list(WT_20kb, KO_20kb))
#'
#' # Visualisation
#' visualise(ara)
#' }
ARA <- function(explist, bed, shift = 1e6,
                size_bin = NULL, size_bp = NULL,
                outlier_filter = c(0, 1),
                anchors = NULL, raw = TRUE) {

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
      mu <- vapply(split(x, i), function(x){mean(x)}, numeric(1))
      mu[i]
    })
    dim(nexp) <- dim(obs)
    results$obsexp <- obs / nexp
  }
  structure(results, 
            class = c("ARA_discovery", "ARMLA_discovery", "discovery"),
            resolution = res, package = "GENOVA")
}


#' Cross spatial chromatin analysis
#'
#' Takes a list of BED-like genomic locations and makes a Hi-C contact analysis
#' for crosswise pairs between elements of that list.
#'
#' @inheritParams PESCAn
#' @param bedlist A named \code{list} of length >= 2 wherein each element is a 
#'   BED-like \code{data.frame} with the following three columns: 
#'   \enumerate{ \item A \code{character} giving the chromosome names. \item An
#'   \code{integer} with start positions. \item An \code{integer} with end
#'   positions. }
#' @param anchors (Optional) A \code{matrix} with two columns containing
#'   pre-computed anchor indices. If this is set, skips calculation of anchor
#'   indices and uses this argument instead. See
#'   \code{\link[GENOVA:anchors]{anchors_CSCAn()}}.
#' @param size_bin,size_bp The size of the lookup regions in bins (i.e. a score
#'   of 21 yields an output with 10 Hi-C bins both up- and downstream of the
#'   anchor). When \code{NULL} (default), it is internally set to \code{21} when
#'   the \code{size_bp} is also \code{NULL}. \code{size_bp} is an alternative
#'   parametrisation for the lookup regions, expressed in basepairs.
#'   \code{size_bp} is not used when the argument \code{size_bin} is set.
#' @param group_direction A \code{logical} of length 1 which when \code{TRUE}
#'   will mirror groups for anchors, where the left anchor location is larger
#'   than the right anchor location. Left and right refer to the \code{bedlist}
#'   elements that generate combinations.
#'
#' @return An \code{CSCAn_discovery} object containing the following slots, 
#'   wherein \code{i} is the number of combinations between the \code{bedlist} 
#'   argument elements:
#'   \describe{ \item{obsexp}{An \code{array} with the dimensions
#'   \code{size_bin} x \code{size_bin} x \code{i} x \code{length(explist)} containing fold
#'   change values for the signal over the median shifted values.}
#'   \item{signal}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{i} x \code{length(explist)} containing mean contact values for
#'   bins surrounding the anchors.} \item{signal_raw}{A \code{list} with
#'   \code{length(explist)} elements for each contacts object, wherein an
#'   element is an n x \code{size_bin} x \code{size_bin} array with contact
#'   values for each anchor. 'n' is the number of non-empty valid anchors.}
#'   \item{shifted}{An \code{array} with the dimensions \code{size_bin} x
#'   \code{size_bin} x \code{i} x \code{length(explist)} containing mean contact values for
#'   bins that are \code{shift} basepairs away from the anchors.}
#'   \item{shifted_raw}{A \code{list} with \code{length(explist)} elements for
#'   each contacts object, wherein an alemeent is an n x \code{size_bin} x
#'   \code{size_bin} array with contact values for each shifted anchors. 'n' is 
#'   the number of non-empty valid (unshifted) anchors.}}
#'   
#' @export
#'
#' @seealso The \code{\link[GENOVA]{rep_mat_lookup}} function that performs the
#'   lookup and summary for the \code{CSCan} function and others. \cr The
#'   \code{\link[GENOVA]{discovery}} class for a general description of
#'   \code{discovery} classes. \cr The \code{\link[GENOVA]{visualise}} function
#'   for visualisation of the results. \cr The \code{\link[GENOVA]{anchors}}
#'   documentation for more information about anchors.
#' @family aggregate repeated matrix lookup analyses
#'
#' @examples
#' \dontrun{
#' cscan <- CSCAn(list(WT_20kb, KO_20kb),
#'                bedlist = list(group_A, group_B))
#'                
#' visualise(cscan)
#' }
CSCAn <- function(explist, bedlist, shift = 1e6L,
                  dist_thres = c(NA, 1e6),
                  size_bin = NULL, size_bp = NULL,
                  outlier_filter = c(0, 1),
                  min_compare = 10,
                  anchors = NULL, raw = TRUE,
                  group_direction = FALSE) {
  explist <- check_compat_exp(explist)
  
  if (!missing(bedlist)) {
    if (length(bedlist) < 2 || !inherits(bedlist, "list")) {
      stop("Less than two 'bedlist' elements found. For self-interaction of a",
           " single BED-like data.frame, see '?PESCAn'.",
           call. = FALSE)
    }
  } else if (missing(anchors)) {
    stop("CSCAn either needs precomputed anchors or the 'bedlist' argument.")
  }

  
  # Initialise parameters
  res <- attr(explist[[1]], "res")
  rel_pos <- parse_rel_pos(res, size_bin, size_bp)
  shift <- round(shift / res)
  
  # Dynamically set lower limit
  if (is.na(dist_thres[1])) {
    dist_thres[1] <- length(rel_pos) * res
  }
  
  # Calculate anchors
  if (is.null(anchors)) {
    anchors <- anchors_CSCAn(
      explist[[1]]$IDX, attr(explist[[1]], "res"),
      bedlist, dist_thres,
      min_compare = min_compare,
      group_direction = group_direction
    )
  }
  
  results <- rep_mat_lookup(explist, anchors,
                            rel_pos = rel_pos,
                            shift = shift,
                            outlier_filter = outlier_filter,
                            raw = raw
  )
  
  structure(results, 
            class = c("CSCAn_discovery", "ARMLA_discovery", "discovery"),
            resolution = res, package = "GENOVA")
}

# Internals ---------------------------------------------------------------

#' Parse size and resolution to relative positions
#'
#' @param res The RES slot of a GENOVA Hi-C experiment.
#' @inheritParams APA
#'
#' @return A \code{integer} vector of relative positions.
#' @keywords internal
#' @noRd
parse_rel_pos <- function(res, size_bin, size_bp) {
  if (!is.null(size_bin) & !is.null(size_bp)) {
    message("Both 'size_bin' and 'size_bp' are set. Continuing with 'size_bin'.")
  }
  if (is.null(size_bin) & !(is.null(size_bp))) {
    # Translate size to relative positions
    rel_pos <- seq_len(((size_bp / res) * 2 + 1))
    # Center size around 0
    rel_pos <- floor(rel_pos - median(rel_pos))
  } else if (!is.null(size_bin)) {
    rel_pos <- seq.int(size_bin)
    rel_pos <- floor(rel_pos - median(rel_pos))
  } else {
    rel_pos <- seq.int(21)
    rel_pos <- floor(rel_pos - median(rel_pos))
  }
  return(rel_pos)
}
