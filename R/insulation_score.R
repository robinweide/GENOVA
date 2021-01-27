#' Calculate insulation scores
#'
#' Insulation scores are calculated by sliding a square along the diagonal of
#' the Hi-C matrix and averaging the values in that square.
#'
#' @param explist Either a single GENOVA \code{contacts} object or list of
#'   GENOVA \code{contact} objects.
#' @param window An \code{integer} of length 1 for the size of the sliding
#'   square.
#' @param norm_to A \code{character} of length 1, either \code{"chromosome"} or
#'   \code{"genome"} noting whether normalisation should occur per chromosome or
#'   for the genome as a whole. Can be set to \code{"none"} to skip
#'   normalisation.
#' @param norm_fun A \code{function} that takes a numeric vector as input and
#'   returns normalised values of equal length. Defaults to calculating the
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} value over the median.
#'
#' @details Typically, when the sliding square passes within a TAD, it yields
#'   high insulation scores whereas in between TADs it yields low insulation
#'   scores. Hence, TAD boundaries can be identified as local minima in the
#'   insulation score.
#'
#'   To follow the Crane \emph{et al}. (2015) strategy for insulation scores,
#'   use 10kb resolution \code{contacts} objects, set the '\code{window}'
#'   argument to 50, set the '\code{norm_to}' argument to \code{"chromosome"}
#'   and the '\code{norm_fun}' argument to \code{log2overmean}.
#'
#' @return An \code{IS_discovery} object containing the following slot: \describe{
#'   \item{insula_score}{A \code{data.table} with genomic locations and
#'   insulation scores for each element in the '\code{explist}' argument.}}
#' @export
#' 
#' @section Resolution recommendation: 10kb-40kb
#'
#' @seealso For calling TADs from insulation scores, see
#'   \code{\link[GENOVA]{call_TAD_insulation}}. For plotting a heatmap of
#'   insulation over genomic locations, see
#'   \code{\link[GENOVA]{tornado_insulation}}.
#'
#' @examples
#' \dontrun{
#' iscore <- insulation_score(list(WT_20kb, KO_20kb), window = 20)
#' }
insulation_score <- function(explist, window = 30,
                             norm_to = c("chromosome", "genome", "none"),
                             norm_fun = log2overmedian) {
  norm_to <- match.arg(norm_to)
  if (norm_to != "none") {
    if (!is.function(norm_fun)) {
      stop("Please provide a function as the `norm_fun` argument.",
           call. = FALSE)
    }
  }
  
  explist <- check_compat_exp(explist)
  expnames <- expnames(explist)
  nexp <- length(expnames)
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Prepare bins
  ext_window <- 2 * window
  bins <- explist[[1]]$IDX$V4
  outer_bins <- CJ(V1 = bins, V2 = (ext_window - 1):(window))
  outer_bins <- outer_bins[, grp := V1 + V2 - window]
  inner_bins <- outer_bins[, list(V1, V2 = V2 - window, grp)]
  bins <- CJ(V1 = bins, V2 = 0:(ext_window - 1))
  
  
  raw <- lapply(explist, function(xp) {
    # Straighten diagonal
    mat <- xp$MAT[V1 > V2 - (ext_window)]
    mat[, V2 := V2 - V1]
    setkeyv(mat, c("V1", 'V2'))
    
    # Pre-calculate cumulative sums
    mat <- mat[bins][is.na(V3), V3 := 0]
    mat[, V3 := cumsum(V3), by = "V1"]
    
    # Calculate outer and inner partial sums
    outer_sum <- mat[outer_bins, nomatch = NULL][, .(V2 = sum(V3)), keyby="grp"]
    inner_sum <- mat[inner_bins, nomatch = NULL][, .(V3 = sum(V3)), keyby="grp"]
    
    # Calculate overall sum by subtracting the inners from the outers
    sums <- outer_sum[inner_sum, list(grp, V2 = V2 - V3)]
    sums
  })
  
  # Combine samples
  insula <- raw[[1]]
  for(i in tail(seq_along(raw), -1)) {
    insula <- merge.data.table(insula, raw[[i]], all = TRUE, by = "grp")
  }
  
  # Match to orginal bins
  insula <- insula[explist[[1]]$IDX, on = c("grp" = "V4")]
  setnames(insula, seq_len(ncol(insula)), 
           new = c("bin", expnames, "chrom", "start", "end"))
  setcolorder(insula, c((nexp + 2):ncol(insula), 1:(nexp + 1)))
  
  # Delete faulty bins
  antibin <- insula[, c(head(bin, 1) + c(0, seq_len(window - 1)), 
                        tail(bin, 1) - c(0, seq_len(window - 1))),
                    by = "chrom"]
  insula <- insula[!antibin, on = c("bin" = "V1")]
  
  # Do normalisation
  if (norm_to == "chromosome") {
    for (i in expnames) {
      xp <- as.symbol(i)
      insula[, as.character(xp) := norm_fun(eval(xp)), by = chrom]
    }
  } else if (norm_to == "genome") {
    for (i in expnames) {
      xp <- as.symbol(i)
      insula[, as.character(xp) := norm_fun(eval(xp))]
    }
  } else {
    # Take window averages instead of sums
    for (i in expnames) {
      xp <- as.symbol(i)
      insula[, as.character(xp) := eval(xp) / (window ^ 2)]
    }
  }
  setkeyv(insula, "bin")
  
  # Format output
  structure(list(insula_score = as.data.frame(insula)),
            package = "GENOVA",
            colours = vapply(explist, attr, character(1L), "colour"),
            class = c("IS_discovery", "genomescore_discovery", "discovery"),
            resolution = attr(explist[[1]], "resolution"),
            window = window)
}

#' Log 2 value over its median
#'
#' A small convenience function to calculate the log 2 of a value divided by the
#' median of this set of values.
#'
#' @param x A \code{numeric} vector
#'
#' @return A \code{numeric} vector with the calculated values.
#' 
#' @details Ignores \code{NA} values in calculation of the median.
#' @export
#'
#' @examples
#' log2overmedian(runif(100))
log2overmedian <- function(x) {
  log2(x / median(x, na.rm = TRUE))
}

#' Log 2 value over its mean
#'
#' A small convenience function to calculate the log 2 of a value divided by the
#' mean of this set of values.
#'
#' @param x A \code{numeric} vector
#'
#' @return A \code{numeric} vector with the calculated values.
#' 
#' @details Ignores \code{NA} values in calculation of the mean
#' @export
#'
#' @examples
#' log2overmean(runif(100))
log2overmean <- function(x) {
  log2(x / mean(x, na.rm = TRUE))
}
