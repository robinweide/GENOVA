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
#'   \code{\link[GENOVA]{heatmap_insulation}}.
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
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  nexp <- length(expnames)
  cols <- vapply(explist, attr, character(1L), "colour")
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Prepare grid
  grid <- CJ(x = seq_len(window), y = seq_len(window))
  grid[, y := y + window]
  grid <- grid - 1 # when added to idx, re-base to 0
  mode(grid$x) <- "integer"
  mode(grid$y) <- "integer"
  grix <- grid$x
  griy <- grid$y
  
  # Prep stuff
  idx <- explist[[1]]$IDX[, head(V4, - window), by = V1]
  setnames(idx, 1:2, c("V1", "V2"))
  
  # Loop over chromosomes so that sorting doesn't take ages
  chroms <- split(idx[, list(V2)], idx$V1)
  insula <- lapply(chroms, function(indices) {
    
    # Setup indices to look up
    looper <- indices[, list(V1 = grix + .BY[[1]], V2 = griy + .BY[[1]]),
                      by = list("id" = V2)]
    setcolorder(looper, c(2,3,1))
    setkey(looper, V1, V2)

    # Grab insulation scores
    ins <- vapply(explist, function(exp) {
      exp$MAT[looper, list(id, V3)][, sum(V3, na.rm = T),  by = id][["V1"]]
    }, numeric(nrow(indices)))
    data.table(ins / nrow(grid))
  })

  # Combine results from all chromosomes
  insula <- rbindlist(insula)
  setnames(insula, seq_len(ncol(insula)), expnames)
  insula[["bin"]] <- idx[["V2"]]
  
  # Match to original bins
  m <- insula[, match(bin + window - 1, explist[[1]]$IDX$V4)]
  insula <- cbind(insula, explist[[1]]$IDX[m, list(V1, V2, V3, V4)])
  insula <- insula[, -"bin", with = FALSE]
  
  # Format output
  setcolorder(insula, c((nexp + 1):ncol(insula), 1:nexp))
  setnames(insula, 1:4, c("chrom", "start", "end", "bin"))
  
  # Do log2 value over median normalisation
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
  }
  insula <- as.data.frame(insula)

  # Format output
  structure(list(insula_score = insula),
            PACKAGE = "GENOVA",
            colours = cols,
            class = "IS_discovery",
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
