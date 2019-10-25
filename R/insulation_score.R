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
#'   for the genome as a whole.
#'
#' @details Typically, when the sliding square passes a TAD it yields high
#'   insulation scores whereas in between TADs it yields low insulation scores.
#'   Hence, TAD boundaries can be identified as local minima in the insulation
#'   score.
#'
#'   The normalisation is performed by taking the
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} value over median,
#'   wherein the median is either calculated per chromosome or across the
#'   genome.
#'
#' @return A \code{IS_discovery} object
#' @export
#'
#' @examples
#' \dontrun{
#' iscore <- insulation_score(explist(WT_20kb, KO_20kb), window = 20) 
#' }
insulation_score <- function(explist, window = 30, norm_to = c("chromosome", "genome")) {
  
  norm_to <- match.arg(norm_to)
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
      insula[, as.character(xp) := log2(eval(xp) / median(eval(xp), na.rm = TRUE)), 
             by = chrom]
    }
  } else {
    for (i in expnames) {
      xp <- as.symbol(i)
      insula[, as.character(xp) := log2(eval(xp) / median(eval(xp), na.rm = TRUE))]
    }
  }

  # Format output
  structure(list(insula_score = insula),
            PACKAGE = "GENOVA",
            colours = cols,
            class = "IS_discovery",
            resolution = attr(explist[[1]], "resolution"),
            window = window)
}