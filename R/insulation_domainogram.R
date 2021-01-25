#' Domainogram of insulation
#'
#' Creates a domainogram of insulation scores for a genomic region of interest
#' by calculating the insulation scores for a range of sliding square sizes.
#'
#' @param explist Either a single GENOVA \code{contacts} object or list of
#'   GENOVA \code{contacts} objects.
#' @param chrom One of the following: \itemize{
#'     \item{A \code{character} of length one indicating a chromosome name.}
#'     \item{A 3-column, 1-row \code{data.frame} in BED-format.}
#'     \item{A single \code{character} describing a locus in UCSC notation, e.g. 
#'     \code{"chr1:30,000,000-40,000,000"}.}
#'   } The latter two options automatically provide the \code{start} and 
#'   \code{end} arguments too.
#' @param start,end An \code{integer} of length 1 with a region's start and end 
#'   position in basepairs.
#' @param window_range An \code{integer} vector of length 2 noting the minimum
#'   and maximum size of the sliding square.
#' @param step An \code{integer} of length 1 with the step size for incrementing
#'   the size of the sliding square.
#'
#' @return A \code{domainogram_discovery} object containing a \code{data.frame}
#'   with the insulation scores for the region of interest at various sizes of
#'   the sliding square.
#' @export
#'
#' @examples
#' \dontrun{
#' # Making a domainogram
#' domgram <- insulation_domainogram(list(WT = WT_20kb, KO = KO_20kb),
#'                                   chrom = "chr7", start = 25e6, end = 30e6)
#'                                  
#' # Plotting the domainogram
#' visualise(domgram)
#' }
insulation_domainogram <- function(
  explist, chrom, start, end,
  window_range = c(5, 101), step = 2
) {
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  loc <- standardise_location(chrom, start, end, singular = TRUE)
  
  # Setup experiments
  explist <- check_compat_exp(explist)
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  
  # Setup windows
  windows <- sort(window_range)
  maxwin <- windows[2]
  windows <- seq(windows[1], windows[2], by = step)
  
  # Find ROI
  bed <- cbind.data.frame(loc$chrom, loc$start, loc$end)
  bed <- sort(c(bed2idx(explist[[1]]$IDX, bed, "start"),
                bed2idx(explist[[1]]$IDX, bed, "end")))
  pos_start <- explist[[1]]$IDX$V3[match(bed[1], explist[[1]]$IDX$V4)]
  res <- attr(explist[[1]], "resolution")
  idx <- seq(bed[1] - maxwin, bed[2] + maxwin)
  idx <- CJ(V1 = idx, V2 = idx)
  
  # Setup grid
  grid <- CJ(V1 = seq_len(maxwin), V2 = seq_len(maxwin))
  grid[, window := pmax(maxwin - V1 + 1, V2)]
  grid[, c("V1", "V2") := list(V1 - 1, V2 + maxwin - 1)]
  grix <- grid[["V1"]]
  griy <- grid[["V2"]]
  griz <- grid[["window"]]
  
  # Setup grid x ROI
  looper <- data.table(id = seq(bed[1] - maxwin, bed[2] - maxwin))
  looper <- looper[, list(V1 = grix + .BY[[1]],
                          V2 = griy + .BY[[1]],
                          window = griz), by = list(id)]
  setcolorder(looper, c(2,3,1))
  setkey(looper, V1, V2)
  
  # Setup data
  dat <- lapply(seq_along(explist), function(i){
    mat <- explist[[i]]$MAT[idx, nomatch = NULL]
    mat[, exp := i]
  })
  dat <- rbindlist(dat)
  setkey(dat, V1, V2)
  dat <- dat[looper, list(id, exp, window, V3), 
             allow.cartesian = TRUE, nomatch = NULL]
  
  # Compute insulation scores
  scores <- lapply(windows, function(w) {
    dat[window <= w, list(ins = sum(V3, na.rm = TRUE), window = w),
        by = c("id", "exp")]
  })
  scores <- rbindlist(scores)
  scores[, ins := ins / window ^ 2]
  scores[, ins := log2(ins / mean(ins)), by = c("exp", "window")]
  
  # Format data
  scores$exp <- expnames[scores$exp]
  scores[, id := (id - bed[1] + maxwin - 1) * res + pos_start]
  scores <- dcast(scores, window + id ~ exp, value.var = "ins")
  scores <- as.data.frame(scores)
  colnames(scores) <- c("window", "position", expnames)
  structure(list(scores = scores), 
            class = c("domainogram_discovery", "genomescore_discovery", 
                      "discovery"),
            package = "GENOVA",
            chrom = loc$chrom, resolution = res)
}

