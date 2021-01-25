#' Call TADs from insulation scores
#'
#' Runs peak detection on the genome wide insulation scores to identify
#' insulation valleys that correspond well to TAD boundaries.
#'
#' @param IS_discovery A GENOVA \code{IS_discovery} object as created by
#'   \code{insulation_score()}.
#' @param method A \code{character} vector of length 1. Currently, only
#'   \code{"crane"} is implemented.
#' @param min_strength A \code{numeric} vector of length 1 setting a threshold
#'   on the boundary strength wherein higher means stronger boundaries.
#'
#' @details The Crane \emph{et al}. (2015) method calculates a delta vector across the
#'   insulation score defined by the difference between 100kb to the right and
#'   left of a central bin. This delta crosses zero at local extrema, of which
#'   the minima are kept. Subsequently, potential boundaries are filtered on the
#'   boundary strength defined as the difference in the delta vector between the
#'   nearest 5' local maximum and 3' local minimum relative to the boundaray.
#'
#' @return A \code{list} containing BED-formatted \code{data.frame}s for each
#'   experiment in the \code{IS_discovery} object.
#' @export
#'
#' @seealso \code{\link[GENOVA]{insulation_score}} for calculating insulation
#'   scores.
#'
#' @examples
#' \dontrun{
#' # Calculating insulation scores
#' ins <- insulation_score(list(WT_20kb, KO_20kb), window = 25)
#' 
#' # Calling TADs from the insulation score
#' tadlist <- call_TAD_insulation(ins)
#' 
#' # Plotting TADs
#' hic.matrixplot(
#'   exp1 = WT_20kb,
#'   chrom = "chr7",
#'   start = 25e6,
#'   end = 30e6,
#'   tads = tadlist$WT_20kb
#' )
#' }
call_TAD_insulation <- function(IS_discovery, method = "crane", 
                                min_strength = 0.1) {
  ins <- as.data.table(IS_discovery$insula_score)
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  expnames <- tail(colnames(ins), -4)
  
  scooch = pmax(100e3 / attr(IS_discovery, "resolution"), 2)
  
  m <- melt.data.table(ins, id.vars = c("chrom", "start", "end", "bin"),
                       measure.vars = expnames)
  m <- m[order(variable, chrom, bin)]
  setnames(m, 5, "exp")
  grouping <- c("chrom", "exp")
  
  # Replace infinite values and NAs
  m[, value := pmin(value, 2 * max(c(0, value[is.finite(value)]))), 
    by = grouping]
  m[, value := pmax(value, min(c(0, value[is.finite(value)]))),
    by = grouping]
  m[is.na(value), value := NaN]
  
  # Standardise insulation scores
  m[, value := scale(value), by = grouping]
  
  # Filter out all-NaN chromosomes
  keep <- m[, list(keep = !all(is.nan(value)), len = .N), 
            by = grouping]
  m <- m[rep(keep$keep, keep$len),]
  
  # Calculate delta
  m[, left := frollmean(value, n = scooch, fill = 0, align = "left"),
    by = grouping]
  m[, right := frollmean(value, n = scooch, fill = 0, align = "right"),
    by = grouping]
  m[, delta := left - right]
  m[, c("left", "right") := NULL]
  
  # Take ends of rles where values are increasing, i.e. valleys
  zerocross <- m[, list(diff = delta > 0), by = grouping]
  zerocross <- zerocross[, rle(diff), by = grouping]
  zerocross <- zerocross[, lengths := cumsum(lengths)]
  zerocross <- zerocross[values == FALSE]
  
  # Search for the truly smallest insulation score within a small window
  shft <- zerocross[, lengths + c(-3:3), by = seq_len(nrow(zerocross))]
  shft[, ins := m[pmax(1, V1), value]]
  shft <- shft[, list(shift = c(-3:3)[which.min(ins)]), by = seq_len]
  zerocross[, lengths := lengths + shft$shift]
  
  # Calculate local extrema in delta space
  loc_extr <- m[, list(diff = diff(c(.Machine$integer.max, delta)) > 0),
              by = grouping]
  loc_extr <- loc_extr[, rle(diff), by = grouping]
  loc_extr <- loc_extr[, lengths := cumsum(lengths)]
  loc_extr <- loc_extr[, list(.SD, type = c("min", "max")[values + 1])]
  setnames(loc_extr, 1:4, c(colnames(zerocross)[1:4]))
  loc_extr <- loc_extr[, score := m[lengths, delta]]
  
  # Split extrema in minima and maxima
  loc_min <- loc_extr[type == "min"]
  loc_max <- loc_extr[type == "max"]
  
  # Finding preceding maxima and succeeding minima
  # Holy data.table roll-join wizardry
  zerocross[, left := loc_max[zerocross, roll = +Inf, 
                              on = c("lengths" = "lengths")]$score]
  zerocross[, right := loc_min[zerocross, roll = -Inf,
                               on = c("lengths" = "lengths")]$score]
  # Technically we're ignoring the subtleties of chromosome starts and ends 
  # in terms of boundary strength by rolling over the entire data, but if 
  # telomeres are not valid boundaries, what are we even doing?
  
  # Filter on boundary strength
  zerocross[, strength := left - right]
  zerocross <- zerocross[strength >= min_strength]
  
  # Translate to TADs
  tads <- m[zerocross$lengths, ]
  tads <- tads[, list(start = c(0, head(end, -1)), end), by = grouping]
  expv <- tads[["exp"]]
  tads[, exp := NULL]
  tads <- split(as.data.frame(tads), expv)
  
  return(tads)
}
