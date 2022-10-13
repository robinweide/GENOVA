#' Empirical centromeres
#'
#' For each chromosome, finds and reports the largest stretch of unused bins as
#' centromeres.
#'
#' @param MAT The \code{MAT} slot of a GENOVA \code{contacts} object.
#' @param IDX The \code{IDX} slot of a GENOVA \code{contacts} object.
#'
#' @details A bin is seen as unused when the index has an entry in the
#'   \code{IDX} argument but does not have entries in the \code{MAT} argument.
#'
#' @return A \code{data.table} with chromosome, start-index and end-index
#'   columns.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' centros <- find_centromeres(exp$MAT, exp$IDX)
#' }
find_centromeres <- function(MAT, IDX) {
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Prevent accidental modify-in-place
  idx <- copy(IDX)
  
  # Find existing bins
  existing_bins <- unique(MAT, by = 1)[["V1"]]
  idx[, V5 := as.numeric(V4 %in% existing_bins)]
  
  # Loop over chromosomes
  idx <- split(idx, idx[["V1"]])
  centros <- lapply(idx, function(chr) {
    
    exists <- chr[["V5"]]
    
    # Find stretches of empty bins
    rle <- rle(cumsum(c(exists[1], diff(exists))))
    if (!any(rle$values == 0)) {
      return(NULL)
    }
    
    # Decide what is the centromere
    is_centro <- with(rle, values == 0 & lengths == max(lengths[values == 0]))
    is_centro <- rle(rep.int(is_centro, rle$lengths))
    
    # Find start and end bins
    end <- cumsum(is_centro$lengths)
    start <- end - is_centro$lengths + 1
    
    bins <- chr[c(start[is_centro$values], end[is_centro$values]), V4]
    data.table(chrom = chr[1, V1], start = bins[1], end = bins[2])
  })
  
  centros <- rbindlist(centros)
  if (is.null(centros)) {
    return(centros)
  }
  if (!(nrow(centros) > 0)) {
    return(NULL)
  }
  setkey(centros, chrom)
  centros
}

#' Partition chromosomes by arm
#'
#' Names the part of a chromosome before the centromere the p-arm and after the
#' centromere the q-arm.
#'
#' @param IDX The \code{IDX} slot of a GENOVA \code{contacts} object.
#' @param centros A \code{data.table} keyed on the first column, containing
#'   chromosome name, start-index and end-index as columns. Typically from the
#'   \code{find_centromeres} or \code{clean_centromeres} functions.
#'
#' @return A run length encoded partitioning of chromosome arms parallel to the
#'   \code{IDX} argument.
#'   
#' @noRd
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' centros <- find_centromeres(exp$MAT, exp$IDX)
#' parts <- partition_chromosomes(exp$IDX, centros)
#' }
partition_chromosomes <- function(IDX, centros) {
  out <- IDX[, paste0(V1, 
                      ifelse(V4 < centros[.BY, start], "p",
                             ifelse(V4 > centros[.BY, end], "q", "centro"))), 
             by = V1]
  rle(out[[2]])
}

#' Clean up centromeres
#'
#' Centromere positions from the cytobands.txt typically have two entries per
#' centromere. This function is just a failsafe such that these are merged
#' before centromere information is added to a \code{contacts} object.
#'
#' @param centros A BED-like \code{data.frame} with centromere positions.
#' @param IDX The \code{IDX} slot of a GENOVA \code{contacts} object.
#'
#' @return A \code{data.table} with chromosome, start-index and end-index
#'   columns.
#'
#' @noMd
#' @noRd
#' @keywords internal
clean_centromeres <- function(centros, IDX) {
  # Essentially does the same as `reduce(a_granges_object, min.gapwidth = resolution)
  
  centros <- as.data.frame(centros)
  
  resolution <- median(diff(IDX$V2))
  
  # Drop unused factor levels
  if (is.factor(centros[, 1])) {
    centros[, 1] <- droplevels(centros[, 1])
  }
  
  # Set column names and order on chromosome and start
  centros <- setNames(centros, paste0("V", 1:3))
  centros <- centros[order(centros[, 1], centros[, 2]), ]
  
  # Test wether adjacent entries should be merged
  centros$merge <- c(FALSE, vapply(seq_len(nrow(centros) - 1), function(i) {
    
    # Is the difference between end and next start sub-resolution?
    test <- abs(centros[i, 3] - centros[i + 1, 2]) < resolution
    
    # Are the entries on the same chromosome?
    test && centros[i, 1] == centros[i + 1, 1]
  }, logical(1)))
  
  # Merge the TRUE entries
  while (any(centros$merge)) {
    # What is the next entry to be merged?
    i <- which(centros$merge)[1]
    
    # Replace end of previous entry with end of current entry
    centros[i - 1, 3] <- centros[i, 3]
    
    # Remove current entry
    centros <- centros[-i, ]
  }
  
  centros <- data.table(
    chrom = centros[, 1],
    start = bed2idx(IDX, centros, "start"),
    end   = bed2idx(IDX, centros, "end"),
    key = "chrom"
  )
  
  # Return centromeres
  centros
}
