#' Intra and inter TAD contacts
#'
#' Calculate the coverage over TADs and between a TAD and its neighbours.
#'
#' @param explist Either a single GENOVA \code{contacts} object or a list of
#'   GENOVA \code{contacts} objects.
#' @param tad_bed A BED-formatted \code{data.frame} with three columns
#'   containing the chromosome name, start- and end-positions of non-overlapping
#'   TADs.
#' @param max_neighbour An \code{integer} of length 1 indicating how many
#'   neighbouring TADs to take into account.
#'
#' @return An \code{IIT_discovery} object containing the following slots:
#'   \describe{ \item{results}{A \code{data.table} with two columns of left and
#'   right TAD identifiers and scores for the experiments in \code{explist}.}
#'   \item{tads}{The \code{tad_bed} argument supplemented with TAD identifiers
#'   used in the \code{results} slot.}}
#' @export
#'
#' @details The TADs are expected to be consecutive: the 3' border of a TAD is
#'   the same position as the 5' border of the subsequent TAD.
#'
#'   Although the \code{explist} argument is allowed to be a different
#'   resolution than whence the TAD calls came, it is discouraged to use higher
#'   basepair resolution for the \code{explist} argument compared to the TAD
#'   calls for accuracy reasons.
#'
#' @section Resolution recommendation: The same resolution as whence the TADs
#'   were called.
#'
#' @seealso For calling TADs with insulation scores, see
#'   \code{\link[GENOVA]{insulation_score}()} and
#'   \code{\link[GENOVA]{call_TAD_insulation}}. For calling TADs with HiCseg,
#'   see \code{\link[GENOVA]{HiCseg.callTAD}}.
#'
#' @examples
#' \dontrun{
#' # Getting TAD calls
#' insula <- insulation_score(WT_20kb)
#' tads <- call_TAD_insulation(insula)[[1]]
#'
#' # Calculating the inter-intra TAD contacts
#' iit <- intra_inter_TAD(list(WT_20kb, KO_20kb), tads)
#' }
intra_inter_TAD <- function(explist, tad_bed, max_neighbour = 5) {
  
  # Verify experiment compatability
  explist <- check_compat_exp(explist)
  
  # Grab some relevant info
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  cols <- vapply(explist, attr, character(1L), "colour")
  
  # Set data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Translate TADs to bins
  tidx <- data.table(
    start = bed2idx(explist[[1]]$IDX, tad_bed, "start"),
    end = bed2idx(explist[[1]]$IDX, tad_bed, "end"),
    id = seq_len(nrow(tad_bed))
  )
  tidx <- tidx[order(start),]
  
  # Setup bins that lie within TADs
  include <- tidx[,seq.int(start, end - 1), by = id][, list(V1 = unique(V1))]
  setkeyv(include, "V1")
  
  results <- lapply(seq_along(explist), function(i) {
    # Include only bins that are in TADs
    mat <- explist[[i]]$MAT[include, nomatch = NULL]
    mat <- mat[include, on = c(V2 = "V1"), nomatch = NULL]
    # Translate bins to TAD order
    mat[, V1 := findInterval(V1, tidx$start)]
    mat[, V2 := findInterval(V2, tidx$start)]
    # Take only neighbouring TADs
    mat <- mat[abs(V2 - V1) <= max_neighbour]
    # Sum contacts by TAD IDs
    mat[, setNames(list(sum(V3)), expnames[i]), by = c("V1", "V2")]
  })
  
  # Format results
  res <- results[[1]]
  for (i in tail(seq_along(results), -1)) {
    res <- merge(res, results[[i]], by = c("V1", "V2"))
  }
  setnames(res, 1:2, c("x", "y"))
  
  # Translate TAD order to TAD id
  res[, c("x", "y") := list(tidx[x, id], tidx[y, id])]
  
  # Exclude trans-interTADs
  res <- res[tad_bed[x, 1, drop = TRUE] == tad_bed[y, 1, drop = TRUE]]
  res <- res[order(x, y)]
  
  # Construct output
  structure(
    list(
      results = res,
      tads = cbind(tad_bed, id = seq_len(nrow(tad_bed)))
    ),
    class = c("IIT_discovery", "discovery"), 
    package = "GENOVA",
    colours = cols,
    resolution = attr(explist[[1]], "resolution")
  )
}






