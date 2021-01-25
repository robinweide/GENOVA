#' Directionality index
#'
#' The directionality index quantifies the degree of bias between upstream and
#' downstream interactions given a bin on the diagonal. Such biases become
#' apparent near the periphery of TADs: the upstream portion of a TAD interacts
#' more with the downstream bins and inversely, the downstream portion of a TAD
#' interacts more with the upstream bins.
#'
#' @param explist Either a single GENOVA \code{contacts} object or list of
#'   GENOVA \code{contacts}.
#' @param range An \code{integer} of length 1 indicating how many bins from the
#'   diagonal should be considered.
#'
#' @details The directionality index is computed as described in Dixon \emph{et
#'   al}. (2012):
#'
#'   \deqn{DI = (\frac{B - A}{|B - A|})}{DI = ((B - A) / (|B - A|)) ((A - E)^2 /
#'   E + (B - E)^2 / E)}
#'
#'   Wherein \eqn{A} and \eqn{B} are the sum of the contacts up- and downstream
#'   of a bin on the diagonal respectively, and \eqn{E = (A + B)/2}.
#'
#'   The first part signs the second part of the equation by the direction of
#'   the effect. The \eqn{\chi^2} test statistics can be recognised in the
#'   second part of the equation with a null-hypothesis that the upstream and
#'   downstream signal is equal.
#'
#'   The authors originally used a 40kb matrix with a 2Mb range, equivalent to
#'   \code{range = 50}.
#'
#' @return A \code{DI_discovery} object containing a directionality index for
#'   every informative bin.
#' @export
#'
#' @examples
#' \dontrun{
#' # As original authors
#' di <- direct_index(list(WT_40kb, KO_40kb), range = 50)
#' 
#' # Plotting the DI
#' visualise(di)
#' }
direct_index <- function(explist, range = 100) {
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Pre-flight checklist
  explist <- check_compat_exp(explist)
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  idx <- copy(explist[[1]]$IDX)
  
  # Get data +- range from diagonal
  mat <- lapply(seq_along(explist), function(i){
    thismat <- explist[[i]]$MAT
    thismat <- thismat[abs(V1 - V2) < range]
    thismat[, exp := i]
  })
  mat <- rbindlist(mat)
  
  # Discard diagonal and trans data
  mat <- mat[V1 != V2]
  keep <- idx[match(mat$V1, V4), V1] == idx[match(mat$V2, V4), V1]
  mat <- mat[keep]
  
  # Calculate upstream and downstream sums
  AA <- mat[, list(A = sum(V3)), keyby = c("V1", "exp")]
  BB <- mat[, list(B = sum(V3)), keyby = c("V2", "exp")]
  
  # Combine upstream and downstream
  setnames(BB, 1, "V1")
  calc <- merge(AA, BB, all = TRUE)
  
  # Replace NAs by 0s
  calc[, A := pmax(0, A, na.rm = TRUE)]
  calc[, B := pmax(0, B, na.rm = TRUE)]
  
  # Calculate directionality index
  calc[, E := (A + B) / 2]
  calc[, term1 := (B - A) / abs(B - A)]
  calc[, term2 := ((A - E) ^ 2) / E + ((B - E) ^ 2) / E]
  calc[, DI := term1 * term2]
  
  # Format stuff
  calc <- calc[, list(V1, exp, DI = -DI)]
  calc$exp <- expnames[calc$exp]
  
  # Combine with genomic coordinates
  calc <- dcast(calc, V1 ~ exp, value.var = "DI")
  setkeyv(calc, "V1")
  setkey(idx, V4)
  calc <- idx[calc]
  
  # Format output
  setnames(calc, 1:4, c("chrom", "start", "end", "bin"))
  
  cols <- vapply(explist, attr, character(1), "colour")
  
  structure(list(DI = as.data.frame(calc)),
            colours = cols,
            class = c("DI_discovery", "genomescore_discovery", "discovery"),
            package = "GENOVA",
            resolution = attr(explist[[1]], "resolution"))
}