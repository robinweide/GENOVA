#' Compartment versus compartment scores
#'
#' Produces a compartment score quantile versus quantile analysis in which
#' observed over expected contacts are compared between the different quantiles
#' of the compartment scores.
#'
#' @param explist Either a single GENOVA '\code{contacts}' object or a list of
#'   GENOVA '\code{contacts}' objects.
#' @param CS_discovery A signed '\code{CS_discovery}' produced by the
#'   \code{compartment_score} function.
#' @param bins An \code{integer} of length 1 setting the number of quantiles the
#'   comparment score should be divided into.
#' @param dist_thres A \code{numeric} of length two noting the lower and upper 
#'   limit of distances in basepairs to consider. Defaults to 
#'   \code{c(-Inf, Inf)} to include all distances.
#'
#' @details Per chromosome arm, compartment scores are divided in quantile bins.
#'   Subsequently, the average observed over expected score is calculated for
#'   every pairwise combination of bins. The observed over expected score is
#'   calculated as the contacts divided by the average of contacts for the same
#'   distance from the diagonal.
#'
#'   Choosing a '\code{bins}' of \code{5} will produce results similar to
#'   Flyamer \emph{et al}. (2007), while setting '\code{bins}' to \code{100}
#'   will produce results similar to Bonev \emph{et al}. (2017)
#'
#' @section Resolution recommendation: The resolution of the
#'   '\code{CS_discovery}' object. Either print the object to see the resolution
#'   or use \code{attr(CS_discovery, "resolution")} for programmatic access.
#'
#' @return A \code{saddle_discovery} object with 1 element:
#' @return \itemize{\item\strong{\code{saddle}}, a \code{data.table} with the following columns: 
#' \describe{
#' \item{\code{exp}}{A \code{character} with the sample names from the 
#' '\code{explist}' argument.}
#' \item{\code{chr}}{A \code{character} with the chromosome names and arms (p or q).}
#' \item{\code{q1}}{An \code{integer} giving the first comparment score quantile bin.
#' Lower values indicate smaller compartment scores than higher values.}
#' \item{\code{q2}}{An \code{integer} giving the second quantile bin.}
#' \item{\code{mean}}{A \code{numeric} with the average observed over expected values
#' at the indicated quantile bins.}
#' }}
#' @export
#'
#' @examples
#' \dontrun{
#' exps <- list(WT_100kb, KO_100kb)
#'
#' # Computing signed compartment scores
#' cs <- comparment_score(exps, bed = H3K4me1_peaks)
#'
#' # Saddle analysis
#' sadl <- saddle(exps, cs)
#'
#' # Visualising results
#' visualise(sadl)
#' }
saddle <- function(explist, CS_discovery, bins = 10L, 
                   dist_thres = c(-Inf, Inf)) {

  explist <- check_compat_exp(explist)
  expnames_list <- names(explist)
  expnames_exp <- vapply(explist, attr, character(1), "samplename")

  scores <- as.data.table(copy(CS_discovery$compart_scores))[order(chrom, start)]
  discnames <- utils::tail(colnames(scores), length(explist))
  
  # Check argument compatability
  if (all(expnames_list == discnames) && !is.null(expnames_list)) {
    expnames <- expnames_list
  } else if (all(expnames_exp %in% discnames) && 
             all(discnames %in% expnames_exp)) {
    order <- match(expnames_exp, discnames)
    explist <- explist[order]
    expnames <- expnames_exp[order]
  } else {
    stop("The samples in 'explist' should match samples in the 'CS_discovery'.")
  }
  
  if (attr(explist[[1]], "resolution") != attr(CS_discovery, "resolution")) {
    stop(paste("The samples in 'explist' should be at the same resolution as",
               "the 'CS_discovery' object."))
  }
  if (!attr(CS_discovery, "signed")) {
    warning(paste0("The 'CS_discovery' object is unsigned. It is strongly ",
                   "recommended that comparment scores are signed with ",
                   "`sign_comparmentscore()`."), call. = FALSE)
    if (ncol(scores) > 5) {
      keep <- !is.na(scores[[5]]) & scores[["chrom"]] == scores[["chrom"]][1]
      cor <- cor(scores[[5]][keep], scores[[6]][keep])
      if (cor < 0) {
        stop("Insufficient correlation between experiments to continue.",
             call. = FALSE)
      }
    }
  }
  
  # Check if distances need to be filtered
  dist_thres <- sort(dist_thres) / resolution(explist[[1]])
  filter_dist <- dist_thres[1] > 0 | dist_thres[2] < Inf
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Assign chromosome arms
  scores[, part := inverse.rle(attr(CS_discovery, "partitioning"))]
  # Remove centromeres and NAs
  scores <- scores[!endsWith(part, "centro")]
  scores <- scores[!is.na(eval(as.symbol(expnames[1]))), ]
  # Remove chromosomes with less than #bins scores
  use_chrom <- scores[, length(start), by = part]
  use_chrom <- use_chrom[["part"]][use_chrom[["V1"]] > bins]
  quants <- scores[part %in% use_chrom,]
  
  # Setup quantile bins
  qbins <- seq(0, 1, length.out = bins + 1)
  
  # Assign quantiles to scores
  for (i in expnames) {
    i <- as.symbol(i)
    quants[, 
           as.character(i) := as.numeric(pmin(
             findInterval(eval(i), 
                          stats::quantile(eval(i), 
                                          qbins, 
                                          na.rm = TRUE),
                          left.open = TRUE, 
                          rightmost.closed = TRUE), 
             bins)), 
           by = "part"]
    
  }
  
  setkey(quants, bin)
  
  # Get all cis-indices
  all_cis <- scores[, CJ(V1 = bin, V2 = bin), by = list(chr = part)]
  all_cis <- all_cis[V2 > V1]
  setkey(all_cis, V1, V2)

  out <- lapply(explist, function(exp) {
    i <- as.symbol(attr(exp, "samplename"))
    quant <- quants[, list(as.integer(bin), binned = as.integer(eval(i)))]
    setkey(quant, V1)
    
    dat <- exp$MAT[all_cis,]
    dat[["V3"]][is.na(dat[["V3"]])] <- 0
    
    # Calculate distances
    dat[["D"]] <- dat[, abs(V1 - V2)]
    if (filter_dist) {
      dat <- dat[D >= dist_thres[1] & D <= dist_thres[2]]
    }
    
    
    # Calculate observed / expected
    dat[, V3 := (V3 / mean(V3)), by = c("chr", "D")]
    dat[["V3"]][!is.finite(dat[["V3"]])] <- 1
    

    # Assign quantiles
    dat[, V1 := quant[list(dat[["V1"]]), binned]]
    dat[, V2 := quant[list(dat[["V2"]]), binned]]
    
    # Sort by triangle
    dat[["q1"]] <- dat[, pmin(V1, V2)]
    dat[["q2"]] <- dat[, pmax(V1, V2)]
    
    setkey(dat, chr, q1, q2)
    
    # Take mean of observed / expected per quantile
    dat[, mean(V3), by = list(chr, q1, q2)]
  })
  
  # Format output
  names(out) <- expnames
  out <- rbindlist(out, use.names = FALSE, idcol = TRUE)
  setnames(out, c(1,5), c("exp", "mean"))
  
  structure(list(saddle = out),
            package = "GENOVA",
            resolution = attr(explist[[1]], "resolution"),
            class = c("saddle_discovery", "discovery"))
}
