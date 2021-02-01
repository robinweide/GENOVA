# Main --------------------------------------------------------------------

#' Compartment score
#'
#' \code{compartment_score} calculates a compartment score per chromosome arm.
#'
#' @param explist Either a single GENOVA \code{contacts} object or list of
#'   GENOVA \code{contacts}.
#' @param ev An \code{integer} of length one, indicating which eigenvector to
#'   retrieve. Useful in case the first eigenvector does not resemble the
#'   compartmentalisation.
#' @inheritParams sign_compartmentscore
#'
#' @details The produre is based on the algorithm of G. Fudenberg, which
#'   performs an eigenvector decomposition of the observed / expected matrix
#'   minus one. Centromeres will be skipped.
#'
#' @note When the '\code{bed}' or '\code{bedgraph}' arguments are both
#'   \code{NULL}, the function returns an unsigned \code{CS_discovery} object
#'   with the attribute '\code{signed = FALSE}'. Signing the compartment scores
#'   is necessary to have positive values correspond to more active A
#'   compartments and negative values to more inactive B compartments.
#'
#' @section Resolution recommendation: 100kb-150kb. Data with <20kb resolution
#'  will return an error due to sparsity.
#'
#' @seealso \code{\link[GENOVA]{sign_compartmentscore}} for getting signed
#'   compartment scores and appropriate '\code{bed}'. and '\code{bedgraph}'
#'   arguments.
#'
#' @return A \code{CS_discovery} object with 1 element.
#' @return \itemize{\item\strong{\code{compart_scores}}, a \code{data.frame} 
#' with the following columns:
#' \describe{
#' \item{\code{chrom}}{A \code{character} with chromosome names.}
#' \item{\code{start}}{An \code{integer} with start positions in the 
#' chromosome.}
#' \item{\code{end}}{An \code{integer} with end positions in the chromosome.}
#' \item{\code{bin}}{An \code{integer} giving the index of the genomic 
#' position.}
#' \item{\emph{samplename_1}}{A \code{numeric}, the calculated compartment score 
#' for the first '\code{explist}' entry. Column name is eponymous with entries 
#' in '\code{explist}'.}
#' \item{\emph{samplename_n} (Optional)}{A \code{numeric}, the calculated 
#' compartment scores for subsequent '\code{explist}' entries.}
#' }
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Proper compartment scores
#' cs <- compartment_score(list(WT_100kb, KO_100kb), bed = H3K4me1_peaks)
#'
#' # Doing the eigenvector decomposition only yields unsigned scores
#' cs <- compartment_score(list(WT_100kb, KO_100kb))
#'
#' # Signing the scores
#' cs <- sign_compartmentscore(cs, bed = H3K4me1_peaks)
#' }
compartment_score <- function(explist, ev = 1, bed = NULL, bedgraph = NULL,
                              ref = 1) {
  
  explist <- check_compat_exp(explist)
  
  if (attr(explist[[1]], "resolution") <= 19999) {
    stop(paste0("Data at this resolution is unfit to compute",
                "meaningful compartment scores. Try data with larger bins."))
  }
  
  znormed <- vapply(explist, attr, logical(1L), "znorm")
  if (any(znormed)) {
    stop(paste0("Compartment scores should not be computed on Z-score ",
                "normalised experiments. Compartment scores are calculated ",
                "on observed over expected matrices and Z-score normalised ",
                "experiments are in essence already observed over expected"),
         call. = FALSE)
  }
  
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Control BLAS threads
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    blas_threads <- RhpcBLASctl::blas_get_num_procs()
    on.exit(RhpcBLASctl::blas_set_num_threads(blas_threads))
    RhpcBLASctl::blas_set_num_threads(1)
  }
  
  # Identify centromeres
  idx <- copy(explist[[1]]$IDX)

  # Parse centromere information for experiments
  centros <- lapply(explist, `[[`, "CENTROMERES")
  widths <- centros[[1]][, end - start]
  names(centros) <- expnames
  centros <- rbindlist(centros, use.names = TRUE)
  centros <- centros[, list(start = min(start), end = max(end)), by = chrom]
  setkey(centros, chrom)
  newwidths <- centros[, end - start]
  if (any(newwidths > 2 * widths)) {
    warning("Centromeres of experiments are probably incompatible.")
  }
  
  partitioning <- partition_chromosomes(idx, centros)
  idx[["V1"]] <- inverse.rle(partitioning)
  idx <- idx[!endsWith(V1, "centro")]
  # Filter out too-small chromosomes
  idx <- idx[, if(.N > 2) .SD, by = V1]
  
  # Select all cis, non-centromere bins
  all_cis <- idx[, CJ(V1 = V4, V2 = V4), by = list(chr = V1)]
  all_cis <- all_cis[V2 > V1]
  setkey(all_cis, V1, V2)

  # Calculate eigenvalues
  eigs <- lapply(explist, function(exp){
    dat <- exp$MAT[all_cis,]
    dat[["V3"]][is.na(dat[["V3"]])] <- 0
    
    # Calculate distances
    dat[["D"]] <- dat[, abs(V1 - V2)]
    
    # Calculate observed / expected - 1
    dat[, V3 := (V3 / mean(V3) - 1), by = c("chr", "D")]
    dat[["V3"]][!is.finite(dat[["V3"]])] <- 0
    
    # Eigenvalue decomposition
    dat <- dat[, uptrimat2eig(V1, V2, V3, ev), by = chr]
    dat
  })
  
  # Combine results
  bins <- eigs[[1]]$bin
  evs <- lapply(eigs, `[[`, "ev")
  evs <- cbind.data.frame(bins, do.call(cbind, evs))
  evs <- as.data.table(evs)
  colnames(evs)[2:(length(expnames) + 1)] <- expnames
  idx <- copy(explist[[1]]$IDX)
  out <- merge(idx, evs, by.x = "V4", by.y = "bins", all = TRUE)
  setcolorder(out, neworder = c(2,3,4,1,5:ncol(out)))
  names(out)[1:4] <- c("chrom", "start", "end", "bin")
  
  cols <- vapply(explist, attr, character(1L), "colour")
  
  out <- structure(
    list(compart_scores = as.data.frame(out)),
    package = "GENOVA",
    colours = cols,
    class = c("CS_discovery", "genomescore_discovery", "discovery"),
    resolution = attr(explist[[1]], "resolution"),
    partitioning = partitioning,
    signed = FALSE
  )
  
  if (!is.null(bedgraph) | !is.null(bed)) {
    out <- sign_compartmentscore(out, bed = bed, bedgraph = bedgraph)
  }
  
  out
}

# Visible Helpers --------------------------------------------------------------

#' Sign compartment score.
#'
#' Takes pre-computed eigenvectors per chromosome arm and adjust the sign based
#' on metrics of active chromatin.
#'
#' @param CS_discovery A \code{CS_discovery} object as produced by the
#'   \code{\link[GENOVA]{compartment_score}} function.
#' @param bed A \code{data.frame} with 3 columns in BED format, containing peaks
#'   of active chromatin marks. Mutually exclusive with the '\code{bedgraph}'
#'   argument.
#' @param bedgraph A \code{data.frame} with 4 columns in bedGraph format,
#'   containing scores correlated with active chromatin. Mutually exclusive with
#'   the '\code{bed}' argument.
#' @param verbose A \code{logical} of length 1. If \code{TRUE}, reports when
#'   '\code{bedgraph}' correlations are weak or uncomputable.
#' @param ref An \code{integer} of length 1 giving a sample index that can be
#'   used as a reference. Set to \code{NULL} to sign every sample individually.
#'
#' @details Signing compartment scores is necessary to have the compartment
#'   scores correspond to A and B compartments.
#'
#'   Some decent '\code{bed}'/'\code{bedgraph}' arguments to use are ChIP-seq
#'   peaks of H3K4me1, gene density, (inverted) Lamina DamID signal, GC content
#'   or early replication timing metrics.
#'
#' @return A \code{CS_discovery} object with \code{signed} attribute set to
#'   \code{TRUE}
#'
#' @seealso \code{\link[GENOVA]{compartment_score}} for computing the
#'   eigenvectors.
#' @export
#'
#' @examples
#' \dontrun{
#' # Doing the eigenvector decomposition only yields unsigned scores
#' cs <- compartment_score(list(WT_100kb, KO_100kb))
#'
#' # Signing the scores
#' cs <- sign_compartmentscore(cs, bed = H3K4me1_peaks)
#' }
sign_compartmentscore <- function(CS_discovery,
                                  bed = NULL,
                                  bedgraph = NULL,
                                  verbose = FALSE,
                                  ref = 1) {
  # Checks
  if (!inherits(CS_discovery, "CS_discovery")) {
    stop(paste0("'sign_compartmentscore()' only works with 'CS_discovery'",
                "objects, which can be generated with the", 
                "'compartment_score()' function."), call. = FALSE)
  }
  if (!is.null(bed) & !is.null(bedgraph)) {
    stop(paste0("Please use either the 'bed' or 'bedgraph' argument."))
  }
  if (is.null(bed) & is.null(bedgraph)) {
    stop("Please supply the 'bed' or 'bedgraph' argument.")
  }
  
  # Setup data
  CS <- as.data.table(CS_discovery$compart_scores)
  CS <- CS[order(chrom, start)]
  expnames <- tail(colnames(CS), -4)
  CS[["part"]] <- inverse.rle(attr(CS_discovery, "partitioning"))
  
  # Setup bed/bedgraph matching
  idx <- CS[, list(chrom, start, end, bin)]
  names(idx) <- paste0("V", 1:4)
  
  if (is.null(ref)) {
    ref <- seq_along(expnames)
  }
  
  if (!is.null(bed)) {
    # Count bed entries per bin
    bed_idx <- table(GENOVA::bed2idx(idx, bed))
    
    # Match bins from counts to compartment scores
    CS[["bedcount"]] <- as.numeric(bed_idx[match(CS[["bin"]], names(bed_idx))])
    CS[["bedcount"]][is.na(CS[["bedcount"]])] <- 0
    part <- as.data.frame(CS)
    part <- part[!is.na(part[, expnames[[1]]]), ]
    
    # Loop over arms
    split <- split(part, part$part, drop = TRUE)
    split <- lapply(split, function(chrpart) {
      count <- chrpart[["bedcount"]]
      # Per experiment, decide to flip
      chrpart[, expnames[ref]] <- lapply(chrpart[, expnames[ref], drop = FALSE], 
                                    function(score) {
        up <- score > 0 & !is.na(score)
        down <- score < 0 & !is.na(score)
        if (sum(up) > 0 & sum(down) > 0) {
          # Is there a difference?
          test1 <- wilcox.test(count[up], 
                               count[down])
          # Is the orientation wrong?
          test2 <- median(count[up]) < median(count[down])
          flip <- test1$p.value < 1e-2 & test2
          if (flip) {
            # Flip score
            return(- score)
          }
        }
        return(score)
      })
      chrpart
    })
  } else if (!is.null(bedgraph)) {
    # Match values to bins
    bedgraph$bin <- GENOVA::bed2idx(idx, bedgraph)
    bedgraph <- as.data.table(bedgraph)
    names(bedgraph)[1:4] <- paste0("V", 1:4)
    
    # Take mean per bin
    bedgraph <- bedgraph[, mean(V4), by = bin]
    
    # Match means to compartment scores
    CS[["graph"]] <- bedgraph[match(CS[["bin"]], bedgraph[["bin"]]), V1]
    part <- as.data.frame(CS)
    part <- part[!is.na(part[, expnames[[1]]]), ]
    
    # Loop over arms
    split <- split(part, part$part, drop = TRUE)
    split <- lapply(split, function(chrpart) {
      chrpart[, expnames[ref]] <- lapply(
        chrpart[, expnames[ref], drop = FALSE], 
        function(score) {
          # Compute correlation
          cor <- suppressWarnings(cor(score, chrpart$graph, 
                                      method = "spearman"))
          if (is.na(cor) & verbose) {
            message(paste0("No correlation could be computed for '", 
                           chrpart$chrom[1], "'."))
            return(score)
          }
          if (abs(cor) < 0.25 & verbose) {
            message(paste0("The correlation with 'bedgraph' on '", 
                           chrpart$chrom[1], "' was found to be less than ",
                           "|0.25|."))
          }
          # Flip if correlation is negative
          return(score * sign(cor))
        })
      return(chrpart)
    })
  } else {
    stop("Some unknown error has occurred", call. = FALSE)
  }
  
  if (length(ref) == 1 && length(expnames) > 1) {
    split <- lapply(split, function(chrpart) {
      refval <- chrpart[[expnames[ref]]]
      chrpart[expnames[-ref]] <- lapply(chrpart[expnames[-ref]], function(val) {
        cor <- suppressWarnings(cor(val, refval))
        return(val * sign(cor))
      })
      chrpart
    })
  }
  
  # Merge data together
  out <- unsplit(split, part$part, drop = TRUE)
  out <- as.data.table(out[, 4:(length(expnames) + 4)])
  names(idx) <- c("chrom", "start", "end", "bin")
  out <- merge(idx, out, by = "bin", all = TRUE)
  setcolorder(out, neworder = c(2,3,4,1,5:ncol(out)))
  
  structure(list(compart_scores = as.data.frame(out)),
            class = c("CS_discovery", "genomescore_discovery", "discovery"),
            colours = attr(CS_discovery, "colours"),
            package = "GENOVA",
            resolution = attr(CS_discovery, "resolution"),
            partitioning = attr(CS_discovery, "partitioning"),
            signed = TRUE)
}

# Utilities ---------------------------------------------------------------


#' Sparse upper triangular symmetric matrices to eigenvector
#'
#' Used to calculate compartment scores within a data.table.
#'
#' @param x row indices of the matrix
#' @param y column indices of the matrix
#' @param value elements of the matrix
#' @param ev A numeric of length 1 indicating which eigenvector to retrieve
#'
#' @details x, y and value should be parallel to oneanother. The offset of x and
#'   y relative to 1 are removed prior to constructing a dense square matrix and
#'   added back in just before returning the output.
#'
#' @return A data.table containing a bin column and eigenvector column.
#' @noRd
#' @noMd
#' @keywords internal
uptrimat2eig <- function(x, y, value, ev = 1) {
  # Symmetric eigen works on lower triangle, so switch x and y
  i <- cbind(y, x)
  # Re-index to start at 1
  i <- i - {mini <- min(i)} + 1
  m <- matrix(0, max(i), max(i))
  m[i] <- pmin(value, quantile(value, 0.995))
  eig <- eigen(m, symmetric = TRUE)
  eig <- eig$vectors[, ev] * sqrt(eig$values[ev])
  data.table(bin = seq_along(eig) + mini - 1,
             ev = eig)
}
