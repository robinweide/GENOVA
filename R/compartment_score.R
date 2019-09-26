compartment_score <- function(explist, ev = 1) {
  
  explist <- GENOVA:::check_compat_exp(explist)
  
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
  cen <- unique(explist[[1]]$MAT, by = 1)[["V1"]]
  idx[, V5 := as.numeric(V4 %in% cen)]
  idx <- split(idx, idx[["V1"]])
  idx <- lapply(idx, function(x) {
    y <- x[["V5"]]
    # Decide what is not centromere
    cs <- cumsum(c(y[1], diff(y)))
    rle <- rle(cs)
    if (!any(rle$values == 0)) {
      return(x)
    }
    attr <- !with(rle, values == 0 & lengths == max(lengths[values == 0]))
    keep <- rep.int(attr, rle$lengths)
    
    # Diversify chromosome names if not centromere
    rle <- rle(keep)
    rle$values <- with(rle, ifelse(values, LETTERS[cumsum(values)], "centro"))
    x[, V1 := paste0(V1, rep.int(rle$values, rle$lengths))]
    x
    # x[, keep := keep]
    # Filter out centromere
    # x <- x[keep,]
  })
  idx <- rbindlist(idx)
  partitioning <- rle(idx[["V1"]])
  idx <- idx[!endsWith(V1, "centro")]
  
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
  structure(list(compart_scores = out), 
            class = c("CS_discovery"),
            res = attr(explist[[1]], "res"),
            partitioning = partitioning,
            signed = FALSE)
}

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

sign_compartmentscore <- function(CS_discovery,
                                  bed = NULL,
                                  bedgraph = NULL) {
  if (!inherits(CS_discovery, "CS_discovery")) {
    stop(paste0("'sign_compartmentscore()' only works with 'CS_discovery'",
                "objects, which can be generated with the", 
                "'compartment_score()' function."), call. = FALSE)
  }
  if (!is.null(bed) & !is.null(bedgraph)) {
    stop(paste0("Please use either the 'bed' or 'bedgraph' argument."))
  }
  CS <- CS_discovery$compart_scores
  CS <- CS[order(chrom, start)]
  expnames <- tail(colnames(CS), -4)
  CS[["part"]] <- inverse.rle(attr(CS_discovery, "partitioning"))
  
  idx <- CS[, list(chrom, start, end, bin)]
  names(idx) <- paste0("V", 1:4)
  
  bed_idx <- table(GENOVA::bed2idx(idx, bed))
  
  CS[["bedcount"]] <- as.numeric(bed_idx[match(CS[["bin"]], names(bed_idx))])
  CS[["bedcount"]][is.na(CS[["bedcount"]])] <- 0
  part <- as.data.frame(CS)
  part <- part[!is.na(part[, expnames[[1]]]), ]
  
  
  split <- split(part, part$part, drop = TRUE)
  split <- lapply(split, function(chrpart) {
    count <- chrpart[["bedcount"]]
    chrpart[, expnames] <- lapply(chrpart[, expnames], function(score) {
      up <- score > 0 & !is.na(score)
      down <- score < 0 & !is.na(score)
      if (sum(up) > 0 & sum(down) > 0) {
        test1 <- wilcox.test(count[up], 
                             count[down])
        test2 <- median(count[up]) < median(count[down])
        flip <- test1$p.value < 1e-5 & test2
        if (flip) {
          return(- score)
        }
      }
      return(score)
    })
    chrpart
  })
  out <- unsplit(split, part$part, drop = TRUE)
  out <- as.data.table(out[, 4:(length(expnames) + 4)])
  names(idx) <- c("chrom", "start", "end", "bin")
  out <- merge(idx, out, by = "bin", all = TRUE)
  
  structure(list(compart_score = out),
            class = "CS_discovery",
            res = attr(CS_discovery, "res"),
            partitioning = attr(CS_discovery, "partitioning"),
            signed = TRUE)
}

bed <- "/DATA/projects/Hap1/ChIP-seq/completeSet_WTWaplMed12DKO/MACS2/WT_K4me1_peaks.narrowPeak"
bed <- read.table(bed)
