select.sub <- function( data, start, end ){
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  x <- rep(start:end, end-start+1)
  y <- rep(start:end, each=end-start+1)
  data.sub <- data[list(x,y)]
  data.sub <- data.sub[!is.na(data.sub$V3)]
  return(data.sub)
}

select.subset.hicseg <- function(EXP, chr, start, end) {
  sel <- EXP$IDX[V1 == chr & V2 >= start & V2 <= end][["V4"]]
  # sel <- EXP$IDX[EXP$IDX[[1]] == chr & EXP$IDX[[2]] >= start & EXP$IDX[[2]] <= end, 4]
  start.i <- min(sel)
  end.i <- max(sel)
  # create matrix
  r1.s <- select.sub(EXP$MAT, start.i, end.i)
  r2.s <- r1.s
  sub.mat <- matrix(0, ncol = length(start.i:end.i), nrow = length(start.i:end.i))
  # fill the matrix with the first replicate
  index <- cbind(r1.s$V1 - start.i + 1, r1.s$V2 - start.i + 1)
  sub.mat[index] <- r1.s$V3
  # fill the matrix with the second replicate
  index <- cbind(r2.s$V2 - start.i + 1, r2.s$V1 - start.i + 1)
  sub.mat[index] <- r2.s$V3
  pos <- EXP$IDX[EXP$IDX[[4]] %in% (start.i:end.i)][[3]]
  # sub.mat
  list(x = pos, y = pos, z = sub.mat)
}

#' HiCseg.callTAD
#'
#' A wrapper with call TADs with HiCseg of Lévy-Leduc et al. (2018)
#'
#' @param experiment The Hi-C experiment object: produced by construct.experiment().
#' @param chromsToUse The chromosome(s) of interest. Omitting this leads to the computation across all chromosomes.
#' @param binPerBorder Minimum number of bins per border: used to set the maximal
#' number of change-points. The default (3) set this maximal number of borders to nbins/3.
#' @param chunk Use chunks of [chunk]bp instead of the whole arm: the whole arm can be a big burden,
#' so it could be beneficial to use chunks of a set size (e.g. 5e6). But be carefull:
#' the segmentation will be different from the full arm!!!
#' @param BEDcolor Color of items in the resulting BEDPE-file
#' @param verbose Bioinformatics can be a bit lonely: set this to true to get a more chatty function.
#' @return A BEDPE-df
#' @export
HiCseg.callTAD <- function(experiment, chromsToUse = NULL, binPerBorder = 3, chunk = NULL, BEDcolor = "255,255,0", verbose = F) {
 
  try_require('HiCseg', 'HiCseg.callTAD')
  
  # check if chunks are given and provide a warning
  if (!is.null(chunk)) {
    message("You chose a chunked approach.
While it is faster, it does give a different result!")
  }

  # if no chroms are given, take all chroms without chrY amd chrMT
  if (is.null(chromsToUse)) {
    chromsToUse <- experiment$CHRS[!grepl(experiment$CHRS, pattern = "[YM]")]
  }
  idx <- experiment$IDX
  cntro <- experiment$CENTROMERES
  experiment$CENTROMERES <- 
    cntro[, .(chrom = chrom, start = idx[.(start), V2, on = c(V4 = "V1")],
              end = idx[.(end), V3, on = c(V4 = "V1")])]

  DFlist <- list()
  for (chr in chromsToUse) {
    if (verbose) {
      message(chr, ": started")
    }
    # get the bed entry of the centromere
    centChrom <- experiment$CENTROMERES[experiment$CENTROMERES[[1]] == chr, ]

    # throw a hissyfit if chomosome is not found
    if (length(centChrom) == 0) {
      message("There is no centromere-information for ", chr, ".
Skipping this chromosome.")
      next()
    }

    # get chromosome-size
    chromSize <- max(experiment$IDX[experiment$IDX$V1 == chr, 3])

    ###
    # first_arm
    ###
    if (verbose) {
      message(chr, ": P-arm")
    }
    # get start and end of arm
    S <- 0
    E <- centChrom[centChrom[[1]] == chr][[2]]

    m <- NULL # m holds alls boundaries (in bps)
    if (!is.null(chunk)) { # check if chunks are given

      # set chunks
      steps <- seq(S, E, by = chunk / 2)

      tmpM <- list()
      for (STEP in steps) {
        hic.mat <- select.subset.hicseg(experiment, chr = chr, start = STEP, end = STEP + chunk)
        res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
        tmpM[[as.character(STEP)]] <- hic.mat$x[res$t_hat]
      }

      # keep boundaries that are found in two chunks
      m <- as.numeric(names(table(unlist(tmpM))[table(unlist(tmpM)) > 1]))
      # also keep boundaries that are found in the non-overlapping first part of the first window
      m <- c(m, tmpM[[1]][tmpM[[1]] < steps[2]])
    } else {
      hic.mat <- select.subset.hicseg(experiment, chr = chr, start = S, end = E)
      res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
      m <- hic.mat$x[res$t_hat]
    }

    DF <- NULL
    if (length(m) < 2) {
      warning(chr, "'s P-arm has no borders")
    } else {
      mm <- unique(m)
      DF <- data.frame(chr, utils::head(mm, -1), utils::tail(mm, -1), chr, utils::head(mm, -1), utils::tail(mm, -1), BEDcolor)
      DFlist[[paste0(chr, "_P")]] <- stats::setNames(DF, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))
    }

    ###
    # second_arm
    ###
    if (verbose) {
      message(chr, ": Q-arm")
    }
    # get start and end of arm
    S <- centChrom[centChrom[[1]] == chr][[3]]
    E <- max(experiment$IDX[experiment$IDX$V1 == chr][[3]])

    m <- NULL # m holds alls boundaries (in bps)
    if (!is.null(chunk)) { # check if chunks are given
      # set chunks
      steps <- seq(S, E, by = chunk / 2)

      tmpM <- list()
      for (STEP in steps) {
        hic.mat <- select.subset.hicseg(experiment, chr = chr, start = STEP, end = STEP + chunk)
        res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
        tmpM[[as.character(STEP)]] <- hic.mat$x[res$t_hat]
      }
      # keep boundaries that are found in two chunks
      m <- as.numeric(names(table(unlist(tmpM))[table(unlist(tmpM)) > 1]))
      # also keep boundaries that are found in the non-overlapping first part of the first window
      m <- c(m, tmpM[[1]][tmpM[[1]] < steps[2]])
    } else {
      hic.mat <- select.subset.hicseg(experiment, chr = chr, start = S, end = E)
      res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
      m <- hic.mat$x[res$t_hat]
    }

    DF <- NULL
    if (length(m) < 2) {
      warning(chr, "'s Q-arm has no borders")
    } else {
      DF <- data.frame(chr, utils::head(unique(m), -1), utils::tail(unique(m), -1), chr, utils::head(unique(m), -1), utils::tail(unique(m), -1), BEDcolor)
      DFlist[[paste0(chr, "_Q")]] <- stats::setNames(DF, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))
    }
  }
  DFlist <- do.call("rbind", DFlist)

  return(DFlist)
}
