select.subset.hicseg <- function(EXP, chrom, start, end) {
  sel <- EXP$ABS[EXP$ABS[, 1] == chrom & EXP$ABS[, 2] >= start & EXP$ABS[, 2] <= end, 4]
  start.i <- min(sel)
  end.i <- max(sel)
  # create matrix
  r1.s <- select.sub(EXP$ICE, start.i, end.i)
  r2.s <- r1.s
  sub.mat <- matrix(0, ncol = length(start.i:end.i), nrow = length(start.i:end.i))
  # fill the matrix with the first replicate
  index <- cbind(r1.s$V1 - start.i + 1, r1.s$V2 - start.i + 1)
  sub.mat[index] <- r1.s$V3
  # fill the matrix with the second replicate
  index <- cbind(r2.s$V2 - start.i + 1, r2.s$V1 - start.i + 1)
  sub.mat[index] <- r2.s$V3
  pos <- EXP$ABS[EXP$ABS[, 4] %in% (start.i:end.i), 3]
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

  DFlist <- list()
  for (chrom in chromsToUse) {
    if (verbose) {
      message(chrom, ": started")
    }
    # get the bed entry of the centromere
    centChrom <- experiment$CENTROMERES[experiment$CENTROMERES[, 1] == chrom, ]

    # throw a hissyfit if chomosome is not found
    if (length(centChrom) == 0) {
      message("There is no centromere-information for ", chrom, ".
Skipping this chromosome.")
      next()
    }

    # get chromosome-size
    chromSize <- max(experiment$ABS[experiment$ABS$V1 == chrom, 3])

    ###
    # first_arm
    ###
    if (verbose) {
      message(chrom, ": P-arm")
    }
    # get start and end of arm
    S <- 0
    E <- centChrom[centChrom[, 1] == chrom, 2]

    m <- NULL # m holds alls boundaries (in bps)
    if (!is.null(chunk)) { # check if chunks are given

      # set chunks
      steps <- seq(S, E, by = chunk / 2)

      tmpM <- list()
      for (STEP in steps) {
        hic.mat <- select.subset.hicseg(experiment, chrom = chrom, start = STEP, end = STEP + chunk)
        res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
        tmpM[[as.character(STEP)]] <- hic.mat$x[res$t_hat]
      }

      # keep boundaries that are found in two chunks
      m <- as.numeric(names(table(unlist(tmpM))[table(unlist(tmpM)) > 1]))
      # also keep boundaries that are found in the non-overlapping first part of the first window
      m <- c(m, tmpM[[1]][tmpM[[1]] < steps[2]])
    } else {
      hic.mat <- select.subset.hicseg(experiment, chrom = chrom, start = S, end = E)
      res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
      m <- hic.mat$x[res$t_hat]
    }

    DF <- NULL
    if (length(m) < 2) {
      warning(chrom, "'s P-arm has no borders")
    } else {
      DF <- data.frame(chrom, head(unique(m), -1), tail(unique(m), -1), chrom, head(unique(m), -1), tail(unique(m), -1), BEDcolor)
      DFlist[[paste0(chrom, "_P")]] <- setNames(DF, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))
    }

    ###
    # second_arm
    ###
    if (verbose) {
      message(chrom, ": Q-arm")
    }
    # get start and end of arm
    S <- centChrom[centChrom[, 1] == chrom, 3]
    E <- max(experiment$ABS[experiment$ABS$V1 == chrom, 3])

    m <- NULL # m holds alls boundaries (in bps)
    if (!is.null(chunk)) { # check if chunks are given
      # set chunks
      steps <- seq(S, E, by = chunk / 2)

      tmpM <- list()
      for (STEP in steps) {
        hic.mat <- select.subset.hicseg(experiment, chrom = chrom, start = STEP, end = STEP + chunk)
        res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
        tmpM[[as.character(STEP)]] <- hic.mat$x[res$t_hat]
      }
      # keep boundaries that are found in two chunks
      m <- as.numeric(names(table(unlist(tmpM))[table(unlist(tmpM)) > 1]))
      # also keep boundaries that are found in the non-overlapping first part of the first window
      m <- c(m, tmpM[[1]][tmpM[[1]] < steps[2]])
    } else {
      hic.mat <- select.subset.hicseg(experiment, chrom = chrom, start = S, end = E)
      res <- HiCseg::HiCseg_linkC_R(nrow(hic.mat$z), floor(nrow(hic.mat$z) / binPerBorder), "P", hic.mat$z, "Dplus")
      m <- hic.mat$x[res$t_hat]
    }

    DF <- NULL
    if (length(m) < 2) {
      warning(chrom, "'s Q-arm has no borders")
    } else {
      DF <- data.frame(chrom, head(unique(m), -1), tail(unique(m), -1), chrom, head(unique(m), -1), tail(unique(m), -1), BEDcolor)
      DFlist[[paste0(chrom, "_Q")]] <- setNames(DF, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))
    }
  }
  DFlist <- do.call("rbind", DFlist)

  return(DFlist)
}
