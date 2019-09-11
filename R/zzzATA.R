#' Aggegrate TAD Analysis
#'
#' Get a z-stack matrix of TADs: this function extracts matrices from a BED-like structure and resizes them,
#' which leads to all start- and end-position of TADs overlapping. Next, it averages all matrices to produce
#' a single Z-stack matrix.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param experiment The Hi-C experiment object of a sample: produced by \code{construct.experiment()}.
#' @param tad.bed Data.frame from a BED-file containing the TAD positions.
#' @param smallTreshold The minimal size of TADs. Too small TADs lead to messy plots.
#' @param rmOutlier Perform outlier-correction by tresholding [outlierCutOff] percentile per pixel.
#' @param verbose Produces a progress-indication.
#' @param outlierCutOff The severity of outliers. We compute the [outlierCutOff] percentile per pixel and set values bigger than that to this value.
#' @return \item{STACK}{a list containing a matrix with the Z-stack scores}
#' @return \item{STACK.raw}{a matrix without the outlier-correction}
#' @return \item{STACK.list}{raw matrices per TAD}
#' @return \item{OUTLIERCORRECTIONSWITCH}{the outlierCutOff if rmOutlier is TRUE, otherwise set to FALSE}
#' @examples
#' ## Not run:
#' ATA_results_of_RAOetal <- ATA(experiment = Rao_20k, tad.bed = TADs)
#' ## End(**Not run**)
#' @export
ATA_old <- function(experiment, tad.bed, smallTreshold = 225000, rmOutlier = F, verbose = F, outlierCutOff = .995) {
  OUTLIERCORRECTIONSWITCH <- ifelse(rmOutlier,
    yes = T,
    no = outlierCutOff
  )

  if (any(tad.bed[, 2] > tad.bed[, 3])) {
    stop("5' TAD border downstream then 3' TAD border for some entries")
  }
  rawMatList <- list()

  bed <- experiment$ABS

  bed[, 1] <- as.factor(bed[, 1])
  resolution <- experiment$RES
  if (is.null(data.table::key(experiment$ICE))) {
    data.table::setkey(experiment$ICE, V1, V2)
  }

  tad.bed <- tad.bed[abs(tad.bed[, 3] - tad.bed[, 2]) >= smallTreshold, ]
  tad.bed <- na.exclude(tad.bed[1:3])
  tad <- tad.bed[, c(1, 2, 3)]
  tad.length <- length(tad[, 1])
  tad[, 1] <- factor(tad[, 1], levels = levels(bed[, 1]))
  outputTAD <- tad


  #####################
  #       vapply      #
  #####################
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(1)

  tad$size <- as.numeric(abs(tad[, 3] - tad[, 2]))
  tad$halftadSize <- round(tad$size / 2)
  tad$regionStart <- tad[, 2] - tad$halftadSize
  tad$regionEnd <- tad[, 3] + tad$halftadSize


  tad.list <- split(tad, seq(nrow(tad)))

  rawMatList <- vapply(tad.list, FUN.VALUE = matrix(NA_real_, nrow = 100, ncol = 100), FUN = function(x) {
    newMat <- select.subset(
      exp = experiment,
      chrom = x[, 1],
      start = x[, 6],
      end = x[, 7]
    )
    resize.mat(newMat$z, c(100, 100))
  })

  data.table::setDTthreads(old_threads)

  # Convert to 3D arra
  sm <- rawMatList # simplify2array(rawMatList)

  #####################
  #  outlier correct  #
  #####################
  if (rmOutlier) {
    sm.bk <- sm
    sm <- outlier3Darray(ARRAY = sm, Q = outlierCutOff)

    ATAraw <- apply(sm.bk, c(1, 2), mean)
    ATA <- apply(sm, c(1, 2), mean)
  } else {
    ATAraw <- apply(sm, c(1, 2), mean)
    ATA <- apply(sm, c(1, 2), mean)
  }


  ##########################
  #       so long and      #
  #  thanks for the fish!  #
  ##########################

  return(list(
    STACK = ATA[1:99, 1:99],
    STACK.raw = ATAraw[1:99, 1:99],
    STACK.list = rawMatList,
    OUTLIERCORRECTIONSWITCH = OUTLIERCORRECTIONSWITCH,
    TADs.bed = outputTAD
  ))
}

outlier3Darray <- function(ARRAY, Q = .995) {
  ARRAY[is.na(ARRAY)] <- 0
  tmp <- apply(ARRAY, c(1, 2), quantile, probs = Q, na.rm = TRUE)
  tmp.bk <- tmp

  for (i in 1:nrow(tmp)) {
    for (j in 1:ncol(tmp)) {
      ARRAY[i, j, ARRAY[i, j, ] > tmp[i, j]] <- tmp[i, j]
    }
  }

  return(ARRAY)
}
