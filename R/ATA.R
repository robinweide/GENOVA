#' Get a z-stack matrix of TADs
#'
#' Extracts matrices from a BED-like structure and resizes them, which leads to all start- and end-position of TADs overlapping.
#' Sums over all matrices to produce a single Z-stack matrix, which is normalised to 100 loops.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param experiment The Hi-C experiment object of a sample: produced by \code{construct.experiment()}.
#' @param tad.bed Data.frame from a Bed file containing the TAD positions.
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @param saveRawList Logical: True will output the raw matrices per TAD.
#' @param saveRaw Logical: True will output an aditional matrix without outlier-correction.
#' @param outlierCutOff The severity of outliers: roughly translates to the amount of MADs above the median.
#' @return \item{STACK}{A list containing a matrix with the Z-stack scores}
#' @return \item{STACK.list}{Optional: raw matrices per TAD}
#' @return \item{STACK.raw}{Optional: a matrix without the outlier-correction}
#' @examples
## Not run:
#' ATA_results_of_RAOetal <- ATA(experiment = Rao_20k,tad.bed = TADs)
## End(**Not run**)
#' @export
ATA <- function (experiment, tad.bed, smallTreshold = 225000, verbose = F,
    saveRaw = T, saveRawList = T,outlierCutOff = 40){
		if(any(tad.bed[,2] > tad.bed[,3])){
			warning("5' TAD border downstream then 3' TAD border for some entries")
		}
    MADTRESHOLD <- outlierCutOff
    rawMatList = list()

    bed <- experiment$ABS

    bed[,1] <- as.factor(bed[, 1])
    resolution <- experiment$RES
    if (is.null(data.table::key(experiment$ICE))) {
        data.table::setkey(experiment$ICE, V1, V2)
    }
    # issue 7: stacked.TAD.R : TAD bed file
    tad.bed <- tad.bed[abs(tad.bed[, 3] - tad.bed[, 2]) >= smallTreshold,]
    tad.bed <- na.exclude(tad.bed[1:3])
    tad <- tad.bed[, c(1, 2, 3)]
    tad.length <- length(tad[, 1])
    tad[,1] <- factor(tad[,1], levels=levels(bed[,1]))

    SL <- 0
    for (i in 1:tad.length) {
        tadSize <- tad[i, 3] - tad[i, 2]
        if (tadSize < smallTreshold) {
            next
        }
        if (verbose) {
            cat(i, " of ", tad.length, " tads.", "\r")
        }
				#select TAD data from select.subset
        halftadSize <- floor((tadSize/resolution)/2) * resolution
				newMat <- select.subset( experiment, tad[i,1], tad[i,2] - halftadSize, tad[i,3] + halftadSize)

        sel.resized <- resize.mat(newMat$z, c(100, 100))
        rawMatList[[i]] <- sel.resized
        if (!exists("results.vector")) {
            results.vector <- sel.resized
        }
        else {
            results.vector <- results.vector + sel.resized
        }
        SL <- SL + 1
    }
    #rawMatListCopy <- rawMatList


      #rawMatList <- rawMatList[!unlist(lapply(rawMatList, is.null))]
      sm <- simplify2array(rawMatList)
      MED <- apply(sm, MARGIN = 1:2, median)
      MAD <- apply(sm, MARGIN = 1:2, mad)
      tres <- MED + (MAD * MADTRESHOLD)
      tres[is.na(tres)] <- 0
      for (i in 1:length(rawMatList)) {
          m <- rawMatList[[i]]

          # get boolMatrix of sites in m bigger than tres:
          biggerThanTres <- m[1:99, 1:99] > tres[1:99, 1:99]
          zeroesInTres <- tres[1:99, 1:99] > 0
          combinedBool <- biggerThanTres + zeroesInTres

          if(!all(ifelse(combinedBool > 1, F, T))) {
              rawMatList[[i]] <- matrix(0, nrow = 100, ncol = 100)
              SL = SL - 1
          }
      }
      STACK.outlierCorrected <- Reduce(as.list(rawMatList), f = "+")
      STACK.rawMatList <- rawMatList

      # check is tres has zeroes. If so, data is too low-qual to do outlier-correction with current set of loops
      if(any(tres == 0)){
        warning("\nThe data is too sparse to do outlier-correction\n\twith current set of TADs.\nOutput will be without outlier-correction")
        STACK.outlierCorrected = results.vector
      }


      if(saveRaw){
        if(saveRawList){
            return(list(STACK = (STACK.outlierCorrected/SL)[1:99, 1:99],
                      STACK.raw = (results.vector/SL)[1:99, 1:99],
                      STACK.list = STACK.rawMatList))
        } else {
            return(list(STACK = (STACK.outlierCorrected/SL)[1:99, 1:99],
                      STACK.raw = (results.vector/SL)[1:99, 1:99]))
          }
      } else{
        if(saveRawList){
          return(list(STACK = (STACK.outlierCorrected/SL)[1:99, 1:99],
                      STACK.list = STACK.rawMatList))
        } else {
          return(list(STACK = (STACK.outlierCorrected/SL)[1:99, 1:99]))
        }
      }
  }
