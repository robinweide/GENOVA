#' Get a z-stack matrix of TADs
#'
#' Extracts matrices from a BED-like structure and resizes them, which leads to all start- and end-position of TADs overlapping.
#' Sums over all matrices to produce a single Z-stack matrix, which is normalised to 100 loops.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param tad.bed Data.frame from a Bed file containing the TAD positions.
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @param saveRawList Logical: True will output the raw matrices per TAD.
#' @param saveRaw Logical: True will output an aditional matrix without outlier-correction.
#' @param outlierCutOff The severity of outliers: roughly translates to the amount of MADs above the median.
#' @return A list containing a matrix with the Z-stack scores (`STACK`). Optionally, it can contain the raw matrices per TAD (`STACK.list`) and a matrix without the outlier-correction (`STACK.raw`).
stackR.old <- function (experiment, tad.bed, smallTreshold = 225000, verbose = F,
    saveRaw = T, saveRawList = T,outlierCutOff = 8){
		if(any(tad.bed[,2] > tad.bed[,3])){
			warning("5' TAD border larger then 3' TAD border for some entries")
		}
    MADTRESHOLD <- outlierCutOff
    rawMatList = list()
    hicdata <- experiment$ICE
    bed <- experiment$ABS
    bed.p <- paste0(bed[,1], ":", bed[,2])
    bed[,1] <- as.factor(bed[, 1])
    resolution <- experiment$RES
    if (is.null(data.table::key(hicdata))) {
        data.table::setkey(hicdata, V1, V2)
    }
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
        ##
        if(verbose){cat(paste0(i, ' of ', tad.length, ' tads.'), "\r")}
        # Determine size of tad
        tadSize <-  tad[i,3] - tad[i,2]
        #Elzo you removed them on line 21
        if(tadSize < smallTreshold){next}
        halftadSize <- floor((tadSize/resolution)/2)*resolution
        # Get start positions of the tad, floored to 10kb windows, and make chr:pos index
        start <- paste0(tad[i,1], ":", as.integer(resolution*floor(  (tad[i,2] - halftadSize  )/resolution)))
        end <- paste0(tad[i,1], ":", as.integer(resolution*floor(  (tad[i,3] + halftadSize   )/resolution)))
        # Bugfix: in some cases paste0 yields an extra whitespace
        start <- gsub(" ", "", start, fixed = TRUE)
        end <- gsub(" ", "", end, fixed = TRUE)
        # Get HiC-pro indexes of start and end of virt4C-plot
        start.pos <- bed[match(start,bed.p),4]
        end.pos <- bed[match(end,bed.p),4]
        # Check if valid indeces (Missing or Mitochondrial stuff)
        if(!all(is.finite(start.pos),is.finite(end.pos))){next}
        # Extract data from HiC-pro matrix
        #Elzo why not select.subset?
        x <- rep(start.pos:end.pos, each=length(start.pos:end.pos))
        y <- rep((start.pos):(end.pos), length(start.pos:end.pos))
        sel <- hicdata[list(x,y)]
        sel$V3[is.na(sel$V3)] <- 0
        # Skip low-coverage sites
        if(nrow(sel) < 10){next}
        # Fastest way!
        newMat <- reshape2::acast(sel, V1~V2, value.var="V3")
        newMat[lower.tri(newMat)] <- t(newMat)[lower.tri(newMat)]
        # Rescale to 100 breaks
        sel.resized <- resize.mat(newMat, c(100,100))
        ##
        rawMatList[[i]] <- sel.resized
        if (!exists("results.vector")) {
            results.vector <- sel.resized
        }
        else {
            results.vector <- results.vector + sel.resized
        }
        SL <- SL + 1
    }

      rawMatList <- rawMatList[!unlist(lapply(rawMatList, is.null))]
      sm <- simplify2array(rawMatList)
      MED <- apply(sm, MARGIN = 1:2, median)
      MAD <- apply(sm, MARGIN = 1:2, mad)
      tres <- MED + (MAD * MADTRESHOLD)
      tres[is.na(tres)] <- 0
      for (i in 1:length(rawMatList)) {
          m <- rawMatList[[i]]
          if (any(m[1:99, 1:99] > tres[1:99, 1:99] * 5)) {
              rawMatList[[i]] <- matrix(0, nrow = 100, ncol = 100)
              SL - 1
          }
      }
      STACK.outlierCorrected <- Reduce(rawMatList, f = "+")
      STACK.rawMatList <- rawMatList
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
