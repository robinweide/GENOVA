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
#' @export
stackR <- function (experiment, tad.bed, smallTreshold = 225000, verbose = F,
    saveRaw = T, saveRawList = T,outlierCutOff = 8){
		if(any(tad[,2] > tad[,5])){
			stop("5' TAD border larger then 3' TAD border for some entries")
		}
    MADTRESHOLD <- outlierCutOff
    rawMatList = list()
    hicdata <- experiment$ICE
    bed <- experiment$ABS
    resolution <- experiment$RES
    if (is.null(data.table::key(hicdata))) {
        data.table::setkey(hicdata, V1, V2)
    }
    tad.bed <- tad.bed[abs(tad.bed[, 6] - tad.bed[, 2]) >= smallTreshold,
        ]
    tad.bed <- na.exclude(tad.bed[1:6])
    tad <- tad.bed[, c(1, 2, 6)]
    tad.length <- length(tad[, 1])
    levels(tad$V1) <- levels(bed[, 1])
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
				newMat <- select.subset( hicdata, tad[i,1], tad[i,2] - halftadSize, tad[i,3] + halftadSize, bed )

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
    if (saveRaw) {
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
        STACKoutlierr <- Reduce(rawMatList, f = "+")
        return(list(STACK = (results.vector/SL)[1:99, 1:99],
            RAW = rawMatList, STACKoutlier = (STACKoutlierr/SL)[1:99,
                1:99]))
    }
    else {
        return((results.vector/SL)[1:99, 1:99])
    }
}
