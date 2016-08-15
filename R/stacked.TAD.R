#' Get a z-stack matrix of TADs
#'
#' Extracts matrices from a BED-like structure and resizes them, which leads to all start- and end-position of TADs overlapping.
#' Sums over all matrices to produce a single Z-stack matrix, which is normalised to 100 loops.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param tad.bed Bedpe file containing the TAD positions: produced by read.bedpe(). Please make a bedpe, by copying tad-bed columns to columns 4,5 and 6!
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @return A matrix containing the Z-stack scores.
#' @export
stacked.TAD <- function(experiment, tad.bed, smallTreshold = 225e3, verbose = F){
  hicdata <- experiment$ICE
  bed <- experiment$ABS
  resolution <- experiment$RES
  # Check for setkey
  if(is.null(data.table::key(hicdata))){data.table::setkey(hicdata, V1, V2)}
  # Make chr:pos index of HiC-index
  bed.p <- paste0(bed[,1], ":", bed[,2])
  # Remove smaller tads
  tad.bed <- tad.bed[abs(tad.bed[,6]-tad.bed[,2]) >= smallTreshold ,]
  # Is BED 1 upstream of BED2?
  for(i in 1:length(tad.bed[,1])){
    if( tad.bed[i,2] > tad.bed[i,5]){
      tad.bed <- tad.bed[i,c(4,5,6,1,2,3)]
    }
  }
  # Prune NA-rows
  tad.bed <- na.exclude(tad.bed[1:6])
  # Make a bed from bedpe, with cols 1,2 and 6
  tad <- tad.bed[,c(1,2,6)]
  # Loop trough tads and sum over scorematrices
  tad.length <- length(tad$V1)
  # Add missing levels to Chromosomes
  levels(tad$V1) <- levels(bed$V1)
  # Successfull tads:
  SL <- 0
  for( i in 1:tad.length){
    if(verbose){cat(paste0(i, ' of ', tad.length, ' tads.'), "\r")}
    # Determine size of tad
    tadSize <-  tad[i,3] - tad[i,2]
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
    if(!exists("results.vector")){
      results.vector <- sel.resized
    } else {
      results.vector <- results.vector + sel.resized
    }
    SL <- SL + 1
  }
  return( (results.vector/SL)[1:99,1:99])
}
