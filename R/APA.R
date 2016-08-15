#' Get a z-stack matrix of 2D-regions (e.g. loops).
#'
#' Extract matrices around a defined set of pixels, like possible loops.
#' Sums over all matrices to produce a single Z-stack matrix, which is normalised to 100 loops.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param loop.bed Bedpe file containing the loop positions: produced by read.bedpe().
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @return A tidy data.frame containing the Z-stack scores.
#' @import data.table
#' @export
APA <- function(experiment, loop.bed, smallTreshold = 225e3, verbose = F){
  hicdata <- experiment$ICE
  bed <- experiment$ABS
  resolution <- experiment$RES
  # Check for setkey
  if(is.null(data.table::key(hicdata))){data.table::setkey(hicdata, V1, V2)}
  # Make chr:pos index of HiC-index
  bed.p <- paste0(bed[,1], ":", bed[,2])
  # Remove smaller loops
  loop.bed <- loop.bed[abs(loop.bed[,6]-loop.bed[,2]) >= smallTreshold ,]
  # Is BED 1 upstream of BED2?
  for(i in 1:length(loop.bed[,1])){
    if( loop.bed[i,2] > loop.bed[i,5]){
      loop.bed <- loop.bed[i,c(4,5,6,1,2,3)]
    }
  }
  # Prune NA-rows
  loop.bed <- na.exclude(loop.bed[1:6])
  # Get start positions of the loop, floored to 10kb windows, and make chr:pos index
  loop.bed1.p <- paste0(loop.bed[,1], ":", as.integer(resolution*floor(  ( (abs( loop.bed[,3] - loop.bed[,2] ) /2 )+ loop.bed[,2])   /resolution)))
  loop.bed2.p <- paste0(loop.bed[,4], ":", as.integer(resolution*floor(  ( (abs( loop.bed[,6] - loop.bed[,5] )/2)  + loop.bed[,5])   /resolution)))
  # Bugfix: in some cases paste0 yields an extra whitespace
  loop.bed1.p <- gsub(" ", "", loop.bed1.p, fixed = TRUE)
  loop.bed2.p <- gsub(" ", "", loop.bed2.p, fixed = TRUE)
  # Check if loops are all found:
  if(!all(loop.bed1.p %in% bed.p)){
    if(!all(loop.bed2.p %in% bed.p) ){
      stop('not all loop-anchors can not be found in HiC-file.\n# Are you sure that both are from the same complete reference?')
    }
    else{
      stop(paste0(table(foo %in% bar)[1],' upstream loop-anchors can not be found in HiC-file.\n# Are you sure that both are from the same complete reference?'))}
  }else{
    if(!all(loop.bed2.p %in% bed.p) ){
      stop(paste0(table(foo %in% bar)[1],' downstream loop-anchors can not be found in HiC-file.\n# Are you sure that both are from the same complete reference?'))
    }
  }
  # Get anchor HiC-indexes
  x.pos <- bed[match(loop.bed1.p,bed.p),4]
  y.pos <- bed[match(loop.bed2.p,bed.p),4]
  # Initialise matrix
  score.matrix <- matrix(0, ncol=21, nrow=21)
  # Loop trough loops and sum over scorematrices
  lx <- length(x.pos)
  for( i in 1:lx){
    if(verbose){cat(paste0(i, ' of ', lx, ' loops.'), "\r")}
    # Get Indeces of 10 up/downstream of anchor
    sel.x <- (x.pos[i]-10):(x.pos[i]+10)
    sel.y <- (y.pos[i]-10):(y.pos[i]+10)
    # Exract scores from HiC
    s <- select.sub.2D(hicdata, sel.x,sel.y)
    s$V1 <- (s$V1 - min(sel.x)) + 1
    s$V2 <- (s$V2 - min(sel.y)) + 1
    # Make matrix
    # Fastest way!
    s.mat <- reshape2::acast(s, V1~V2, value.var="V3")
    # # Old way
    # s.mat <- matrix(nrow=21, ncol=21,
    #                 dimnames=list(1:21, 1:21))
    # s.mat[cbind(s$V1, s$V2)] <- s$V3
    s.mat[is.na(s.mat)] <- 0
    # Add to Z-stack
    score.matrix <- score.matrix + s.mat
  }
  # Normalise to 100 loops
  norma_loopCounts <- (score.matrix/as.numeric(unname(table(loop.bed1.p %in% bed.p)) ) ) * 100
  # Rotate 90CW, so that diagonal of HiC-matrix is in bottomleft
  norma_loopCounts <- t(apply(norma_loopCounts, 2, rev))
  colnames(norma_loopCounts) <- 1:21
  return(norma_loopCounts)
}
