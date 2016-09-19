#' Get a z-stack matrix of 2D-regions (e.g. loops).
#'
#' Extract matrices around a defined set of pixels, like possible loops.
#' Sums over all matrices to produce a single Z-stack matrix.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param loop.bed Bedpe file containing the loop positions: produced by read.bedpe().
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor). 
#' @param saveRaw Logical: True will output the raw matrices per loop and performs outlier-detection.
#' @param outlierCutOff The severity of outliers: roughly translates to the amount of MADs above the median.
#' @return A list of a matrix containing the Z-stack scores (APA), the raw matrices (rawMatList) and the outlier-removed matrix containing the Z-stack scores (APAoutlier).
#' @import data.table
#' @export
APA <- function(experiment, loop.bed, smallTreshold = 225e3, size = 21, verbose = F, saveRaw = T, outlierCutOff = 8, ...){
  MADTRESHOLD <- outlierCutOff
  if(((size-1) /2 )%%2 != 0){stop("Size should be an even number +1")}
  size.offset = (size-1)/2
  if(saveRaw){
    rawMatList = list()
  }
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
  #Elzo maybe do something like this
  loop.bed[loop.bed[,2] > loop.bed[,5],] <- loop.bed[loop.bed[,2] > loop.bed[,5],c(4,5,6,1,2,3)]
  
  for(i in 1:length(loop.bed[,1])){
    if( loop.bed[i,2] > loop.bed[i,5]){
      loop.bed <- loop.bed[i,c(4,5,6,1,2,3)]
    }
  }
  # Prune NA-rows
  #Elzo comma before 1:6
  loop.bed <- na.exclude(loop.bed[,1:6])
  #loop.bed <- na.exclude(loop.bed[1:6])
  # Get start positions of the loop, floored to 10kb windows, and make chr:pos index
  loop.bed1.p <- paste0(loop.bed[,1], ":", as.integer(resolution*floor(  ( (abs( loop.bed[,3] - loop.bed[,2] ) /2 )+ loop.bed[,2])   /resolution)))
  loop.bed2.p <- paste0(loop.bed[,4], ":", as.integer(resolution*floor(  ( (abs( loop.bed[,6] - loop.bed[,5] )/2)  + loop.bed[,5])   /resolution)))
  
  #Elzo why not the following
  #just average start and end (middle), divide by the resolution, floor and multiply by the resolution
  #Robin: because that doesn't work: *e*-ids throughout
  # loop.bed1.p <- paste0(loop.bed[,1], ":", resolution*floor( ( (loop.bed[,2]+loop.bed[,3])/2)/resolution) )
  # loop.bed2.p <- paste0(loop.bed[,4], ":", resolution*floor( ( (loop.bed[,5]+loop.bed[,6])/2)/resolution) )
  
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
  score.matrix <- matrix(0, ncol=size, nrow=size)
  # Loop trough loops and sum over scorematrices
  lx <- length(x.pos)
  for( i in 1:lx){
    if(verbose){cat(paste0(i, ' of ', lx, ' loops.'), "\r")}
    # Get Indeces of 10 up/downstream of anchor
    sel.x <- (x.pos[i]-size.offset):(x.pos[i]+size.offset)
    sel.y <- (y.pos[i]-size.offset):(y.pos[i]+size.offset)
    # Exract scores from HiC
    s <- select.sub.2D(hicdata, sel.x,sel.y)
    s$V1 <- (s$V1 - min(sel.x)) + 1
    s$V2 <- (s$V2 - min(sel.y)) + 1
    # Check for 21x21 mat:
    if(length(unique(s$V2)) < size){
      cols <- 1:size
      toADD <- cols[!cols %in% unique(s$V2)]
      newdf <- data.frame(unique(cbind(rep(1:size, each = size), rep(toADD, length = size),0)))
      colnames(newdf) <- c("V1","V2","V3")
      s <- rbind(s, newdf)
    }
    if(length(unique(s$V1)) < size){
      cols <- 1:size
      toADD <- cols[!cols %in% unique(s$V1)]
      newdf <- data.frame(unique(cbind(rep(toADD, length = size),rep(1:size, each = size),0)))
      colnames(newdf) <- c("V1","V2","V3")
      s <- rbind(s, newdf)
    }
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
    # Add to rawMatList
    if(saveRaw){
      rawMatList[[i]] <- s.mat
    }
  } 
  #Elzo you tested for the fact that all the loop.bed1.p are in bed.p right?
  #why the match, cant you do.
  norma_loopCounts <- (score.matrix/length(loop.bed1.p)) 
  SL <- length(loop.bed1.p)
  # Rotate 90CW, so that diagonal of HiC-matrix is in bottomleft
  norma_loopCounts <- t(apply(norma_loopCounts, 2, rev))
  colnames(norma_loopCounts) <- 1:size
  if(saveRaw){
    #return(list(STACK = (results.vector/SL)[1:99,1:99],RAW = simplify2array(rawMatList)))
    rawMatList <- rawMatList[!unlist(lapply(rawMatList, is.null))]
    sm <- simplify2array(rawMatList)
    MED <- apply(sm,MARGIN = 1:2, median)
    MAD <- apply(sm,MARGIN = 1:2, mad)
    tres <- MED+(MAD*MADTRESHOLD)
    tres[is.na(tres)] <- 0
    for(i in 1:length(rawMatList)){
      #rawMatList[[i]][rawMatList[[i]][1:99,1:99] > tres[1:99,1:99]] <- 0 #tres[rawMatList[[i]][1:99,1:99] > tres[1:99,1:99]]
      m <- rawMatList[[i]]
      #cat(i, "\n")
      if(any(m[1:size,1:size] > tres[1:size,1:size]*5)){
        rawMatList[[i]] <- matrix(0, nrow = size, ncol=size)
        SL - 1
      }
    }
    STACKoutlierr <- Reduce(rawMatList, f = '+')
    norma_loopCountss <- (STACKoutlierr/SL) 
    norma_loopCountss <- t(apply(norma_loopCountss, 2, rev))
    colnames(norma_loopCountss) <- 1:size
    
    return(list(APA = norma_loopCounts,rawMatList = rawMatList, APAoutlier =  norma_loopCountss   ))
    #return(list(APA = norma_loopCounts, rawMatList = rawMatList))
  } else {
    return(norma_loopCounts)
  }
}
