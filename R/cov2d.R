#' cov2d
#'
#' Get scores around given positions
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param add Add constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor). 
#' @return A list containing a score-matrix and a count-variable of the amount of used bed-entries.
#' @import data.table
cov2d <- function( experiment, bed, minDist=5e6, size=500e3, add=0){
  #read data from experiment list
  data = experiment$ICE
  window = experiment$RES 
  idx = experiment$ABS
  
  minDist <- minDist/window
  middle <- (bed[,2]+bed[,3])/2
  bed[,2] <- middle-size+add
  bed[,3] <- middle+size+add
  
  p1 <- paste0(bed[,1], ":", as.integer(window*floor(bed[,2]/window)) )
  p2 <- paste0(bed[,1], ":", as.integer(window*floor(bed[,3]/window)) )
  idx.p <- paste0(idx[,1], ":", idx[,2])
  start <- idx[match(p1,idx.p),4]
  end   <- idx[match(p2,idx.p),4]
  
  #make sure the bed positions are found in the genome/Hi-C pro file
  sel <- !is.na(start) & !is.na(end)
  bed <- bed[sel,]
  start <- start[sel]
  end <- end[sel]
  score.matrix <- matrix(0, ncol=max(end-start)+1, nrow=max(end-start)+1)
  count <- 0
  #double looping over the data
  for( i in 1:(length(start)-1)){
    #cat(i, "\r")
    for( j in (i+1):length(start)){
      if(start[j]-start[i] < minDist)
        next
      
      x <- rep(start[i]:end[i],end[j]-start[j]+1)
      y <- rep(start[j]:end[j],each=end[i]-start[i]+1)
      
      s <- data[list(x,y)]
      
      #NOTE! this is currently hard-coded, how can we make this a dynamic value based
      #on the number of reads we normalize to (this should then be stored in the experiment list)
      s <- s[!is.na(s$V3) & s$V3 < 10]
      
      #create a sub matrix to add to the actual matrix
      score.matrix[cbind(s$V1-start[i]+1,s$V2-start[j]+1)] <- score.matrix[cbind(s$V1-start[i]+1,s$V2-start[j]+1)] + s$V3
      count <- count+1
    }       
    
  }
  list(score=score.matrix, count=count)
}
