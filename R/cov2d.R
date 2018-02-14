#' cov2d
#'
#' Get scores around given positions
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param add Add constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor).
#' @param outlierCutOff The severity of outliers: roughly translates to the amount of MADs above the median.
#' @param rmOutlier Try to perform outlier-correction
#' @return A list containing a score-matrix and a count-variable of the amount of used bed-entries.
#' @import data.table
cov2d <- function( experiment, bed, minDist=5e6, size=500e3, add=0, outlierCutOff = 40, rmOutlier = F){
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
  # cat(sel,'\n')
  bed <- bed[sel,]
  start <- start[sel]
  end <- end[sel]
  score.matrix <- matrix(0, ncol=max(end-start)+1, nrow=max(end-start)+1)
  count <- 0
  rawMatList <- list()
  #double looping over the data
  for( i in 1:(length(start)-1)){
    #cat(i, "\r")
    for( j in (i+1):length(start)){
      if(start[j]-start[i] < minDist)
        next

      x <- rep(start[i]:end[i],end[j]-start[j]+1)
      y <- rep(start[j]:end[j],each=end[i]-start[i]+1)

      s <- data[list(x,y)]
      s.mat <- reshape2::acast(s, V1~V2, value.var="V3")

      # Remove NA
      s.mat[is.na(s.mat)] <- 0

      # Add to rawMatList
      rawMatList[[paste(i,j,sep = "_")]] <- s.mat

      count <- count+1
    }

  }

  if(length(rawMatList) < 10){
    warning(paste0("There are not a lot of comparisons for chromosome ",as.character(unique(bed[,1])),". Please tred carefully!\n"))
  }
  rawMatList <- rawMatList[!unlist(lapply(rawMatList, is.null))]
  if(rmOutlier){
    sm <- simplify2array(rawMatList)
    sm[sm == 0] <- NA
    MED <- apply(sm,MARGIN = 1:2, FUN = function(x) median(x,na.rm = T))
    MAD <- apply(sm,MARGIN = 1:2, FUN = function(x) mad(x,na.rm = T))
    tres <- MED+(MAD*outlierCutOff)
    tres[is.na(tres)] <- 0
    for(i in 1:length(rawMatList)){
      #rawMatList[[i]][rawMatList[[i]][1:99,1:99] > tres[1:99,1:99]] <- 0 #tres[rawMatList[[i]][1:99,1:99] > tres[1:99,1:99]]
      m <- rawMatList[[i]]
      #cat(i, "\n")
      if(any(m > tres)){
        rawMatList[[i]] <- matrix(0, nrow = dim(sm)[1], ncol=dim(sm)[2])
        count <- count - 1
      }
    }
    STACKoutlierr <- Reduce(rawMatList, f = '+')
  } else {
    STACKoutlierr <- Reduce(rawMatList, f = '+')
  }
  return(list(score=STACKoutlierr, count=count))
}
