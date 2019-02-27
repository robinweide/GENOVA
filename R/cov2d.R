outlier3Darray = function(ARRAY, Q = .995){
  ARRAY[is.na(ARRAY)] = 0
  tmp = apply(ARRAY, c(1,2), quantile, probs = Q, na.rm = TRUE)
  tmp.bk = tmp

  for(i in 1:nrow(tmp)){
    for(j in 1:ncol(tmp)){
      ARRAY[i,j,ARRAY[i, j , ] > tmp[i,j]] = tmp[i,j]
    }
  }

  return(ARRAY)
}

#' cov2d
#'
#' Get scores around given positions
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param minDist The maximal distance
#' @param add Shift pairs with a constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor).
#' @param outlierCutOff The severity of outliers: roughly translates to the amount of MADs above the median.
#' @param rmOutlier Try to perform outlier-correction
#' @param verbose Verbose?
#' @return A list containing a score-matrix and a count-variable of the amount of used bed-entries.
#' @import data.table
cov2d <- function( experiment, bed, minDist=5e6, maxDist = Inf, size=500e3, add=0,
                   outlierCutOff = 40, rmOutlier = F, verbose = T){
  #read data from experiment list
  data = experiment$ICE
  window = experiment$RES
  idx = experiment$ABS

  minDist <- minDist/window
  maxDist <- maxDist/window
  middle <- (bed[,2]+bed[,3])/2
  bed[,2] <- (middle-size)+add
  bed[,3] <- (middle+size)+add

  # wrap around pos, which exceed chromosome-length
  chromEnd                      <- max(idx[idx[,1] == unique(bed[,1]),3])
  bed[bed[,3] > chromEnd , 2:3] <- bed[bed[,3] > chromEnd , 2:3] - chromEnd


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
  rawMatList <- list()
  #double looping over the data


  if(length(start) < 2 ){
    return(NULL)
  }

  combos  <- t(combn(length(start), 2))
  combos[combos[,1] > combos[,2], 1:2] <- combos[combos[,1] > combos[,2], 2:1]
  combosD <- abs(start[combos[,1]] - start[combos[,2]])
  TF = combosD >= minDist  & combosD <= maxDist
  combos  <- combos[TF,]

  if(sum(TF) >= 2 ){
    for( idx in 1:nrow(combos)){
      #cat(i, "\r")

      i = combos[idx, 1]
      j = combos[idx, 2]

      if(start[i] > start[j]){
        i = combos[idx, 2]
        j = combos[idx, 1]
      }

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

  # Convert to 3D array
  if(length(rawMatList) > 0){
    rawMatList2 <- rawMatList[!unlist(lapply(rawMatList, is.null))]
    sm <- simplify2array(rawMatList2)

    #####################
    #  outlier correct  #
    #####################
    DATA = NULL
    if(rmOutlier){
      sm = outlier3Darray(ARRAY = sm, Q = outlierCutOff)
      DATA = apply(sm, c(1,2), sum, na.rm = T)
    } else {
      DATA = apply(sm, c(1,2), sum, na.rm = T)
    }

    return(list(score=DATA, count=count))
  } else {
    return(NULL)
  }
}
