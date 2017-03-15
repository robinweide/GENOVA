#' Get a per-chromosome scaling plot
#'
#' Produce a dataframe with the probabilities of contacts in a set of distance-bins.
#' Bins are created on a log scale, which leads to equal amounts of datapoints per bin.
#'
#' @param experiment List of experiment-objects from `construct.experiment()`.
#' @param chromsToUse A vector containing the chromosome-names of interest.
#' @param maxDistance The maximal distance to calculate scalings for.
#' @param scaling Generate scaling instead of RCP?
#' @param verbose Produces a progress-indication.
#' @return A data_frame with distance-bin and probabilities.
#' @export
RCP <- function(experimentList, chromsToUse = NULL, maxDistance = 5e08, verbose = F,scaling = T){
  #find the largest stretch of 0s, which is
  #most likely the centromere
  largest.stretch <- function( x ){
    temp <- cumsum(c(1,diff(x) - 1))
    temp2 <- rle(temp)
    x[which(temp == with(temp2, values[which.max(lengths)]))]
  }

  #select a matrix of interactions for between two chromosomes
  selectTransData <- function(exp, chrom1, chrom2){
    bed <- exp$ABS
    data <- exp$ICE
    X <- bed[bed[,1]==chrom1,4]
    Y <- bed[bed[,1]==chrom2,4]
    #the order of the chromosomes matters for the analysis
    #make sure that X is smaller than Y, otherwise switch
    #them around
    if(X[1] > Y[1]){
      temp <- X
      X <- Y
      Y <- temp
      temp <- chrom1; chrom1 <- chrom2; chrom2 <- temp; #switch the chromosomes around as well
    }
    #create x and y vectors that contain the positions of the
    #entries in the matrix that we are creating
    x <- rep(X[1]:X[length(X)],        tail(Y, n=1) - Y[1] + 1)
    y <- rep(Y[1]:Y[length(Y)], each = tail(X, n=1) - X[1] + 1)
    data.sub <- data[base::list(x, y)]
    data.sub <- data.sub[!is.na(data.sub$V3)]
    #create an empty matrix, that has as many rows as the 'X' chromosome has
    #windows and as many columns as the 'Y' chromosome has windows
    mat <- matrix(0, ncol=tail(Y, n=1) - Y[1] + 1, nrow=tail(X, n=1) - X[1] + 1)
    mat[cbind(data.sub$V1-min(X)+1, data.sub$V2-min(Y)+1)] <- data.sub$V3
    x.pos <- bed[bed[,1]==chrom1,2]
    y.pos <- bed[bed[,1]==chrom2,2]
    #create a list that is compatible with the image function
    mat <- list(x=x.pos, y=y.pos, z=mat)
    mat
  }


  amountOfSamples <- length(experimentList)
  exp.names <- c()
  if(is.null(chromsToUse)){
    chromsToUse <-  unique(experimentList[[1]]$ABS[,1])
  }

  # Check for N chromosomes: more than 100: please specify chroms!
  if(length(chromsToUse) > 75){stop("Please restrict the amount of chromosomes with chromsToUse")}


  #check whether experiment names have been declared uniquely
  #otherwise use standard names for RCP
  for( i in 1:length(experimentList)){
    exp.names <- c(exp.names, experimentList[[i]]$name)
  }
  standard <- FALSE
  if(length(exp.names) != length(unique(exp.names))){
    warning("Experiment names have not been declared uniquely, using standard names")
    standard <- TRUE
  }

  d <- dplyr::data_frame(distance = integer(),
                  prob = numeric(),
                  sample = integer())
  dat <- NULL
  for(Ci in 1:length(chromsToUse)){
    chrom <- chromsToUse[Ci]
    for(i in 1:amountOfSamples){
      if(verbose){cat(paste0('Chromosome ', Ci,' of ', length(chromsToUse), ' chromosomes\r'))}
      dat <- NULL
      resolu <- experimentList[[i]]$RES
      BED_ABS <- data.table::data.table(experimentList[[i]]$ABS)
      data.table::setkey(BED_ABS, V1)
      idx1 <-BED_ABS[list(chrom)]

      if(chrom %in% names(experimentList[[i]]$CENTROMERES)){

        centromerixIDXes <- experimentList[[i]]$CENTROMERES[[chrom]]
        idx1 <- idx1[!V4 %in% centromerixIDXes]

      } else {
        # find biggest chrom to compare with
        allChroms <- table(data.table::data.table(experimentList[[i]]$ABS)$V1)
        if(length(allChroms) > 1){ # is there are more than 1 chromosome in the data, remove centromeres

          sizeChrom <- allChroms[names(allChroms)== chrom]
          allChroms <- allChroms[!names(allChroms) == chrom] # not the same chrom!
          toCompare <- names(which(allChroms == max(allChroms)))
          siztoCompare <- allChroms[names(allChroms)== toCompare]

          chromTheBiggest <-siztoCompare < sizeChrom
          trans.mat <-NULL
          if(chromTheBiggest){
            trans.mat <- selectTransData(experimentList[[i]], chrom1 = chrom , chrom2 = toCompare)
          } else {
            trans.mat <- selectTransData(experimentList[[i]], chrom2 = chrom , chrom1 =  toCompare)
          }

          #set the q.top highest values to this value
          th <- quantile(trans.mat$z, 1-(1e-5))
          trans.mat$z[trans.mat$z > th] <- th

          #get the centromere positions emprically
          centromerixIDXes <- NULL
          if(chromTheBiggest){ # chrom is chrom1
            cent1 <- which(apply(trans.mat$z,1,sum)==0)
            centromerixIDXes <- largest.stretch(cent1)
          } else {
            cent2 <- which(apply(trans.mat$z,2,sum)==0)
            centromerixIDXes <- largest.stretch(cent2)
          }
          experimentList[[i]]$CENTROMERES[[chrom]] <- centromerixIDXes
          idx1 <- idx1[!V4 %in% centromerixIDXes]
        }
      }




      x <- rep(idx1$V4, length(idx1$V4))
      y <- rep(idx1$V4, each=length(idx1$V4))
      xydf <- dplyr::data_frame(x,y)
      xydf.b <- dplyr::filter(xydf,abs(x-y) <= (maxDistance/resolu))

      x <- as.numeric(xydf$x)
      y <- as.numeric(xydf$y)

      breaks <- 10**seq(4,log(maxDistance, base = 10),length.out=81)

      m.chrom <- experimentList[[i]]$ICE[list(x,y)]

      if(nrow(m.chrom) == 0){ next }
      if(scaling == T){
        m.chrom <- dplyr::filter(m.chrom, V1 < V2)

        m.chrom$V3[is.na(m.chrom$V3)] <- 0


        distance = resolu*(m.chrom$V2-m.chrom$V1)
        cbb <- cut(distance, breaks,labels = FALSE)
        rcp <- tapply(m.chrom$V3, cbb, mean)
        rcp <- rcp[rcp != 0 & !is.na(rcp)]
        # SEM Calc
        robin_sem <- function(x) sqrt(var(x)/length(x))
        rcpsem <- tapply(m.chrom$V3, cbb, robin_sem)
        rcpsem <- rcpsem[names(rcpsem) %in% names(rcp)]
        rcpsem[is.na(rcpsem)] <- 0
        ###
        # divided by the number of non-zero fields of bin X
        tccb <- table(cbb)

        ourRCP <- rcp/unname(tccb[names(tccb) %in% names(rcp)])
        attributes(ourRCP) <- NULL

        ourSEM <- rcpsem/unname(tccb[names(tccb) %in% names(rcpsem)])
        attributes(ourSEM) <- NULL

        #m.sum <- sum(m.chrom$V3)
        dat <- dplyr::data_frame(breaks[as.numeric(names(rcp))],
                                 ourRCP,
                                 ourSEM)
        colnames(dat) <- c("distance", "prob", "SEM")
      } else {
        m.chrom <- dplyr::filter(m.chrom, V1 < V2)
        m.chrom$V3[is.na(m.chrom$V3)] <- 0
        distance = resolu*(m.chrom$V2-m.chrom$V1)
        cbb <- cut(distance, breaks,labels = FALSE)
        rcp <- tapply(m.chrom$V3, cbb, mean)
        rcp <- rcp[rcp != 0 & !is.na(rcp)]
        # SEM Calc
        robin_sem <- function(x) sqrt(var(x)/length(x))
        rcpsem <- tapply(m.chrom$V3, cbb, robin_sem)
        rcpsem <- rcpsem[names(rcpsem) %in% names(rcp)]
        rcpsem[is.na(rcpsem)] <- 0
        ###
        m.sum <- sum(m.chrom$V3)
        dat <- dplyr::data_frame(breaks[as.numeric(names(rcp))],
                                 rcp/m.sum,rcpsem/m.sum)
        colnames(dat) <- c("distance", "prob", "SEM")
      }
      if(nrow(dat) == 0){next()}
      #add names to data.frame
      if(standard){
        dat$sample <- paste("Exp.", i)
      }else{
        dat$sample <- experimentList[[i]]$NAME
      }

      dat$sample <- experimentList[[i]]$NAME
      dat$chrom <- chrom
      dat$color <- experimentList[[i]]$COL
      d <- rbind(d, dat)
    }
  }
  d$sample <- factor(d$sample)
  d$color <- as.character(d$color)
  experimentList[[i]]$ABS <- as.data.frame(experimentList[[i]]$ABS)
  return(d)
}
