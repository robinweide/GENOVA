#' Relative Contact Probabilities
#'
#' Produce a dataframe with the probabilities of contacts in a set of distance-bins.
#' Bins are created on a log scale, which leads to equal amounts of datapoints per bin.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param experiment List of experiment-objects from `construct.experiment()`.
#' @param chromsToUse A vector containing the chromosome-names of interest. If bedList is also given, only BED-entries in chromsToUse will be used.
#' @param bedList A named list of BED-like dataframse of the regions of interest. RCP will intersect this with the Hi-C bin by intersecting the middle of both ranges.
#' @param colors Override colors of experiments. Will be used for BEDs if only one experiment is given.
#' @param maxDistance The maximal distance to calculate RCPs for.
#' @param outlierCutoff Remove the [outlierCutoff] top and bottom percentages.
#' @param ignoreLengthWarning Force RCP to use maxDistance instead of longest chromosome if that is smaller than maxDistance.
#' @param verbose Produces a progress-indication.
#' @return A data.frame, containing:
#' @return \item{distance}{the distance-bin}
#' @return \item{prob}{the relative contact probability}
#' @return \item{SEM}{the SEM of the RCP}
#' @return \item{BED}{if and which BED is used}
#' @return \item{color}{the corresponding color of the lines, taken from exp$COL as default}
#' @return \item{sample}{the corresponding sample-name}
#' @return \item{chrom}{the corresponding}
#' @examples
#' # Calculate the RCP of chromosome 1
#' RCP_out = RCP(experimentList = list('WT' = WT_1MB),
#' chromsToUse = 'chr1')
#'
#' # Plot the RCP
#' visualise.RCP.ggplot(RCP_out)
#' @export
RCP <- function(experimentList, chromsToUse = NULL,  bedList = NULL, colors = NULL, maxDistance = NULL, ignoreLengthWarning = F, outlierCutoff = 1,verbose = F){

  BEDCOL = NULL
  if(is.null(colors)){
    colors = c()
    for(i in 1:length(experimentList)){
      colors = c(colors, unique(experimentList[[i]]$COL))
    }
    if(length(experimentList) == 1 & length(bedList) > 0){
      colors = rainbow(length(bedList))
    }
  } else {
    if(length(experimentList) == 1){
      BEDCOL = colors
    } else {
      colors = colors
    }
  }

  # check if all have the same resolution
  RESSES = c()
  for(i in 1:length(experimentList)){
    RESSES = c(RESSES , experimentList[[i]]$RES)
  }
  if( length(unique(RESSES)) != 1 ){warning("The resolutions of the different experiments are not the same")}

  amountOfSamples <- length(experimentList)
  exp.names <- c()
  if(is.null(chromsToUse)){
    chromsToUse <-  unique(experimentList[[1]]$ABS[,1])
  }

  # Compare largest chromosome to maxDistance
  largestChrom = max(experimentList[[1]]$ABS[experimentList[[1]]$ABS$V1 %in% chromsToUse,3])
  if(is.null(maxDistance)){
    maxDistance = largestChrom
  }
  if(maxDistance > largestChrom & !ignoreLengthWarning){
    warning(paste0("MaxDistance is -unnecessarily- too large.\nLargest chromosome is ", largestChrom, " long.\nSetting MaxDistance to ",largestChrom, "."))
    maxDistance = largestChrom
  }


  # Check for N chromosomes: more than 50: please specify chroms!
  if(length(chromsToUse) > 50){stop("Please restrict the amount of chromosomes with chromsToUse")}


  #check whether experiment names have been declared uniquely
  #otherwise use standard names for RCP
  exp.names = c()
  for( i in 1:length(experimentList)){
    exp.names <- c(exp.names, experimentList[[i]]$NAME)
  }
  standard <- FALSE
  if(length(exp.names) != length(unique(exp.names))){
    warning("Experiment names have not been declared uniquely, using standard names")
    standard <- TRUE
    exp.names = paste("Exp.", 1:length(experimentList))
  }

  d <- dplyr::data_frame(distance = integer(),
                  prob = numeric(),
                  sample = integer())
  dat <- NULL
  for(Ci in 1:length(chromsToUse)){
    chrom <- unique(chromsToUse[Ci])
    # check if we have a bed and, if yes, entries on this chomosome
    if( !is.null(bedList) ) { # bed present
      bedCHROMS = as.character(unique(unlist(lapply(bedList, function(x) unique(x[,1])) ) ))
      if( ! chrom %in% bedCHROMS   ){ # chrom not in bed

        next()
      }
    }
    for(i in 1:amountOfSamples){
      if(verbose){cat(paste0('Chromosome ', Ci,' of ', length(chromsToUse), ' chromosomes\r'))}
      dat <- NULL
      resolu <- experimentList[[i]]$RES
      BED_ABS <- data.table::data.table(experimentList[[i]]$ABS)
      data.table::setkey(BED_ABS, V1)
      idx1 <- BED_ABS[V1 == list(chrom)]

      chromLength = max(BED_ABS[chrom,3])
      centromerixIDXes = c()
      if(chrom %in% experimentList[[i]]$CENTROMERES[,1]){
        centroBed <- experimentList[[i]]$CENTROMERES[experimentList[[i]]$CENTROMERES[,1] == chrom,]
        centromerixIDXes =  unname(unlist(idx1[V2 >= centroBed[,2] & V3 <= centroBed[,3],4]))
        idx1 <- idx1[!V4 %in% centromerixIDXes]

      } else {
        # find biggest chrom to compare with
        allChroms <- table(data.table::data.table(experimentList[[i]]$ABS)$V1)
        if(length(allChroms) > 1){ # if there are more than 1 chromosome in the data, remove centromeres

          sizeChrom <- allChroms[names(allChroms)== chrom]
          allChroms <- allChroms[!names(allChroms) == chrom] # not the same chrom!
          toCompare <- names(which(allChroms == max(allChroms)))
          siztoCompare <- allChroms[names(allChroms)== toCompare]

          chromTheBiggest <- siztoCompare < sizeChrom
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

          idx1 <- idx1[!V4 %in% centromerixIDXes]
        }
      }

      centStart = min(centromerixIDXes)
      centEnd = max(centromerixIDXes)

      x <- rep(idx1$V4, length(idx1$V4))
      y <- rep(idx1$V4, each=length(idx1$V4))
      xydf <- dplyr::data_frame(x,y)
      xydf.b <- dplyr::filter(xydf,abs(x-y) <= (maxDistance/resolu))

      x <- as.numeric(xydf$x)
      y <- as.numeric(xydf$y)

      breaks <- 10**seq(4,log(maxDistance, base = 10),length.out=81)

      m.chrom <- experimentList[[i]]$ICE[list(x,y)]

      m.chrom_UP = m.chrom[V1 < centStart & V2 < centStart, ]
      m.chrom_DOWN = m.chrom[V1 > centEnd & V2 > centEnd, ]

      dat = data.frame()
      message(BEDCOL)
      if(!is.null(bedList)){
        if(!is.list(bedList)){stop('bedList must be a named list of dataframes!')}
        ABSchrom = experimentList[[i]]$ABS[experimentList[[i]]$ABS[,1] == chrom,]

        for(BLidx in 1:length(bedList)){
          bed = bedList[[BLidx]]
          BEDchrom = bed[bed[,1] == chrom,1:3]

          BEDchrom = BEDchrom[ order(BEDchrom[,2], BEDchrom[,3]), ]
          BEDchrom[ BEDchrom[,2] > BEDchrom[,3] , ] = BEDchrom[ BEDchrom[,2] > BEDchrom[,3] , c(1,3,2)]

          # which bins overlap with bed-entries?
          ii <- findInterval( apply(BEDchrom[,2:3], 1 ,mean) , apply(ABSchrom[,2:3], 1 ,mean) ) +1

          found = ABSchrom[unique(ii),]

          # plot(PCA[PCA$V1 == 'chr1', 2], PCA[PCA$V1 == 'chr1', 4], t = 'h')
          # points(found$V2, rep(60, nrow(found)), col = 'red', pch = 20)
          # points(BEDchrom$V2, rep(50, nrow(BEDchrom)), col = 'blue', pch = 20)

          bedi = found[,4]
          if(length(bedi) > 0) {
            tmp = coreRCP(rbind(m.chrom_DOWN,m.chrom_UP), resolu, breaks, outlierCutoff, bedi = bedi)
            if(nrow(tmp) > 0){
              tmp$BED = names(bedList)[BLidx]

              # if multiple samples, do not do this!
              if(length(experimentList) > 1){
                tmp$color <- colors[i]
              }else {
                tmp$color <- colors[BLidx]
              }

              dat = rbind(dat, tmp)
            }
          }

        }


      } else {
        dat = coreRCP(rbind(m.chrom_DOWN,m.chrom_UP), resolu, breaks, outlierCutoff)
        if(nrow(dat) > 0){
          dat$BED = NA
        }
        dat$color <- colors[i]
      }


      if(nrow(dat) == 0){next()}

      #add names to data.frame
      if(standard){
        dat$sample <- paste("Exp.", i)
      }else{
        dat$sample <- experimentList[[i]]$NAME
      }

      dat$chrom <- unique(chrom)


      d <- rbind(d, dat)

    }
  }


  d$sample <- factor(d$sample, levels = exp.names)
  #d$color <- as.character(d$color)
  if(!is.na(unique(d$BED)[1])){
    d$BED = factor(d$BED, levels = names(bedList))
  }
  return(d)
}


#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
  temp <- cumsum(c(1,diff(x) - 1))
  temp2 <- rle(temp)
  x[which(temp == with(temp2, values[which.max(lengths)]))]
}

outlierRemovedMean = function(x, Q = 0.005){ # Q is percentage
  x[is.na(x)] = 0
  boundaries = quantile(x, c(Q/100, 1-(Q/100)))
  x[x < boundaries[1]] = boundaries[1]
  x[x > boundaries[2]] = boundaries[2]
  return(mean(x))
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

robin_sem <- function(x) sqrt(var(x)/length(x))

# centroBed <- experimentList[[i]]$CENTROMERES[experimentList[[i]]$CENTROMERES[,1] == chrom,]
# centromerixIDXes =  unname(unlist(idx1[V2 >= centroBed[,2] & V3 <= centroBed[,3],4]))
# idx1 <- idx1[!V4 %in% centromerixIDXes]

# if bed is given, use that. otherwise YOLOOO
coreRCP = function(m.chrom, resolu, breaks, outlierCutoff, bedi = NULL){
  m.chrom <- dplyr::filter(m.chrom, V1 < V2)
  m.chrom$V3[is.na(m.chrom$V3)] <- 0
  if(!is.null(bedi)){

    # filter m.chrom
    m.chrom <- dplyr::filter(m.chrom, V1 %in% bedi)
  }


  distance = resolu*(m.chrom$V2-m.chrom$V1)
  cbb <- cut(distance, breaks,labels = FALSE)
  rcp <- tapply(m.chrom$V3, cbb, FUN = function(x){outlierRemovedMean(x, outlierCutoff )})
  rcp <- rcp[rcp != 0 & !is.na(rcp)]
  # SEM Calc
  rcpsem <- tapply(m.chrom$V3, cbb, robin_sem)
  rcpsem <- rcpsem[names(rcpsem) %in% names(rcp)]
  rcpsem[is.na(rcpsem)] <- 0
  ###
  m.sum <- sum(m.chrom$V3)
  dat <- dplyr::data_frame(breaks[as.numeric(names(rcp))],
                           rcp/m.sum,rcpsem/m.sum)
  colnames(dat) <- c("distance", "prob", "SEM")
  return(dat)
}
