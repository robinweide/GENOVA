#' Get a per-chromosome scaling plot
#'
#' Produce a dataframe with the probabilities of contacts in a set of distance-bins.
#' Bins are created on a log scale, which leads to equal amounts of datapoints per bin.
#'
#' @param experiment List of experiment-objects from `construct.experiment()`.
#' @param chromsToUse A vector containing the chromosome-names of interest.
#' @param maxDistance The maximal distance to calculate scalings for.
#' @param verbose Produces a progress-indication.
#' @return A data_frame with distance-bin and probabilities.
#' @export
RCP.region <- function(experimentList, bed , verbose = F){
  maxDistance = 5e06
  amountOfSamples <- length(experimentList)
  exp.names <- c()
  chromsToUse <-  unique(experimentList[[1]]$ABS[,1])
  chromsToUse <- chromsToUse[-1]
  chromsToUse <- chromsToUse[-24]
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

  for(Ci in 1:length(chromsToUse)){
    chrom <- chromsToUse[Ci]
    for(i in 1:amountOfSamples){
      if(verbose){cat(paste0('Chromosome ', Ci,' of ', length(chromsToUse), ' chromosomes\n' ,  'Sample ', i,' of ', amountOfSamples, ' samples'), "\r")}
      #experiment <- experimentList[[i]]
      #m <- experimentList[[i]]$ICE
      resolu <- experimentList[[i]]$RES


      BED_ABS <- data.table::data.table(experimentList[[i]]$ABS)
      data.table::setkey(BED_ABS, V1)
      idx1 <-BED_ABS[list(chrom)]
      # if(!is.null(idxToUse)){
      #   idx1 <- idx1[idx1$V4 %in% idxToUse,]
      # }
      x <- rep(idx1$V4, length(idx1$V4))
      y <- rep(idx1$V4, each=length(idx1$V4))
      xydf <- dplyr::data_frame(x,y)
      xydf.b <- dplyr::filter(xydf,abs(x-y) <= (maxDistance/resolu))

      x <- as.numeric(xydf$x)
      y <- as.numeric(xydf$y)

      breaks <- 10**seq(4,8,length.out=81)

      m.chrom <- experimentList[[i]]$ICE[list(x,y)]

      ### Does this ^^^ include 0's?
      # head(m.chrom)
      #   V1     V2 V3
      #   1: 278172 278172 NA
      #   2: 278173 278172 NA
      #   3: 278174 278172 NA
      #   4: 278175 278172 NA
      #   5: 278176 278172 NA
      #   6: 278177 278172 NA
      # head(m.chrom[! is.na(m.chrom$V3)])
      #   V1     V2       V3
      #   1: 279113 279113 256.3800
      #   2: 279113 279114 367.5193
      #   3: 279114 279114 288.4004
      #   4: 279113 279115 161.1691
      #   5: 279114 279115 245.6197
      #   6: 279115 279115 209.5638
      m.chrom <- dplyr::filter(m.chrom, V1 < V2)
      if(nrow(m.chrom) == 0){ next }
      m.chrom$V3[is.na(m.chrom$V3)] <- 0
      # head(m.chrom)
      #   V1     V2 V3
      #   1: 278172 278172 0
      #   2: 278173 278172 0
      #   3: 278174 278172 0
      #   4: 278175 278172 0
      #   5: 278176 278172 0
      #   6: 278177 278172 0

      # Get IDX's from bed-file
      sel <- NULL
      for(h in 1:nrow(bed)){
        sel <- c(sel,experimentList[[i]]$ABS[experimentList[[i]]$ABS[,1]==bed[h,1] & experimentList[[i]]$ABS[,2] >= bed[h,2] & experimentList[[i]]$ABS[,2] <= bed[h,3],4])
      }

      m.chrom <- dplyr::filter(m.chrom, V1 %in% sel, V2 %in% sel:sel+(maxDistance/resolu))

      m.chrom <- dplyr::filter(m.chrom, V1 < V2)

      distance = resolu*(m.chrom$V2-m.chrom$V1)
      cbb <- cut(distance, breaks,labels = FALSE)
      rcp <- tapply(m.chrom$V3,cbb,mean)
      rcp <- rcp[rcp!=0 & !is.na(rcp)]
      #m.sum <- sum(m.chrom$V3)
      dat <- dplyr::data_frame(breaks[as.numeric(names(rcp))],rcp)
      colnames(dat) <- c('distance',  'prob')
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
  dat$sample <- factor(dat$sample)
  dat$color <- as.character(dat$color)
  experimentList[[i]]$ABS <- as.data.frame(experimentList[[i]]$ABS)
  return(d)
}
