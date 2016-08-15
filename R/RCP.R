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
RCP <- function(experimentList, chromsToUse, maxDistance = 1e09, verbose = F){
  amountOfSamples <- length(experimentList)

  d <- data_frame(distance = integer(),
                  prob = numeric(),
                  sample = integer())

  for(Ci in 1:length(chromsToUse)){
    chrom <- chromsToUse[Ci]
    for(i in 1:amountOfSamples){
      if(verbose){cat(paste0('Chromosome ', Ci,' of ', length(chromsToUse), ' chromosomes\n' ,  'Sample ', i,' of ', amountOfSamples, ' samples'), "\r")}
      experiment <- experimentList[[i]]
      BED_ABS <- experiment$ABS
      m <- experiment$ICE
      resolu <- experiment$RES

      BED_ABS <- data.table::data.table(BED_ABS)
      data.table::setkey(BED_ABS, V1)
      idx1 <-BED_ABS[list(chrom)]

      x <- rep(idx1$V4, length(idx1$V4))
      y <- rep(idx1$V4, each=length(idx1$V4))
      xydf.a <- dplyr::data_frame(x,y)
      xydf.b <- dplyr::filter(xydf.a,abs(x-y) <= (maxDistance/resolu))

      x <- as.numeric(xydf$x)
      y <- as.numeric(xydf$y)

      breaks <- 10**seq(4,8,length.out=81)

      m.chrom <- m[list(x,y)]
      m.chrom <- dplyr::filter(m.chrom, V1 < V2)
      m.chrom$V3[is.na(m.chrom$V3)] <- 0
      distance = resolu*(m.chrom$V2-m.chrom$V1)
      cbb <- cut(distance, breaks,labels = FALSE)
      rcp <- tapply(m.chrom$V3,cbb,mean)
      rcp <- rcp[rcp!=0 & !is.na(rcp)]
      m.sum <- sum(m.chrom$V3)
      dat <- dplyr::data_frame(breaks[as.numeric(names(rcp))],rcp/m.sum)
      colnames(dat) <- c('distance',  'prob')
      dat$sample <- experiment$NAME
      dat$chrom <- chrom
      dat$color <- experiment$COL
      d <- rbind(d, dat)
    }
  }
  dat$sample <- factor(dat$sample)
  dat$color <- as.character(dat$color)
  return(d)
}
