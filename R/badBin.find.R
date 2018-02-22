#' Find all bad bins
#'
#' Builds a representative RCP and tests per Hi-C bin, whether its log-linear fitted data correlates.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param p.treshold Cutoff of spearmans correlation.
#' @param getBed Produce a bed-like object.
#' @param plotBad Plot some of the bad bins for visual inspection.
#' @return A list containing a vector with the bad bin-id's (`badBins`). Optionally, it can output a bed-file (`bed`).
#' @examples
#' # Identify and plot the bad bins
#' bb_out <- badBin.find(experiment = WT_100kb, p.treshold = 0.1, getBed = T, plotBad = T)
#'
#' # plot the region from within hic.matrixplot
#' hic.matrixplot(exp1 = WT_100kb, chrom = 'chr7',  start = 26.75e6,  end=28.5e6, chip = list(bb_out$bed))
badBin.find <- function(experiment, p.treshold = 0.1, plotBad = F, getBed = T){
  ## INIT-phase
  cat("Initiation...\n")
  allbins <- unique(experiment$ICE$V1)
  maxbin <- max(allbins)
  oneToMaxBin <- 1:maxbin

  res <- experiment$RES
  ### Sample 1000 bins (more, if runtime is low enough)
  sampledBINs <- sample(unique(experiment$ICE$V1), size = 1000)

  cat("Generating base RCP... 0 percent\t\t\t\t\r\r")
  ### Make mean RCP of sampled bins -> RCP_base
  RCP_base <- rep(0, (1e6 / res)/2)
  for(i in 1:length(sampledBINs)){
    if(i%% 10 == 0){
      cat("Generating base RCP...", round((i/1000)*100, digits = 2), "percent\t\t\t\t\r")
    }
    bin <- sampledBINs[i]
    binStop <- (1e6 / res)/2 + bin
    x <- rep(bin, (1e6 / res)/2)
    y <- seq(bin,binStop-1)
    vals <- experiment$ICE[list(x,y)]
    vals[is.na(vals$V3)]$V3   <- 0
    ### Make sum RCP of sampled bins -> RCP_base
    RCP_base <- RCP_base + vals$V3
  }
  RCP_base <- RCP_base / 1000
  RCP_base.model <- lm(log(c(1:((1e6 / res)/2)))~ RCP_base)
  #
  cat("\n")
  cat("Iterating over matrix... 0 percent\r")
  COR.list <- list()
  for(I in oneToMaxBin){
    if(I%% 1000 == 0){
      cat("Iterating over matrix...", round((I/maxbin)*100, digits = 2), "percent\t\t\r")
    }
    bin <- I
    binStop <- (1e6 / res)/2 + bin
    x <- rep(bin, (1e6 / res)/2)
    y <- seq(bin,binStop-1)
    vals <- experiment$ICE[list(x,y)] # maybe also the other way around?
    if( sum(is.na(vals$V3)) == ((1e6 / res)/2)){# saves time: no need to do cor-test on NA's
      COR.list[[I]] <- c(bin, 0, 1)
    } else{
      for (i in seq_along(vals)) data.table::set(vals, i=which(is.na(vals[[i]])), j=i, value=0)
      exponential.model <- lm(log(c(1:((1e6 / res)/2)))~ vals$V3)
      b <- cor.test(fitted.values(exponential.model), fitted.values(RCP_base.model), method = 'spearman')
      COR.list[[I]] <- c(bin, b$estimate, b$p.value)
    }
  }
  COR_df <- data.frame(abs(matrix(unlist(COR.list), ncol = 3, byrow = T)))
  COR_df.ordered <- COR_df[order(COR_df$X3),]
  badBins <- COR_df.ordered[COR_df.ordered$X3 > p.treshold,1]
  # This doesn't store it...
  #experiment$MASK <- c(experiment$MASK,badBins)
  if(plotBad){
    par(mfrow = c(4,5))
    toPlot <- sample(badBins, size = 20)
    for(i in 1:20){
      badBin.plot(hicMat  = experiment$ICE, res = res, RCP_base = RCP_base, bin = toPlot[i])
    }
  }
  if(getBed){
    bed <- data.frame(seqnames = character(),
                      start = numeric(),
                      stop = numeric())
    for(j in badBins){
      bed <- rbind(bed, experiment$ABS[j,c(1,2,3)])
    }
    return(list(badBins = badBins, bed = bed))
  } else {
    return(list(badBins = badBins))
  }

}
