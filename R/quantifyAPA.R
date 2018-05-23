getPixelMean <- function(MAT, pixWidth = 3, npix, NAasZero = T){
  low  <- ((npix-1 ) /2 ) - ((pixWidth-1 ) /2 )
  high <- ((npix-1 ) /2 ) + ((pixWidth-1 ) /2 )
  subMat <- MAT[low:high, low:high]
  if(NAasZero){
    subMat[is.na(subMat)] <- 0
  }
  return(mean(subMat))
}


#' Get statistics from the centers of APA-results
#'
#' Takes the input from APA() and produces a list of average centers and
#' comparitory statistics, where centers is a square region of
#' pixWidth x pixWidth in the middle.
#'
#' @param APAlist A named list of outputs from APA().
#' @param pixWidth The width of the square to use. 1 will give you just the
#' center-point
#' @return Alist of two dataframes: data (per-loop and -sample average scores)
#' and stats (all-vs-all wilcoxon.test p-values and log2-foldchanges.)
#' @examples
#' # run APA() on sample(s)
#' APA_WT <- APA(experiment = WT, loop.bed = loopsWT)
#' APA_WA <- APA(experiment = WA, loop.bed = loopsWT)
#'
#' # run quantifyAPA()
#' APAstats <- quantifyAPA(list("WT" = APA_WT, "MED12" = APA_MD), pixWidth = 3)
#'
#' # get the statistics
#' APAstats$stats
#'
#' # plot boxplots
#' ggplot(APAstats$data, aes(x = sample, y = value)) +
#'   geom_boxplot()
#' @export
quantifyAPA <- function(APAlist, pixWidth = 3){

  outDF <- list() # make a df with a line per loop (add column for color and name)

  for(i in 1:length(APAlist)){
    apaout <- APAlist[[i]]$rawMatList
    npix   <- nrow(apaout[[1]])
    avgs   <- lapply(apaout, getPixelMean, pixWidth = pixWidth, npix = npix)
    avgs   <- unlist(avgs)
    df     <- data.frame(sample = names(APAlist)[i],
                         loopIDX = 1:length(apaout),
                         value = avgs)
    outDF[[i]] <- df
  }

  # stats
  statsDF <- NULL
  for(i in 1:length(outDF)){
    for(j in 1:length(outDF)){
      if(i == j){next()}
      if(!any(statsDF$sampleA == levels(outDF[[j]]$sample) & statsDF$sampleB == levels(outDF[[i]]$sample))){
        tst <- wilcox.test(outDF[[i]]$value, outDF[[j]]$value)$p.value
        df <- data.frame(sampleA = unique(outDF[[i]]$sample),
                         sampleB = unique(outDF[[j]]$sample),
                         Pwilcox = tst,
                         log2_BoverA  = log2(mean(outDF[[j]]$value) /
                                               mean(outDF[[i]]$value) ))
        statsDF <- rbind(statsDF, df)
      }
    }
  }

  outDF <- data.table::rbindlist(outDF)

  return(list(data = outDF, stats = statsDF))
}



