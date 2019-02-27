getPixelMean <- function(MAT, pixWidth = 3, npix, NAasZero = T){
  low  <- ((npix+1 ) /2 ) - ((pixWidth-1 ) /2 )
  high <- ((npix+1 ) /2 ) + ((pixWidth-1 ) /2 )
  subMat <- MAT[low:high, low:high]
  if(NAasZero){
    subMat[is.na(subMat)] <- 0
  }
  return(mean(subMat))
}

getDonutMean <- function(MAT, pixWidth = 3, npix, NAasZero = T){
#  if(!is.null(pixWidth)){warning('Pixwidth is decrepated.')}

  # cut to smaller mat
  pixelLoc   <- ((npix-1 )/2)+1
  cutPix <- ((pixWidth +2) * 2)/2
  MAT <- MAT[(pixelLoc-cutPix):(pixelLoc+cutPix),(pixelLoc-cutPix):(pixelLoc+cutPix)]
  MAT[is.na(MAT)] <- 0

  pixelLoc <- ((ncol(MAT) -1 )/ 2)+1
  SCORE_pixel <- MAT[pixelLoc,pixelLoc]

  # remove middle 5x5
  extendedPWwidth = ((pixWidth-1)/2)*2

  MAT[(pixelLoc-extendedPWwidth):(pixelLoc+extendedPWwidth),
      (pixelLoc-extendedPWwidth):(pixelLoc+extendedPWwidth)] <- NA

  # score H
  SCORE_h <- MAT[pixelLoc, ]
  SCORE_h <- mean(SCORE_h[!is.na(SCORE_h)])

  # score V
  SCORE_v <- MAT[,pixelLoc]
  SCORE_v <- mean(SCORE_v[!is.na(SCORE_v)])

  # remove H and V
  MAT[, pixelLoc] <- NA
  MAT[pixelLoc, ] <- NA

  # score quantile
  SCORE_q <- MAT
  SCORE_q <- mean(SCORE_q[!is.na(SCORE_q)])

  # score daig. quantile
  MAT[1:pixelLoc,] <- NA
  SCORE_qd <- MAT[,1:pixelLoc]
  SCORE_qd <- mean(SCORE_qd[!is.na(SCORE_qd)])

  SCORES <- c(SCORE_pixel, SCORE_qd, SCORE_h, SCORE_v, SCORE_q)

  BOOLOUTCOME <- ((SCORES[1] - SCORES[-1])/SCORES[-1]) > .5
  MEDIANOUTCOME <- SCORES[1] / median(SCORES[-1])


  return(c(mean(BOOLOUTCOME), MEDIANOUTCOME))
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
#' @param enrichment Calculate log2-enrichment instead of contacts?
#' @param enrichmentType meanScore gives the pixel/mean(backgroundRegions) score.
#' medianBool gives the fraction of background-regions less
#' than 50 percent of the pixel (1 = a true loop).
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
quantifyAPA <- function(APAlist, enrichment = F, pixWidth = 3, speudoCount = 1,
                        enrichmentType = 'medianScore'){
  if((pixWidth %% 2) == 0){
    stop('Please use an uneven number for pixWidth')
  }
  if(!enrichmentType %in% c('meanBool', 'medianScore')){
    stop('enrichmentType must be either meanBool or medianScore.')
  }
  outDF <- list() # make a df with a line per loop (add column for color and name)

  resOut <- NULL
  for(i in 1:length(APAlist)){
    apaout <- APAlist[[i]]$rawMatList
    npix   <- nrow(apaout[[1]])
    avgsPix   <- lapply(apaout, getPixelMean, pixWidth = pixWidth, npix = npix)
    avgsPix   <- unlist(avgsPix)

    avgs <- NULL
    if(enrichment){
      tmp = unlist(lapply(apaout, min))
      if(is.null(speudoCount)){
        speudoCount = min(tmp[tmp != 0])
      }
      donut_out <- lapply(apaout, getDonutMean, npix = npix, pixWidth = pixWidth)
      donut_out <- unlist(donut_out)
      donut_out <- matrix(unlist(donut_out), ncol = 2, byrow = T)
      if(enrichmentType == 'meanBool'){
        avgs <- donut_out[,1]
      } else {
        avgs <- donut_out[,2]
      }

    } else {
      avgs <- avgsPix
    }


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



