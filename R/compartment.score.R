obs.exp.matrix <- function( mat, correct=F, outlier.correct = 0.995, lowess = F ){
  mat$z[is.na(mat$z)] <- 0
  pos <- which(mat$z > -1, arr.ind=T)
  val <- mat$z[mat$z > -1]
  d <- abs(pos[,1]-pos[,2])
  if(correct){
    norm.factor <- tapply(val, d, function(x) mean( x[x < quantile(x, outlier.correct)] ) )
  }else{
    norm.factor <- tapply(val, d, mean)
  }
  if(lowess){
    x <- as.numeric(names(norm.factor)); y <- norm.factor
    norm.lowess <- lowess(x, y, f = 0.01, iter=20)
    norm.factor <- norm.lowess$y
    names(norm.factor) <- norm.lowess$x
  }

  val <- val/norm.factor[as.character(d)]
  obs.exp <- matrix(val, ncol = length(val)**0.5)
  list(x=mat$x, y=mat$y, z=obs.exp)
}


#input is an observed over expected matrix
get.eigen <- function( mat, with.eigen.value = T, outlier.correct = 0.995 ){
  oe <- obs.exp.matrix( mat )
  #set non-finite values to 1
  oe$z[!is.finite(oe$z)] <- 1
  #remove outliers
  th <- quantile(oe$z, outlier.correct)
  oe$z[oe$z > th] <- th

  ev <- eigen(oe$z - 1)
  ev
}

#' Get compartment-scores for a chromosome
#'
#' This uses the per-arm compartment-score algorithm of G. Fudenberg.
#' The procedure is an eigenvector decomposition of the observed/expected matrix.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom The chromosome of interest.
#' @param comparableTrack An optional vector, which will be used to flip the A/B signal.
#' For example, one can use GC-content to ensure that regions with positive scores are postiviely correlated with GC-content.
#' @param shuffle Will shuffle the Hi-C data if set to true. Useful for showing compartment-scores if no A/B structure.
#' @return A BED-like dataframe with a fourth column, that holds the compartment-scores, and a fifth column, which holds the name of the experiment.
#' @export
#'
compartment.score <- function(exp, chrom = "chr2", comparableBed = NULL, comparableTrack = NULL, shuffle = F){

  if(is.null(exp$CENTROMERES)){
    stop('Please store centromere-data in construct.experiment.\n')
  }
  chromSize <- max(exp$ABS[exp$ABS$V1 == chrom,3])

  ###
  # first_arm
  ###
  C = chrom
  S = 0
  E = exp$CENTROMERES[exp$CENTROMERES$V1 == chrom,2]

  # get all bins of this region
  first_startLocVector <-  exp$ABS[exp$ABS$V1 == C & exp$ABS$V2 <= E,3] - exp$RES

  ss <- select.subset(exp, C, S, E)
  if(shuffle){
    ss <- suppressMessages(shuffleHiC(ss))
  }
  ss_eigen <- get.eigen(ss)
  compScore <- ss_eigen$vectors[,1] * (ss_eigen$values[1]  ** 0.5)

  # If there is a comparison-vector, now is the chance to use it...

  ##### making a GC-comparison-vector
  # GCdf <- read.delim('~/Data/References/zebraFish/40kb.gc', h = F,stringsAsFactors = F) # to be made with GCContentByInterval from GATK
  # library(reshape2)
  # GGC <- colsplit(GCdf$V1,pattern = '[-:]', names = c('seqnames','start','end'))
  # GGC$perc <- GCdf$V2
  # comparableTrack <- GGC[GGC$seqnames ==C,4]
  ####

  if(!is.null(comparableBed)){

    comparableBed <- comparableBed[order(comparableBed[,1],comparableBed[,2]),]
    comparableBed <- comparableBed[comparableBed[,1]==C,]

    #overlap the peaks with the windows
    i <- findInterval(comparableBed[,2], first_startLocVector)
    comparableBed.window <- table(i+1)
    i.up <- which(compScore > 0)
    i.down <- which(compScore < 0)

    if(wilcox.test(comparableBed.window[as.character(i.up)], comparableBed.window[as.character(i.down)])$p.value < 1e-5){
      if(median(comparableBed.window[as.character(i.up)], na.rm=T) < median(comparableBed.window[as.character(i.down)], na.rm=T)){
        compScore <- -compScore
      }
    }

  }

  if(!is.null(comparableTrack)){
    SMPLCOR <- tail(comparableTrack,length(compScore))
    CORTEST <- cor(compScore,SMPLCOR, method = 'spearman')
    if(CORTEST < 0){
      compScore <- compScore * -1

    }
  }

  firstArm <- compScore

  ###
  # second_arm
  ###
  C = chrom
  S = exp$CENTROMERES[exp$CENTROMERES$V1 == chrom,2]
  E = max(exp$ABS[exp$ABS$V1 == C,3])

  # get all bins of this region
  second_startLocVector <-  exp$ABS[exp$ABS$V1 == C & exp$ABS$V2 >= S,3] - exp$RES

  ss <- select.subset(exp, C, S, E)
  if(shuffle){
    ss <- suppressMessages(shuffleHiC(ss))
  }
  ss_eigen <- get.eigen(ss)
  compScore <- ss_eigen$vectors[,1] * (ss_eigen$values[1]  ** 0.5)

  # If there is a comparison-vector, now is the chance to use it...

  ##### making a GC-comparison-vector
  # GCdf <- read.delim('~/Data/References/zebraFish/40kb.gc', h = F,stringsAsFactors = F) # to be made with GCContentByInterval from GATK
  # library(reshape2)
  # GGC <- colsplit(GCdf$V1,pattern = '[-:]', names = c('seqnames','start','end'))
  # GGC$perc <- GCdf$V2
  # comparableTrack <- GGC[GGC$seqnames ==C,4]
  ####

  if(!is.null(comparableBed)){

    comparableBed <- comparableBed[order(comparableBed[,1],comparableBed[,2]),]
    comparableBed <- comparableBed[comparableBed[,1]==C,]

    #overlap the peaks with the windows
    i <- findInterval(comparableBed[,2], first_startLocVector)
    comparableBed.window <- table(i+1)
    i.up <- which(compScore > 0)
    i.down <- which(compScore < 0)

    if(wilcox.test(comparableBed.window[as.character(i.up)], comparableBed.window[as.character(i.down)])$p.value < 1e-5){
      if(median(comparableBed.window[as.character(i.up)], na.rm=T) < median(comparableBed.window[as.character(i.down)], na.rm=T)){
        compScore <- -compScore
      }
    }

  }

  if(!is.null(comparableTrack)){
    SMPLCOR <- tail(comparableTrack,length(compScore))
    CORTEST <- cor(compScore,SMPLCOR, method = 'spearman')
    if(CORTEST < 0){
      compScore <- compScore * -1

    }
  }
  secondArm <- compScore

  ###
  # return
  ###
  outDF <- data.frame(seqnames = C,
             start = c(first_startLocVector,second_startLocVector),
             end = c(first_startLocVector,second_startLocVector)+exp$RES,
             compartmentScore = c(firstArm,secondArm),
             name = exp$NAME)

  return(outDF)
}
