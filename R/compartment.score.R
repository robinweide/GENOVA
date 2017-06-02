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
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom The chromosome of interest.
#' @param comparableTrack An optional vector, which will be used to flip the A/B signal.
#' For example, one can use GC-content to ensure that regions with positive scores are postiviely correlated with GC-content.
#' @return A vector with the compartment-score for every bin in the Hi-C matrix.
#' @export
#'
compartment.score <- function(exp, chrom = "chr2", comparableTrack = NULL){

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

  ss <- select.subset(exp$ICE, C, S, E, exp$ABS)
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

  if(!is.null(comparableTrack)){
    SMPLCOR <- head(comparableTrack,length(compScore))
    if(cor(compScore,SMPLCOR) < 0){
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

  ss <- select.subset(exp$ICE, C, S, E, exp$ABS)
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

  if(!is.null(comparableTrack)){
    SMPLCOR <- tail(comparableTrack,length(compScore))
    if(cor(compScore,SMPLCOR) < 0){
      compScore <- compScore * -1
    }
  }
  secondArm <- compScore

  ###
  # return
  ###
  # merged vector of P and Q, incrusing zeroes at centromere
  CentSize <- diff(unlist(exp$CENTROMERES[exp$CENTROMERES$V1 == chrom,2:3]))
  centAB <- rep(0,ceiling(CentSize/exp$RES))
  ABtrack <- c(firstArm,centAB,secondArm)

  return(ABtrack)
}
