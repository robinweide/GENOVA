#select a matrix of interactions for between two chromosomes
selectData_saddle <- function (exp, chrom1, chrom2){
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

#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
  temp <- cumsum(c(1,diff(x) - 1))
  temp2 <- rle(temp)
  x[which(temp == with(temp2, values[which.max(lengths)]))]
}




#' Get compartment-scores
#'
#' This uses the per-arm compartment-score algorithm of G. Fudenberg.
#' The procedure is an eigenvector decomposition of the observed/expected matrix minus one.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chromsToUse The chromosome(s) of interest. Omitting this leads to the computation of all chromosomes.
#' @param empericalCentromeres Emperically determine the centromeres (i.e. based on the largest white stripe in the Hi-C matrix).
#' @param comparableBed An optional dataframe of active regions (e.g. H3K27ac-peaks), which will be used to flip the A/B signal.
#' @param comparableTrack An optional earlier output of compartment.score, which will be used to flip the A/B signal. Or one can use GC-content to ensure that regions with positive scores are postiviely correlated with GC-content.
#' @param shuffle Will shuffle the Hi-C data if set to true. Useful for showing compartment-scores if no A/B structure.
#' @param pEV,qEV Sometimes, the second or third EV will better resemble the compartmentalisation. In this case, you can choose to get the 2nd or 3rd EV per chromosome-arm.
#' @return A BED-like dataframe with a fourth column, that holds the compartment-scores, and a fifth column, which holds the name of the experiment.The centromere will be skipped.
#' @examples
#' # Compute the compartment-scores of chromosome 20.
#' WT_40kb_CS = compartment.score(WT_40kb, chromsToUse = 'chr20')
#'
#' # Plot chromosome 20 with skipAnn = T to allow for an
#' # insertion op a plot on top.
#' hic.matrixplot(WT_40kb, chrom = 'chr20', start = 0,
#' end = 63025520, skipAnn = T, cut.off = 25)
#'
#' # Plot the compartment-score as a histogram
#' plot(WT_40kb_CS$start, WT_40kb_CS$compartmentScore,
#' t = 'h', axes = F); axis(2)
#' @export
compartment.score <- function(exp, chromsToUse = NULL, empericalCentromeres = T,
                              pEV = 1, qEV =1, comparableBed = NULL,
                              comparableTrack = NULL, shuffle = F, verbose = T){
  if(is.null(chromsToUse)){
    chromsToUse <- exp$CHRS
  }

  outList <- list()
  for(chrom in chromsToUse){
    if(verbose){message(chrom)}
    tmp <- compartment.score.chr(exp = exp, chrom = chrom,
                                 empericalCentromeres = empericalCentromeres,
                                 pEV = pEV, qEV =qEV,
                                 comparableBed = comparableBed,
                                 comparableTrack = comparableTrack,
                                 shuffle = shuffle)
    outList[[chrom]] <- tmp
  }
  outDF <- data.table::rbindlist(outList)

  return(outDF)
}

compartment.score.chr <- function(exp, chrom = "chr2", empericalCentromeres = T, pEV = 1, qEV =1, comparableBed = NULL, comparableTrack = NULL, shuffle = F){
  centChrom = NULL

  if(empericalCentromeres){
    mat <- selectData_saddle( exp, chrom, chrom)

    cent <- which(apply(mat$z,1,sum)==0)
    cent <- largest.stretch(cent)
    centChrom <- data.frame(chrom = chrom,
                            startC = min(cent)*exp$RES,
                            endC = max(cent)*exp$RES)
  } else if(!is.null(exp$CENTROMERES) & chrom %in% exp$CENTROMERES[,1]){
    centChrom = exp$CENTROMERES[exp$CENTROMERES[,1] == chrom,]
  } else {
    stop("User doesn't want empericalCentromeres,\nbut chromosome is not found in exp$CENTROMERES.\nQuitting...")
  }

  chromSize <- max(exp$ABS[exp$ABS$V1 == chrom,3])
  first_startLocVector <- NULL
  second_startLocVector <- NULL

  ###
  # first_arm
  ###
  C = chrom
  S = 0
  E = centChrom[centChrom[,1] == chrom,2]

  if(centChrom[1,2] > ( exp$RES*10)){ # otherwise too small to do
    # get all bins of this region
    first_startLocVector <-  exp$ABS[exp$ABS$V1 == C & exp$ABS$V2 <= E,3] - exp$RES

    ss <- select.subset(exp, C, S, E)
    if(shuffle){
      ss <- suppressMessages(shuffleHiC(ss))
    }
    ss_eigen <- get.eigen(ss)
    compScore <- ss_eigen$vectors[,pEV] * (ss_eigen$values[pEV]  ** 0.5)

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
        if(abs(CORTEST) > 0.25){
          compScore <- compScore * -1
        } else {
          warning("The correlation with the comparable track is lower than |0.25| on the p-arm of ", chrom, ".")
          compScore <- compScore * -1
        }
      }
    }
    firstArm <- compScore
  } else {
    firstArm <- NULL
  }


  ###
  # second_arm
  ###
  C = chrom
  S = centChrom[centChrom[,1] == chrom,3]
  E = max(exp$ABS[exp$ABS$V1 == C,3])

  if( (chromSize - centChrom[1,3]) > ( exp$RES*10)){ # otherwise too small to do
    # get all bins of this region
    second_startLocVector <-  exp$ABS[exp$ABS$V1 == C & exp$ABS$V2 >= S,3] - exp$RES
    ss <- select.subset(exp, C, S, E)
    if(shuffle){
      ss <- suppressMessages(shuffleHiC(ss))
    }
    ss_eigen <- get.eigen(ss)
    compScore <- ss_eigen$vectors[,qEV] * (ss_eigen$values[qEV]  ** 0.5)

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
      i <- findInterval(comparableBed[,2], second_startLocVector)
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
        if(abs(CORTEST) > 0.25){
          compScore <- compScore * -1
        } else {
          warning("The correlation with the comparable track is lower than |0.25| on the q-arm of ", chrom, ".")
          compScore <- compScore * -1
        }
      }
    }
    secondArm <- compScore
  } else {
    secondArm <- NULL
  }
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
