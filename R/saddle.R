
#create an observed over expected matrix
obs.exp.matrix_saddle <- function( mat, correct=F, outlier.correct = 0.995, lowess = F ){
  pos <- which(mat > -1, arr.ind=T)
  val <- mat[mat > -1]
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
  return(obs.exp)
}


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

#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch_saddle <- function( x ){
  temp <- cumsum(c(1,diff(x) - 1))
  temp2 <- rle(temp)
  x[which(temp == with(temp2, values[which.max(lengths)]))]
}

saddleBinsPerChrom = function(exp, nBins, errEVcorrect = T, ChIP = NULL, chrom, start = NULL, end = NULL, cutOffQuantile =.995, closeCis = NULL, CS = NULL){
  # if start or end is NULL, use whole chrom
  if(is.null(start) | is.null(end)){
    start = 0
    end = max(exp$ABS[exp$ABS$V1 == chrom,3])
  }

  # get matrix
  inMat = select.subset(exp, chrom = chrom, start = start,  end =end)

  # bin ChIP
  peaks = NULL
  if(!is.null(ChIP)){
    ChIP[,1] = as.character(ChIP[,1])
    peaks = ChIP[ChIP[,1] == chrom & ChIP[,2] >= min(inMat$x) & ChIP[,3] <= max(inMat$x),]

    cuts = as.numeric(cut(rowMeans(peaks[,2:3]),breaks =  inMat$x, include.lowest = T))

    peaks = table( factor(cuts, levels = 1:max(ncol(inMat$z))))
  }

  # get bins with zero-sum
  xaxis = inMat$x
  inMat = as.matrix(inMat$z)
  ZSMAT = inMat[colSums(inMat) != 0,colSums(inMat) != 0]
  xaxis = xaxis[colSums(inMat) != 0]
  if(!is.null(ChIP)){peaks = peaks[colSums(inMat) != 0]}

  # get EV of ZSMAT
  OE = obs.exp.matrix_saddle(ZSMAT, correct = F, lowess = T)
  OE[!is.finite(OE)] <- 1
  Q = quantile(OE, c(1-cutOffQuantile, cutOffQuantile))
  OE[OE < Q[1]] = Q[1]
  OE[OE > Q[2]] = Q[2]

  # get compartment-score
  EV = NULL
  if(is.null(CS)){
    EV <- eigen(OE-1)
    EV = EV$vectors[,1]
  } else {
    #  xaxis: middel of eacht bin

    # use cut to get mean CS per hi-c bin.
    cuts = as.numeric(cut(rowMeans(CS[,2:3]),breaks = xaxis, include.lowest = T))

    CSs = tapply(CS[,4], cuts, mean, na.rm = T)
    tmp = which(!1:ncol(OE) %in%  rownames(CSs))
    tmp = c(1:ncol(OE))[tmp]
    CSs = c(CSs, setNames(rep(0, length(tmp)), tmp))
    CSs = CSs[order(names(CSs))]
    EV = CSs
  }


  if(errEVcorrect){
    tmp = length(EV[EV < 0]) /  length(EV)
    if(tmp < 0.1 | tmp > 0.9){
      #message(paste(exp$NAME,chrom,start))
      return(list(dat = NULL, EVbins = NULL))
    }
  }

  # flip EV on basis of chip-peaks
  if(!is.null(ChIP)){
    TEST = median(peaks[EV > 0]) > median(peaks[EV < 0])
    #message( TEST )
    if(TEST & !is.na(TEST)){
      tmp = EV
      EV = tmp*-1
    }
  } else {
    warning('No ChIP was given.\nIt is highly recommend that you\nadd some active regions such as H3K27ac!')
  }

  # set close-cis to NA
  if(!is.null(closeCis)){
    OE[is.na(OE)] = 1
    NAwidth = closeCis / exp$RES
    size = nrow(OE)
    idxDF = expand.grid(seq(1, size), seq(1, size))
    idxDF$D = abs(idxDF[,1] - idxDF[,2])
    idxDF = idxDF[idxDF$D <= NAwidth,]
    OE[cbind(idxDF[,1], idxDF[,2])] = NA
  }

  # find Qbins
  # 1 is A, 5 is B
  QBINS = quantile(EV*-1, seq(0,1,length.out = nBins+1))
  QCUTS = as.numeric(cut(EV*-1, QBINS, include.lowest = T) )

  EVbins = rev(tapply(EV, QCUTS, FUN = mean))

  df = data.frame()
  for(i in 1:nBins){
    for(j in 1:nBins){
      score = OE[QCUTS == i,QCUTS == j]
      score = mean(score, na.rm = T)
      df = dplyr::bind_rows(df, data.frame(i,j, score))
    }
  }

    return(list(dat = df, EVbins = EVbins))

}

#' Compute compartment-v-compartment scores
#'
#' Splits the range of compartment-scores in nBins bins (bin one having the highest compartment-scores and bin n the lowest).
#' Generates an average O/E contact-score for each compartment-score bin.
#' Sorts matrix on compartment-score.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param ChIP BED-dataframe containing active sites (e.g. H3K27ac-peaks).
#' @param chromsToUse Do the computation for a subset of chromosomes.
#' @param verbose Produces a progress-indication.
#' @param nBins The number of bins to split the compartment-score.
#' @param closeCis Do not take into account interactions close-by (10Mb)
#' @param CS Use a bedgraph-df of compartment-scores instead of creating it on the fly.The resolution should be the same as the Hi-C.
#' @note
#' A binsize of 5 will produce a plot similar to FLyamer et al. (2017), while 100 will produce a plot similar to Bonev et al. (2017). The increase in time with more nBins is in exponential time at this point in time. Future versions will tackle this.
#' @examples
#' # run saddleBins on all chromosomes with 25 bins.
#' saddle_WT = saddleBins(exp = Hap1_WT_1MB, ChIP = H3K27acPeaks,nBins = 25)
#'
#' # plot saddle-plot
#' visualise.saddle(list(saddle_WT), crossLines = T, addText = T)
#'
#' # plot compartment-strengths
#' CS = visualise.compartmentStrength(list(saddle_WT))
#' median(tmp$score)
#'
#' @return A log2(O/E) matrix and a DF of compartment-scores.
#' @import data.table
#' @export
saddleBins = function(exp, ChIP = NULL, chromsToUse = NULL, nBins =5, CS = NULL, closeCis = NULL, verbose = T){

  if(is.null(chromsToUse)){
    chromsToUse = exp$CHRS
  }

  df = data.frame()
  evDF = data.frame()
  for( C in chromsToUse){
    if(verbose){
      message( paste0(exp$NAME, ": ", C))
    }
    mat <- selectData_saddle( exp, C, C)

    cent <- which(apply(mat$z,1,sum)==0)
    cent <- largest.stretch_saddle(cent)
    chromStructure <- data.frame(chrom = C,
                                 start = 0,
                                 startC = min(cent)*exp$RES,
                                 endC = max(cent)*exp$RES,
                                 end= max(exp$ABS[exp$ABS$V1 == C,3]))



    if( chromStructure[1,3] < 1e6){

      # only Q
      Q = saddleBinsPerChrom(exp = exp,
                             ChIP = ChIP,
                             chrom = chromStructure[1,1],
                             start = chromStructure[1,4],
                             end = chromStructure[1,5],
                             nBins = nBins,
                             closeCis = closeCis, CS = CS)

      df_this = data.frame(sample = exp$NAME,
                           chrom = C,
                           arm = 'Q',
                           Q$dat)

      df = rbind(df, df_this)
      evDF = rbind(evDF,
                   data.frame(sample = exp$NAME,
                              chrom = C,
                              arm = 'Q',
                              bin = 1:length(Q$EVbins),
                              EV =  Q$EVbins))
    } else {
      # Q and P
      P = saddleBinsPerChrom(exp,
                             ChIP = ChIP,
                             chrom = chromStructure[1,1],
                             start = chromStructure[1,2],
                             end = chromStructure[1,3],
                             nBins = nBins,
                             closeCis = closeCis, CS = CS)

      df_this = data.frame(sample = exp$NAME,
                           chrom = C,
                           arm = 'P', P$dat)

      df = rbind(df, df_this)
      evDF = rbind(evDF, data.frame(sample = exp$NAME,
                                    chrom = C,
                                    arm = 'P',
                                    bin = 1:length(P$EVbins),
                                    EV = P$EVbins))

      Q = saddleBinsPerChrom(exp = exp,
                             ChIP = ChIP,
                             chrom = chromStructure[1,1],
                             start = chromStructure[1,4],
                             end = chromStructure[1,5],
                             nBins = nBins,
                             closeCis = closeCis, CS = CS)

      df_this = data.frame(sample = exp$NAME,
                           chrom = C,
                           arm = 'Q',
                           Q$dat)

      df = rbind(df, df_this)
      evDF = rbind(evDF, data.frame(sample = exp$NAME,
                                    chrom = C,
                                    arm = 'Q',
                                    bin = 1:length(Q$EVbins),
                                    EV = Q$EVbins))
    }

  }
  df$score = log2(df$score)
  #message(evDF)
  return(list(matrix = df, EV = evDF))

}

#' Plot compartment-strengths
#'
#' Takes the input from saddleBins and produces a boxplot of compartment-strenghts.
#'
#' @param SBoutList A list of outputs trom saddleBins.
#' @return A plot and an invisible underlying dataframe containing per-sample and -chromosome-arm the compartment-score.
#' @examples
#' # run saddleBins on all chromosomes with 25 bins.
#' saddle_WT = saddleBins(exp = Hap1_WT_1MB, ChIP = H3K27acPeaks,nBins = 25)
#'
#' # plot compartment-strengths
#' CS = visualise.compartmentStrength(list(saddle_WT))
#' median(tmp$score)
#'
#' @export
visualise.compartmentStrength = function(SBoutList){
  strengthDF = data.frame()
  for(i in 1:length(SBoutList)){
    dat = SBoutList[[i]]
    dat$matrix$CC = 'XX'
    MAXbin = max(dat$matrix$i)
    # how many binsare in 20%
    binsTOse = floor(MAXbin * .2)
    dat$matrix$unLog = 2 ** dat$matrix$score
    dat$matrix[dat$matrix$i <= binsTOse & dat$matrix$j <= binsTOse,"CC"] = "AA"
    dat$matrix[dat$matrix$i <= binsTOse & dat$matrix$j >= MAXbin-binsTOse+1,"CC"] = "AB"
    dat$matrix[dat$matrix$i >= MAXbin-binsTOse+1 & dat$matrix$j >= MAXbin-binsTOse+1,"CC"] = "BB"
    dat$matrix = dat$matrix[dat$matrix$CC != 'XX',]

    tmp = dplyr::summarise(dplyr::group_by(dat$matrix,
                                           sample,
                                           chrom,
                                           arm,
                                           CC),score = mean(unLog))

    for(S in unique(tmp$sample)){

      # get chroms in this sample
      meta = unique(tmp[tmp$sample == S,c("chrom","arm")] )

      for(C in unique(meta$chrom)){


        for(A in unique(unname(unlist(meta[meta$chrom == C,"arm"])))){
          tmpi = tmp[tmp$sample == S & tmp$chrom == C & tmp$arm == A, ]
          strength = log(tmpi[tmpi$CC == 'AA',"score"] * tmpi[tmpi$CC == 'BB',"score"] / tmpi[tmpi$CC == 'AB',"score"]**2)
          strengthDF = rbind(strengthDF, data.frame(S, C, A, strength))
        }
      }
    }

  }
  RANGE = range(strengthDF$score, na.rm = T)
  RANGE[1] = floor(RANGE[1]) *.95
  RANGE[2] = ceiling(RANGE[2])*1.05

  boxplot(split(strengthDF$score,
                strengthDF$S),
          ylim = RANGE,
          ylab = "compartment strength")
  invisible(strengthDF)
}

rotate <- function(x) t(apply(x, 2, rev))

#' Plot saddle-plots
#'
#' Takes the input from saddleBins and produces saddleplots.
#'
#' @param SBoutList A list of outputs trom saddleBins.
#' @param addText Add annotations in the corners for legibility.
#' @param crossLines Plot crosslines in matrix. Crossline denote border of- A and B-bins. If not set, will plot them when nBins  >= 10.
#' @param zlim Set the z-lims of the matrix.
#' @param EVlim If set, use a different set of z-lims for the compartment-scores.
#' @param square Set this to FALSE if plots are unaligned: will not produce a rigidly square matrix.
#' @return A saddle-plot
#' @examples
#' # run saddleBins on all chromosomes with 25 bins.
#' saddle_WT = saddleBins(exp = Hap1_WT_1MB, ChIP = H3K27acPeaks,nBins = 25)
#'
#' # plot saddle-plot
#' visualise.saddle(list(saddle_WT), crossLines = T, addText = T)
#' @export
visualise.saddle = function(SBoutList, addText = T, zlim = c(-0.5,0.5), EVlim = NULL, square = T, crossLines = NULL){
  par_temp = par()
  if(is.null(EVlim)){
    EVlim = zlim
  }

  df_mat = NULL
  df_ev = NULL
  NSAMPLE = length(SBoutList)
  SAMPLENAMES = c()
  for(i in 1:NSAMPLE){
    SAMPLENAMES = c(SAMPLENAMES, as.character(unique(SBoutList[[i]]$matrix$sample) ))
    df_mat = rbind(df_mat, SBoutList[[i]]$matrix)
    df_ev = rbind(df_ev, SBoutList[[i]]$EV)
  }

  if(is.null(crossLines)){
    if(max(df_mat$i) < 10){
      warning('Number of bins is lower that 10: not showing crossLines')
      crossLines = F
    } else {
      crossLines = T
    }
  }

  # make canvas
  lay = matrix(0, nrow = 4, ncol =  NSAMPLE)
  lay[1,] = 1:NSAMPLE
  for(i in 1:NSAMPLE){
    lay[2:4,i] = i + NSAMPLE
  }
  layout(lay)
  par(mar=rep(1,4), xaxs="i", yaxs="i")


  # plot ev-tracks
  EV0 = c()
  forEVPLot = aggregate(df_ev$EV, by = list(df_ev$sample, df_ev$bin), mean)
  forEVPLot[,3] = forEVPLot[,3] * -1
  for(S in SAMPLENAMES){
    tmp = setNames(forEVPLot[forEVPLot[,1] == S, 2:3], c("x","y"))


    tmp.bk = tmp
    crossPoint = approx(y = tmp.bk$x, x = tmp.bk$y, xout = 0)$y
    EV0 = c(EV0, crossPoint/max(forEVPLot[,2]))
    X <- tmp.bk$x
    y.low <- rep(0, length(tmp.bk$x))
    y.high <- tmp.bk$y

    plot(X,y.high,type = 'n', ylim = EVlim, axes = F, main = S)
    lines(X, y.low, col = 'black')
    lines(X, y.high, col = 'black')
    polygon(c(X, rev(X)), c(y.high, rev(y.low)),
            col = "black", border = NA); box()

    if(crossLines){abline(v = crossPoint, lty = 2)}

  }

  # plot matrices
  if(square){
    par(pty = 's')
  }
  PAL = colorRampPalette(c('#0571b0','white', "#ca0020"))
  for(Si in 1:NSAMPLE){
    S = SAMPLENAMES[Si]
    tmp = df_mat[df_mat$sample == S,]
    tmp = setNames(aggregate(tmp$score, by = list(tmp$i, tmp$j), mean), c("x",'y','z') )
    M = matrix(0, nrow = max(tmp$x), ncol = max(tmp$y))
    M[cbind(tmp$x,tmp$y)] = tmp$z

    M[M < zlim[1]] = zlim[1]
    M[M > zlim[2]] = zlim[2]

    image(rotate(M), zlim = zlim, col = PAL(11), axes = F); box()

    if(addText){
      text(x = .05, y = .95, label = 'AA')
      text(y = .05, x = .95, label = 'BB')
      text(y = .05, x = .05, label = 'AB')
      text(y = .95, x = .95, label = 'BA')
    }
    if(crossLines){
      abline(v = (EV0[Si])-(.5/max(tmp$y)) ,lty = 2)
      abline(h = 1-((EV0[Si])-(.5/max(tmp$y))) ,lty = 2)
    }

  }
  if(square){
    par(pty = 'm')
  }
  par(par_temp)

}
