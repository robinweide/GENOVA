# SADDLE

#===============================================================================
# helper-functions
#===============================================================================
#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
  temp <- cumsum(c(1,diff(x) - 1))
  temp2 <- rle(temp)
  x[which(temp == with(temp2, values[which.max(lengths)]))]
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

#create an observed over expected matrix
obs.exp.matrix <- function(mat, correct=F, outlier.correct = 0.995, lowess = F){
  pos <- which(mat$z > -1, arr.ind=T)
  val <- mat$z[mat$z > -1]
  d <- abs(pos[,1]-pos[,2])
  if(correct){
    norm.factor <- tapply(val, d, function(x) mean( x[x < quantile(x,
                                                          outlier.correct)] ) )
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

#calculate the compartment score
eigen.struct <- function( mat, outlier.correct = 0.995 ){
  oe <- obs.exp.matrix ( mat )
  #set non-finite values to 1
  oe$z[!is.finite(oe$z)] <- 1
  #remove outliers
  th <- quantile(oe$z, outlier.correct)
  oe$z[oe$z > th] <- th

  ev <- eigen(oe$z - 1)

  oe$ev1 <- ev$vector[,1]*(ev$value[1]**0.5)
  oe$ev2 <- ev$vector[,2]*(ev$value[2]**0.5)

  oe
}

switch.EV <- function( ev.data, chip, chrom ){
  start <- min(ev.data$x); end <- max(ev.data$x)
  sub.chip <- chip[chip[,1]==chrom & chip[,2] > start & chip[,2] < end,]
  window.cnt <- findInterval(sub.chip[,2], ev.data$x)
  window.cnt <- table(factor(window.cnt,1:length(ev.data$x)))
  #return true or false depending on whether the median value of
  #ChIP in the up or down compartment is highests (i.e. true if
  #down scores have the highest number of ChIP peaks)
  #up comp.                             down comp.
  median(window.cnt[ev.data$ev1 > 0]) < median(window.cnt[ev.data$ev1 < 0])
}
#===============================================================================
# core saddle function
#===============================================================================
saddle.core <- function(exp, chip, chrom, start, end, CS , nBins ){
  #select a region of the genome
  mat <- select.subset(exp, chrom, start, end)

  ev <- NULL
  if(is.null(CS)){
    #calculate the compartment score
    ev <- eigen.struct( mat )

    #orient the eigen vector properly for A and B
    if(!is.null(chip)){
      if(switch.EV( ev, chip, chrom ) ){
        ev$ev1 <- -ev$ev1
      }
    }
    ev <- ev$ev1
  } else {
    ev <- CS
  }

  #create the breaks for the eigen vectors
  # 0 is B, 1 is A
  br <- quantile( ev, seq(0,1,len=nBins+1) )
  bins <- cut(ev, br, include.lowest=T)

  unique.bins <- sort(unique(bins))

  #create a matrix of n x n bins
  ave.matrix <- matrix(0, ncol=nBins, nrow=nBins )

  oe.mat <- obs.exp.matrix( mat )
  oe.mat$z[!is.finite(oe.mat$z)] <- 1
  #loop over the bins
  for( i in 1:nBins ){
    for( j in 1:nBins ){
      sel.x <- which(bins==unique.bins[i]); sel.y <- which(bins==unique.bins[j])
      ave.matrix[i,j] <- mean(oe.mat$z[sel.x,sel.y])
    }
  }

  # return list of ev/cs and matrices in three-column form
  MAT       <- reshape2::melt(ave.matrix)
  MAT$chrom <- chrom

  ev        <- sort(ev)
  EV        <- data.frame(bin = 1:length(ev), ev = ev)
  EV$chrom  <- chrom
  return(list(MAT = MAT, EV = EV ))
}



#===============================================================================
# chromosome-wide wrapper function (run saddle.core for one or both arms)
#===============================================================================

saddle.chr <- function(exp, chip, chrom, CS, nBins){

  out     <- NULL
  out$MAT <- NULL
  out$EV  <- NULL
  ##################################################################    get arms
  chromStructure = NULL
  if(is.null(exp$CENTROMERES)){
    mat <- selectData_saddle( exp, chrom, chrom)

    cent <- which(apply(mat$z,1,sum)==0)
    cent <- largest.stretch(cent)

    chromStructure <- data.frame(chrom  = chrom,
                                 start  = 0,
                                 startC = min(cent)*exp$RES,
                                 endC   = max(cent)*exp$RES,
                                 end    = max(exp$ABS[exp$ABS$V1 == chrom, 3]))

  } else if(!is.null(exp$CENTROMERES) & chrom %in% exp$CENTROMERES[,1]){
    centChrom = exp$CENTROMERES[exp$CENTROMERES[,1] == chrom,]
    chromStructure <- data.frame(chrom  = chrom,
                                 start  = 0,
                                 startC = centChrom[,2],
                                 endC   = centChrom[,3],
                                 end    = max(exp$ABS[exp$ABS$V1 == chrom, 3]))
  } else {
    stop("Chromosome is not found in exp$CENTROMERES.\nQuitting...")
  }

  # check if an arm is smaller than Nbins*10 (eg: at least 10 bins per bin)
  # 1
  doArm1 = T
  if(!(((chromStructure[,3] - chromStructure[,2])/exp$RES)/nBins) > 10){
    doArm1 = F
  }
  doArm2 = T
  if(!(((chromStructure[,5] - chromStructure[,4])/exp$RES)/nBins) > 10){
    doArm2 = F
  }

  #################################################################    run arm 1
  if(doArm1){
    ## get CS of arm?
    ev <- NULL
    if(!is.null(CS)){
      CS_C  <- CS[CS[,1] == chrom,]
      CS_CA <- CS_C[CS_C[,3] <= chromStructure$startC, 4]
      # find out if the length of CS_CA is the same as the #bins in the arm-mat
      len_abs <- length(exp$ABS[exp$ABS$V1 == chrom &
                                  exp$ABS$V3 <= chromStructure$startC,4])
      len_CS  <- length(CS_CA)
      if(len_abs == len_CS){
        # cool. the lengths are the same. the same resolution is used.
        ev <- CS_CA
      } else {

        # try to find the end of arm of to that of CS
        altEnd = max(CS_C[CS_C[,3] <= chromStructure$startC, 2])

        CS_C  <- CS[CS[,1] == chrom,]
        CS_CA <- CS_C[CS_C[,3] <= altEnd, 4]
        # find out if the length of CS_CA is the same as the #bins in the arm-mat
        len_abs <- length(exp$ABS[exp$ABS$V1 == chrom &
                                    exp$ABS$V3 <= altEnd,4])
        len_CS  <- length(CS_CA)
        if(len_abs == len_CS){
          # cool. the lengths are now the same. the same resolution is used.
          ev <- CS_CA
        } else {
          warning('The amount on CS-scores and the amount of Hi-C bins
                  is not the same. \nEither call CS with the same resolution
                  or set CS to NULL.\nContinuing with NULL.')
          ev <- NULL
        }
      }
    } else {
      ev <- NULL
    }
    S_out <- saddle.core(exp = exp,
                         chip = chip,
                         chrom = chrom,
                         start = chromStructure$start,
                         end = chromStructure$startC,
                         CS = ev, nBins = nBins)
    S_out$MAT$arm <- "p"
    S_out$EV$arm  <- "p"

    out$MAT <- rbind(out$MAT, S_out$MAT)
    out$EV  <- rbind(out$EV, S_out$EV)
  }

  #################################################################    run arm 2
  if(doArm2){
    ## get CS of arm?
    ev <- NULL
    if(!is.null(CS)){
      CS_C  <- CS[CS[,1] == chrom,]
      CS_CA <- CS_C[CS_C[,2] >= chromStructure$endC &
                      CS_C[,3] <= chromStructure$end, 4]
      # find out if the length of CS_CA is the same as the #bins in the arm-mat
      len_abs <- length(exp$ABS[exp$ABS$V1 == chrom &
                                  exp$ABS$V2 >= chromStructure$endC &
                                  exp$ABS$V3 <= chromStructure$end,4])
      len_CS  <- length(CS_CA)
      if(len_abs == len_CS){
        # cool. the lengths are the same. the same resolution is used.
        ev <- CS_CA
      } else {

        # try to extend start of arm to that of CS
        altStart = min(CS_C[CS_C[,2] >= chromStructure$endC &
                              CS_C[,3] <= chromStructure$end, 2])
        CS_C  <- CS[CS[,1] == chrom,]
        CS_CA <- CS_C[CS_C[,2] >= altStart &
                        CS_C[,3] <= chromStructure$end, 4]
        # find out if the length of CS_CA is the same as the #bins in the arm-mat
        len_abs <- length(exp$ABS[exp$ABS$V1 == chrom &
                                    exp$ABS$V2 >= altStart &
                                    exp$ABS$V3 <= chromStructure$end,4])
        len_CS  <- length(CS_CA)

        if(len_abs == len_CS){
          # cool. the lengths are the same. the same resolution is used.
          ev <- CS_CA
        } else {warning('The amount on CS-scores and the amount of Hi-C bins
                is not the same. \nEither call CS with the same resolution
                or set CS to NULL.\nContinuing with NULL.')
        ev <- NULL
        }
      }
    } else {
      ev <- NULL
    }
    S_out <- saddle.core(exp = exp,
                         chip = chip,
                         chrom = chrom,
                         start = chromStructure$endC,
                         end = chromStructure$end,
                         CS = ev, nBins = nBins)

    S_out$MAT$arm <- "q"
    S_out$EV$arm  <- "q"

    out$MAT <- rbind(out$MAT, S_out$MAT)
    out$EV  <- rbind(out$EV, S_out$EV)
    }

  # MAT <- out$MAT
  # MAT[MAT$value < 0.5, "value"] <- 0.5
  # MAT[MAT$value > 2, "value"] <- 2
  # ggplot(MAT, aes(x = Var1, y = Var2, fill = value)) +
  #   geom_tile() +
  #   facet_grid(chrom ~ start+end) +
  #   scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b',
  #                        midpoint = 1) + coord_equal(expand = F) +
  #   GENOVA_THEME() + labs(x = "B <> A", y = "B <> A")

  # return list of ev/cs and matrices in three-column form
  return(list(MAT = out$MAT, EV = out$EV ))
}
#===============================================================================
# genome-wide wrapper function
#===============================================================================

#' Compute compartment-v-compartment scores
#'
#' Splits the range of compartment-scores in nBins bins (bin one having the highest compartment-scores and bin n the lowest).
#' Generates an average O/E contact-score for each compartment-score bin.
#' Sorts matrix on compartment-score.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chip BED-dataframe containing active sites (e.g. H3K27ac-peaks). If CS is NULL, this should be filled in!
#' @param chromsToUse Do the computation for a subset of chromosomes.
#' @param nBins The number of bins to split the compartment-score.
#' @param CS Use a bedgraph-df of compartment-scores instead of creating it on the fly.The resolution should be the same as the Hi-C (therefore, use compartment.score()). Using this will speed things up greatly.
#' @note
#' A binsize of 5 will produce a plot similar to FLyamer et al. (2017), while 100 will produce a plot similar to Bonev et al. (2017). The increase in time with more nBins is in exponential time at this point in time. Future versions will tackle this.
#' @examples
#' # run saddle() on all chromosomes with 25 bins.
#' saddle_WT = saddle(exp = Hap1_WT_1MB, chip = H3K27acPeaks, nBins = 25)
#'
#' # plot saddle-plot
#' visualise.saddle(list(saddle_WT), crossLines = T, addText = T)
#'
#' # plot compartment-strengths
#' visualise.compartmentStrength(list(saddle_WT))
#'
#' @return A log2(O/E) matrix and a DF of compartment-scores.
#' @import data.table
#' @export
saddle <- function(exp, chip = NULL, CS = NULL, chromsToUse = NULL, nBins = 10){
  if(!is.null(CS)){                                                         # CS
    if(!is.null(chip)){                                                # CS chip
      stop("Either CS or chip must be used!")
    } else {                                                           #CS !chip
      # no chip, so CS will just be used
    }

  } else {                                                               # ! CS
    if(is.null(chip)){
      stop("Either CS or chip must be used!")
    } else {
      # no CS, so CS will just be made
    }
  }

  if(!is.null(CS)){
    CS <- as.data.frame(CS)
  }
  if(is.null(chromsToUse)){
    chromsToUse <- exp$CHRS
  }

  out <- NULL
  for(C in chromsToUse){
    tmp = saddle.chr(exp = exp, chip = chip, chrom = C,  CS = CS, nBins = nBins)
    out$MAT <- rbind(out$MAT, tmp$MAT)
    out$EV <- rbind(out$EV, tmp$EV)

  }

  out$MAT$sample <- exp$NAME
  out$MAT$color <- exp$COL
  out$EV$sample <- exp$NAME
  out$EV$color <- exp$COL
  # return list of ev/cs and matrices in three-column form
  return(list(MAT = out$MAT, EV = out$EV ))
}


#' Plot compartment-strengths
#'
#' Takes the input from saddle() and produces a boxplot of compartment-strenghts.
#'
#' @param SBoutList A list of outputs trom saddle().
#' @param showInteractions Instead of compartment-strength, plot the average AA, BB and AB interactions.
#' @return A plot and an invisible underlying dataframe containing per-sample and -chromosome-arm the compartment-score.
#' @examples
#' # run saddle() on all chromosomes with 25 bins.
#' saddle_WT = saddle(exp = Hap1_WT_1MB, chip = H3K27acPeaks,nBins = 25)
#'
#' # plot compartment-strengths
#' CS_out = visualise.compartmentStrength(list(saddle_WT))
#'
#' @export
visualise.compartmentStrength = function(SBoutList, showInteractions = F){
  require(ggplot2)
  strengthDF = data.frame()
  interactionDF = data.frame()
  namesVector <- c()
  for(i in 1:length(SBoutList)){
    dat = SBoutList[[i]]
    namesVector <- c(namesVector, unique(dat$MAT$sample))
    dat$MAT$CC = 'XX'
    MAXbin = max(dat$MAT$Var1)
    # how many binsare in 20%
    binsTOse = floor(MAXbin * .2)
    binsTOse = max(1, binsTOse)
    dat$MAT$unLog = 2 ** dat$MAT$value
    dat$MAT[dat$MAT$Var1 <= binsTOse & dat$MAT$Var2 <= binsTOse,"CC"] = "BB"
    dat$MAT[dat$MAT$Var2 <= binsTOse & dat$MAT$Var1 >= MAXbin-binsTOse+1,"CC"] = "AB"
    dat$MAT[dat$MAT$Var1 >= MAXbin-binsTOse+1 & dat$MAT$Var2 >= MAXbin-binsTOse+1,"CC"] = "AA"
    dat$MAT = dat$MAT[dat$MAT$CC != 'XX',]

    tmp = dplyr::summarise(dplyr::group_by(dat$MAT,
                                           color,
                                           sample,
                                           chrom,
                                           arm,
                                           CC),score = mean(unLog))
    interactionDF = rbind(interactionDF, as.data.frame(tmp))
    for(S in unique(tmp$sample)){

      # get chroms in this sample
      meta = unique(tmp[tmp$sample == S,c("chrom","arm")] )

      for(C in unique(meta$chrom)){


        for(A in unique(unname(unlist(meta[meta$chrom == C,"arm"])))){
          tmpi = tmp[tmp$sample == S & tmp$chrom == C & tmp$arm == A, ]
          strength = log(tmpi[tmpi$CC == 'AA',"score"] * tmpi[tmpi$CC == 'BB',"score"] / tmpi[tmpi$CC == 'AB',"score"]**2)
          strengthDF = rbind(strengthDF, data.frame(S, C, A, strength, unique(tmp[tmp$sample == S, "color"])))
        }
      }
    }

  }
  RANGE = range(strengthDF$score, na.rm = T)
  RANGE[1] = floor(RANGE[1]) *.95
  RANGE[2] = ceiling(RANGE[2])*1.05




  if(showInteractions != TRUE){
    coltmp <- col2rgb(levels(as.factor(strengthDF$color)), alpha = T)/255
    coltmp[4, ] <- coltmp[4,] *0.85
    cols <- rgb(red = coltmp[1,], green = coltmp[2,], blue = coltmp[3,], alpha = coltmp[4,])

    boxplot(split(strengthDF$score,
                  strengthDF$S),
            ylim = RANGE, col = cols,
            ylab = "compartment strength")
  } else {

    interactionDF$sample <- factor(interactionDF$sample, levels = namesVector)
    P <- ggplot2::ggplot(interactionDF, ggplot2::aes(x = sample,
                                                     y = score,
                                                     fill = color)) +
      ggplot2::facet_wrap("CC") +
      ggplot2::geom_hline(yintercept = 1, lty = 3) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(y = "Average O/E") +
      GENOVA_THEME() +
      ggplot2::scale_fill_identity()+
      ggplot2::guides(fill = F) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
    print(P)
  }

  invisible(list(strengthDF, interactionDF))
}

rotate <- function(x) t(apply(x, 2, rev))

#' Plot saddle-plots
#'
#' Takes the input from saddle() and produces saddleplots.
#' Also plot a lowess-smoothed summary track of the sorted compartment-scores on top.
#'
#' @param SBoutList A list of outputs trom saddle().
#' @param addText Add annotations in the corners for legibility.
#' @param crossLines Plot crosslines in matrix. Crossline denote border of- A and B-bins. If not set, will plot them when nBins  >= 10.
#' @param zlim Set the z-lims of the matrix.
#' @param EVlim Set the z-lims for the average compartment-scores on top.
#' @param square Set this to FALSE if plots are unaligned: will not produce a rigidly square matrix.
#' @return A saddle-plot
#' @examples
#' # run saddle() on all chromosomes with 25 bins.
#' saddle_WT = saddle(exp = Hap1_WT_1MB, chip = H3K27acPeaks, nBins = 25)
#'
#' # plot saddle-plot
#' visualise.saddle(list(saddle_WT), crossLines = T, addText = T)
#' @export
visualise.saddle = function(SBoutList, addText = T,
                            zlim = c(0.5, 2), EVlim = c(-1.5, 1.5),
                            square = T, crossLines = NULL){
  #par_temp = par()
  if(is.null(EVlim)){
    EVlim = zlim
  }

  df_mat = NULL
  df_ev = NULL
  NSAMPLE = length(SBoutList)
  SAMPLENAMES = c()
  for(i in 1:NSAMPLE){
    SAMPLENAMES = c(SAMPLENAMES, as.character(unique(SBoutList[[i]]$MAT$sample) ))
    df_mat = rbind(df_mat, SBoutList[[i]]$MAT)
    df_ev = rbind(df_ev, SBoutList[[i]]$EV)
  }

  if(is.null(crossLines)){
    if(max(df_mat$Var1) < 10){
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

  ##### EV-plot
  # here, i get a lowess-smoothed plot of nBins mean compartment-scores
  # this allows us to not overcome problems with chromosome-arm-length

  newMat <- NULL
  df_ev$ID <- factor(apply(df_ev[,3:5], 1, paste0, collapse = "_"))
  for(i in unique(df_ev$ID)){
    tt <- df_ev[df_ev$ID == i,]
    ttscaled <- as.data.frame(approx(tt$ev, n = 100))
    ttscaled <- suppressWarnings(cbind(ttscaled, unique(df_ev[df_ev$ID == i, c("chrom","arm", "sample", "ID")])))
    ttscaled$x <- 1:100
    newMat <- rbind(newMat, ttscaled)
  }

  toPlotEV <- NULL
  for(s in unique(newMat$sample)){
    tmp <- newMat[newMat$sample == s,]
    loesje <- as.data.frame(loess.smooth(x = tmp$x, y = tmp$y))
    loesje$sample = s
    toPlotEV <- rbind(toPlotEV, loesje)
  }
  toPlotEV <- toPlotEV[,c(3,1,2)]
  #forEVPLot = aggregate(df_ev$ev, by = list(df_ev$sample, df_ev$bin), mean)

  for(S in SAMPLENAMES){
    tmp = setNames(toPlotEV[toPlotEV[,1] == S, 2:3], c("x","y"))



    tmp.bk = tmp
    crossPoint = approx(y = tmp.bk$x, x = tmp.bk$y, xout = 0)$y
    EV0 = c(EV0, crossPoint/max(toPlotEV[,2]))
    X <- tmp.bk$x
    y.low <- rep(0, length(tmp.bk$x))
    y.high <- tmp.bk$y

    plot(X,y.high,type = 'n', ylim = EVlim, axes = F, main = S, ylab = "", xlab = "")
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
  PAL = colorRampPalette(rev(c("#B2182B", "white", "#2166AC" )))
  for(Si in 1:NSAMPLE){
    S = SAMPLENAMES[Si]
    tmp = df_mat[df_mat$sample == S,]


    tmp = setNames(aggregate(tmp$value, by = list(tmp$Var1, tmp$Var2), mean, na.rm = T),
                   c("x",'y','z') )
    M = matrix(0, nrow = max(tmp$x), ncol = max(tmp$y))
    M[cbind(tmp$x,tmp$y)] = tmp$z

    M[M < zlim[1]] = zlim[1]
    M[M > zlim[2]] = zlim[2]

    image(M, zlim = zlim, col = PAL(50), axes = F,
          breaks = 2**(seq(log2(zlim[1]), log2(zlim[2]), length.out = 51))); box()


    if(addText){
      text(x = .05, y = .95, label = 'AB')
      text(y = .05, x = .95, label = 'BA')
      text(y = .05, x = .05, label = 'BB')
      text(y = .95, x = .95, label = 'AA')
    }
    if(crossLines){
      abline(v = (EV0[Si])-(.5/max(tmp$y)) ,lty = 2)
      abline(h = ((EV0[Si])-(.5/max(tmp$y))) ,lty = 2)
    }

  }
  if(square){
    par(pty = 'm')
  }
  #par(par_temp)

}


