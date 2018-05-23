#' Aggegrate Peak Analysis
#'
#' Get a z-stack matrix of 2D-regions (e.g. loops).
#' The function extracts matrices around a defined set of pixels, like possible
#' loops.
#' Next, it averages over all matrices to produce a single Z-stack matrix.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param experiment The Hi-C experiment object of a sample:
#' produced by construct.experiment().
#' @param loop.bed BEDPE-df containing the loop positions.
#' Should have at least six columns, containing the regions of both
#' loop-anchors: chrom1/start1/end1/chrom2/start2/end2.
#' @param smallTreshold The minimal size of loops. Too small loops lead to
#' messy plots, with the diagonal visible.
#' @param verbose Produces a progress-indication.
#' @param rmOutlier Perform outlier-correction by tresholding [outlierCutOff]
#' percentile per pixel.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21
#'  yield an output with 10 Hi-C bins up- and downstream of the anchor).
#' @param outlierCutOff The severity of outliers. We compute the [outlierCutOff]
#'  percentile per pixel and set values bigger than that to this value.
#' @examples
#' # Run APA on WT loops
#' APA_WT <- APA(experiment = WT_10kb, loop.bed = WT_Loops)
#'
#' # Visualise the APA-results
#' visualise.APA.ggplot(APAlist = list('WT' = APA_WT),
#' zTop = c(0,9.5), zBottom = c(-5,5),focus = 1)
#' @return \item{APA}{a matrix with the average Z-stack scores}
#' @return \item{rawMatList}{the raw underlying matrices per loop}
#' @return \item{APAraw}{the uncorrected APA_matrix}
#' @return \item{APAxy}{a version for base::image}
#' @return \item{RES}{the resolution of the used Hi-C matrix}
#' @return \item{OUTLIERCORRECTIONSWITCH}{the outlierCutOff if rmOutlier is
#' TRUE, otherwise set to FALSE}
#' @import data.table
#' @export
APA <- function(experiment, loop.bed, smallTreshold = NULL, rmOutlier = T,
                size = 21, verbose = F,  outlierCutOff = 0.995, ...){

  #########################
  # experiment operations #
  #########################

  hicdata <- experiment$ICE
  # Check for setkey
  if(is.null(data.table::key(hicdata))){data.table::setkey(hicdata, V1, V2)}

  bed <- experiment$ABS
  # Make chr:pos index of HiC-index
  bed.p <- paste0(bed[,1], ":", bed[,2])

  ########################
  # misc initialisations #
  ########################

  if(is.null(smallTreshold )){
    smallTreshold = experiment$RES*  (((size+2)) )
  }

  if(((size-1) /2 )%%1 != 0){stop("Size should be an even number +1")}
  size.offset = (size-1)/2

  resolution <- experiment$RES
  pos <- seq(-(size-1)/2, (size-1)/2)*resolution

  #####################
  #  set output-vars  #
  #####################
  OUTLIERCORRECTIONSWITCH = ifelse(rmOutlier,
                                   yes = T,
                                   no = outlierCutOff)
  APA = NULL
  APAraw = NULL
  APAxy = NULL
  RES =  experiment$RES
  rawMatList = list()

  #####################
  #   prepare loops   #
  #####################

  # Remove smaller loops
  loop.bed <- loop.bed[abs(loop.bed[,6]-loop.bed[,2]) >= smallTreshold ,]

  # Is BED 1 upstream of BED2?
  lb = loop.bed
  loop.bed[loop.bed[,2] > loop.bed[,5],] <- lb[lb[,2] > lb[,5],c(4,5,6,1,2,3)]

  # Prune NA-rows
  loop.bed <- na.exclude(loop.bed[,1:6])

  # Get start positions of the loop, floored to 10kb, and make chr:pos index
  # the tmp is just because Robin is anal about linelengths.
  tmp = ((abs(loop.bed[,3] - loop.bed[,2])/2) + loop.bed[,2])/resolution
  tmp = floor(tmp)
  loop.bed1.p <- paste0(loop.bed[,1], ":", as.integer(resolution*tmp))

  tmp = ((abs(loop.bed[,6] - loop.bed[,5])/2) + loop.bed[,5])/resolution
  tmp = floor(tmp)
  loop.bed2.p <- paste0(loop.bed[,4], ":",as.integer(resolution*tmp ))

  # Bugfix: in some cases paste0 yields an extra whitespace
  loop.bed1.p <- gsub(" ", "", loop.bed1.p, fixed = TRUE)
  loop.bed2.p <- gsub(" ", "", loop.bed2.p, fixed = TRUE)

  # Check if loops are all found:
  if(!all(loop.bed1.p %in% bed.p)){
    if(!all(loop.bed2.p %in% bed.p) ){
      msg = paste0('Not all loop-anchors can be found in HiC-file.\n\t',
                   'Are you sure that both are from the same reference?')
      warning(msg)
    }
  }

  #####################
  # prepare hi-c idx  #
  #####################

  # Get anchor HiC-indexes
  x.pos <- bed[match(loop.bed1.p,bed.p),4]
  y.pos <- bed[match(loop.bed2.p,bed.p),4]

  # check for empty indexes
  na.pos <- unique(c(which(is.na(x.pos)) ,which(is.na(y.pos)) ))
  if(length(na.pos) != 0){
    x.pos <- x.pos[-na.pos]
    y.pos <- y.pos[-na.pos]
  }


  ##########################
  # loop over loops (haha) #
  ##########################
  lx <- length(x.pos)
  for( i in 1:lx){
    if(verbose){cat(paste0(i, ' of ', lx, ' loops.'), "\r")}

    # Get Indeces of 10 up/downstream of anchor
    sel.x <- (x.pos[i]-size.offset):(x.pos[i]+size.offset)
    sel.y <- (y.pos[i]-size.offset):(y.pos[i]+size.offset)

    # Exract scores from HiC
    s <- select.sub.2D(hicdata, sel.x,sel.y)
    s$V1 <- (s$V1 - min(sel.x)) + 1
    s$V2 <- (s$V2 - min(sel.y)) + 1

    # Check for 21x21 mat:
    if(length(unique(s$V2)) < size){
      cols <- 1:size
      toADD <- cols[!cols %in% unique(s$V2)]
      newdf <- data.frame(unique(cbind(rep(1:size, each = size),
                                       rep(toADD, length = size),0)))
      colnames(newdf) <- c("V1","V2","V3")
      s <- rbind(s, newdf)
    }
    if(length(unique(s$V1)) < size){
      cols <- 1:size
      toADD <- cols[!cols %in% unique(s$V1)]
      newdf <- data.frame(unique(cbind(rep(toADD, length = size),
                                       rep(1:size, each = size),0)))
      colnames(newdf) <- c("V1","V2","V3")
      s <- rbind(s, newdf)
    }

    # Make matrix
    s.mat <- reshape2::acast(s, V1~V2, value.var="V3")
    s.mat[is.na(s.mat)] <- 0

    # Add to rawMatList
    rawMatList[[i]] <- s.mat

  }

  # Convert to 3D array
  rawMatList <- rawMatList[!unlist(lapply(rawMatList, is.null))]
  sm <- simplify2array(rawMatList)

  #####################
  #  outlier correct  #
  #####################
  if(rmOutlier){
    sm.bk = sm
    sm = outlier3Darray(ARRAY = sm, Q = outlierCutOff)
    # par(mfrow = c(1,2), pty = 's')
    # image(apply(sm.bk, c(1,2), mean) , zlim = c(0,50), main = 'raw')
    # image(apply(sm, c(1,2), mean) , zlim = c(0,50), main = 'corrected')
    APAraw = apply(sm.bk, c(1,2), mean)
    APA = apply(sm, c(1,2), mean)
  } else {
    APAraw = apply(sm, c(1,2), mean)
    APA = apply(sm, c(1,2), mean)
  }


  #####################
  #  misc operations  #
  #####################

  # Rotate matrices 90CW, so that diagonal of HiC-matrix is in bottom-left
  APAraw <- t(apply(APAraw, 2, rev))
  colnames(APAraw) <- 1:size

  APA <- t(apply(APA, 2, rev))
  colnames(APA) <- 1:size

  # make a image-ready list
  APAxy=list(x=pos, y=pos, z=APA)

  ##########################
  #       so long and      #
  #  thanks for the fish!  #
  ##########################

  return(list(APA = APA,
              rawMatList = rawMatList,
              APAraw = APAraw,
              APAxy= APAxy,
              RES = RES,
              OUTLIERCORRECTIONSWITCH = OUTLIERCORRECTIONSWITCH ))
}
