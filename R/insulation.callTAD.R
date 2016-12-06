#' insulation.callTAD
#'
#' Call TADs
#'
#' @param exp1 The Hi-C experiment object: produced by construct.experiment().
#' @param BEDCOLOR Color of items in the resulting BEDPE-file
#' @return A BEDPE-df
#' @note Call insulation scores first and store these in exp$INSULATION
#' @export
insulation.callTAD <- function(exp,  BEDCOLOR = "127,201,127"){
  if(is.null(exp$INSULATION)){ stop("Call insulation score first and store in exp$INSULATION")}
  res <- exp$RES
  scooch <- floor(100e3 / res)
  entries <- list()
  df = NULL
  CHROMS <- unique(exp$INSULATION[,1])
  
  exp$INSULATION$V2 <- exp$INSULATION[,2] + exp$RES
  
  for(CCC in CHROMS){
    #car("Starting chromosome",CCC, "\n")
    INSU <- exp$INSULATION[exp$INSULATION[,1] == CCC ,]
    insCol <- INSU[,4]
    insCol[!is.finite(unlist(insCol))] <- 0
    INSU <- cbind(INSU[,1:3],scale(insCol, center = TRUE, scale = TRUE))
    colnames(INSU)[4] <- "V4"
    if( length(INSU$V4) == length(which(is.nan(INSU$V4)))){break}
    #car("Computing the delta-vector...\n")
    
    add <- 1:scooch
    i <- (scooch+1):(nrow(INSU)-scooch)
    i.rep <- rep(i, each=length(add))
    right <- matrix(INSU$V4[i.rep+add], ncol=scooch, byrow=T)
    left <- matrix(INSU$V4[i.rep-add], ncol=scooch, byrow=T)
    right <- apply(right, 1, mean)
    left <- apply(left, 1, mean)
    delta = left - right
    deltaDF <- data.frame(INSU[i,1:4], delta)
    colnames(deltaDF)[1:4] <- paste0("V", 1:4)
    colnames(deltaDF)[5] <- "delta"
    deltaDF <- dplyr::arrange(deltaDF,  V1,V2)
    deltaDF$ID <- 1:nrow(deltaDF)
    
    ####
    # First find peaks
    ###
    VALLEYS <- diff(c(.Machine$integer.max, deltaDF$V4)) > 0L
    VALLEYS <- cumsum(rle(VALLEYS)$lengths)
    VALLEYS <- VALLEYS[seq.int(1L, length(VALLEYS), 2L)]
    if (deltaDF$V4[[1]] == deltaDF$V4[[2]]) {
      VALLEYS <- VALLEYS[-1]
    }
    #car("Calling borders...\n")
    boundaryCalls <- NULL
    VALLEYS <- sort(unique(c(VALLEYS+1, VALLEYS, VALLEYS-1)))
    VALLEYS <- VALLEYS[VALLEYS > 2]
    for(j in VALLEYS[2:length(VALLEYS)]){
      dat <- deltaDF[(j-1):(j+1),]
      normalOrder <- order(dat$delta)
      dat0 <- dat
      dat0$delta[2] <- 0
      zeroOrder <- order(dat0$delta)
      if(! all(complete.cases(dat))){next}
      if(! all(normalOrder == c(3,2,1))){next}
      if(!all(zeroOrder == normalOrder)){next}
      UP <- max(deltaDF[(j-1),5],na.rm = T)
      DOWN <- min(deltaDF[(j+1),5],na.rm = T)
      if((UP - DOWN) < 0.1){next}
      if(is.null(boundaryCalls)) {
        boundaryCalls <- rbind(boundaryCalls,dat[2,])
      } else{
        nearbyPoints <- dplyr::filter(boundaryCalls, V2  >= (dat[2,2] - (res*scooch)) , V2 <= (dat[2,3] + (res*scooch) ))  # res*scooch
        if(nrow(nearbyPoints) == 0){
          boundaryCalls <- rbind(boundaryCalls,dat[2,])
        } else {
          closeToZeroDelta <- min(abs(c(nearbyPoints$delta, dat[2,]$delta)))
          winner <- rbind(nearbyPoints, dat[2,])[which(abs(c(nearbyPoints$delta, dat[2,]$delta)) == closeToZeroDelta),]
          losers <- rbind(nearbyPoints, dat[2,])[-which(abs(c(nearbyPoints$delta, dat[2,]$delta)) == closeToZeroDelta),]
          boundaryCalls <- rbind(boundaryCalls,dat[2,])
          boundaryCalls <- dplyr::anti_join(boundaryCalls, losers, by = c("V1", "V2", "V3", "V4", "delta"))
        }
      }
    }
    
    for(i in 1:nrow(boundaryCalls)){
      if( boundaryCalls[i,3] < boundaryCalls[i,2]  ){
        tmp3 <- boundaryCalls[i,3]
        tmp2 <- boundaryCalls[i,2]
        boundaryCalls[i,3] <- tmp2
        boundaryCalls[i,2] <- tmp3
      }
    }
    
    boundaryCalls <- boundaryCalls[with(boundaryCalls, order(V1, V2)), ]

    #car("Generating bedgraph...\n")
    
    for(i in 2:nrow(boundaryCalls)){
      if(!boundaryCalls[i-1,1] == boundaryCalls[i,1]){next}
      prev <- boundaryCalls[i-1,3]
      now <- boundaryCalls[i,1:2]
      now[,2] <- now[,2]
      prev <- prev
      ddd <-cbind(now, prev,now,prev, BEDCOLOR)[c(1,3,2,1,3,2,7)]
      colnames(ddd) <- c("a", 'b', 'c', 'd', 'e', 'f', 'g')
      df <- rbind(df, ddd)
    }
  }
  
  return(bedgraph = df)
}
