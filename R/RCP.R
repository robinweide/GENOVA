
#' Relative Contact Probabilities
#'
#' Produce a dataframe with the probabilities of contacts in a set of distance-bins.
#' Bins are created on a log scale, which leads to equal amounts of datapoints per bin.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param explist List of  GENOVA contacts-objects from `loadContacts()`.
#' @param bedlist A named list of BED-like dataframes of the regions of interest. 
#' @param chromsToUse A vector containing the chromosome-names of interest. 
#' @param maxDistance The maximal distance to calculate RCPs for.
#' @param genomeWide Normalise with genome-wide counts or per-chromosome (default).
#' 
#' @return An RCP_discovery object containing:
#' @return \item{raw}{a per-distance probability}
#' @return \item{smooth}{a log10-mean smoothed probability}
#' 
#' @section Resolution recommendation: 40kb-500kb
#' 
#' @examples
#' # Calculate the RCP of chromosome 1
#' \dontrun{
#' RCP_out = RCP(experimentList = list('WT' = WT_1MB), chromsToUse = 'chr1')
#'
#' # Plot the RCP
#' visualise(RCP_out)
#' }
#' @export
RCP = function(explist, bedlist = NULL, chromsToUse = NULL, maxDistance = NULL, genomeWide = NULL){
  
  data.table::setDTthreads(threads = 1)
  ##############################################################################
  ############################################################ Verify experiment 
  ##############################################################################
  if (any(c("MAT", "IDX") %in% names(explist))) {
    explist <- list(explist)
  }
  
  
  if(!is.null(bedlist)){
    if(!inherits(bedlist, 'list')){
      
      if(inherits(bedlist, 'data.frame')){
        
        bedlist = list(bedlist)
      }
      
    }
  }
  
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)

  ##############################################################################
  ######################################################################## setup
  ##############################################################################
  if (is.null(chromsToUse)) { # set the chromosomes to use
    if(is.null(genomeWide)){genomeWide = T}
    chromsToUse <- unique(as.vector(sapply(explist, function(x){x$CHRS})))
    
    chromsToUse = unique(unlist(chromsToUse))
    
  } else{
    if(is.null(genomeWide)){genomeWide = F} 
  }
  
  
  ##############################################################################
  ################################################################## maxDistance
  ##############################################################################
  
  largestChrom <- max(explist[[1]]$IDX[ explist[[1]]$IDX$V1 %in% chromsToUse ,3])
  if (is.null(maxDistance)) {
    maxDistance <- largestChrom
  }
  maxDistance = 10^(ceiling(log10(maxDistance)))

    ##############################################################################
  ##################################################################### main RCP
  ##############################################################################
  
  # choose between reiong or chom-based
  if(is.null(bedlist)){
    RCP_out = RCPchrom(explist, chromsToUse, genomeWide)
  } else {
    RCP_out = RCPbed(explist, bedlist, chromsToUse)
  }
  RCP_out$samplename = factor(RCP_out$samplename, 
                              levels = unname(unlist(lapply(explist, function(x){attr(x, "samplename")}))))
  
  ##############################################################################
  ####################################################################### smooth
  ##############################################################################
  
  breaks <- 10**seq(from = log10(max(sapply(explist, function(x){attr(x, 'resolution')}))),
                    to = log10(maxDistance), length.out = 101)
  breaks = c(0,breaks)
  
  splitted = split(RCP_out, list(RCP_out$region, RCP_out$samplename))
  smthd = lapply(splitted ,function(x){
    cuts = cut(x$distance, breaks, include.lowest = T)
    splt = split(x, cuts, drop = T)
    
    lp = lapply(splt ,function(chunk){
      data.frame(distance = mean(chunk$distance), P = mean(chunk$P))
    })
    lp = rbindlist(lp)
    
    lp$samplename = x[1,'samplename']
    lp$colour = x[1,'colour']
    lp$region = x[1,'region']
    lp[,c(5,1,2,3,4)]
  })
  smthd = rbindlist(smthd)
  smthd[smthd$distance == min(smthd$distance), 'distance'] = 0
  ##############################################################################
  ########################################################################## out
  ##############################################################################
  
  RCP_out = structure(
    list('raw' = RCP_out, smooth = smthd),
    class = c("RCP_discovery", "discovery"), 
    package = "GENOVA",
    resolution = resolution(explist)[[1]]
  )
  
  if(!is.null(bedlist)){
    attr(RCP_out, 'norm') = 'local'
  } else  if(genomeWide){
    attr(RCP_out, 'norm') = 'genome-wide'
  } else {
    attr(RCP_out, 'norm') = 'chromosome'
  }
  
  return(RCP_out)
  
}

RCPchrom = function(explist, chromsToUse, genomeWide){
  
  # init
  .        <- NULL
  V1       <- NULL
  V3       <- NULL
  V4       <- NULL
  C        <- NULL
  distance <- NULL
  
  RCP_out = lapply(explist, function(x){
    
    SIG = x$MAT
    
    chromRange = x$IDX[ , .(first = min(V4)), by = V1]
    chromRange = chromRange[order(chromRange$first),]
    
    F1 = findInterval(SIG$V1, chromRange$first)
    F2 = findInterval(SIG$V2, chromRange$first)
    
    F1[F1 != F2] = NA
    
    SIG$C = chromRange$V1[F1]
    SIG = SIG[!is.na(SIG$C),]
    
    D = abs(SIG$V2-SIG$V1) * attr(x, 'resolution')
    
    VALS = NULL
    if(genomeWide){
      
      SIG$D = D
      VALS = SIG[ , sum(V3, na.rm = T), by = list(D)]

      VALS$P = VALS$V1/sum(SIG$V3)
      VALS$V1 = NULL
      colnames(VALS) = c('distance', 'P')
      VALS$region = 'genome-wide'
    } else { # no GW
      
     SIG$D = D
     VALS = SIG[ , sum(V3, na.rm = T), by = list(C, D)]
     
      VALS = VALS[VALS$C %in% chromsToUse,]
      VALS$P = VALS$V1/sum(SIG$V3)
      VALS$V1 = NULL
      colnames(VALS) = c('region','distance', 'P')

    }
    
    VALS$colour = attr(x, 'colour')
    VALS$samplename = attr(x, 'samplename')

    
    data.table::setkey(VALS, distance)
    VALS
  })
  
  RCP_out = data.table::rbindlist(RCP_out)
  
  data.table::setkey(RCP_out, distance)

  RCP_out[,c(3,1,2,5,4)]
  

  
}



RCPbed = function(explist, bedlist, chromsToUse){
  
  # init
  .          <- NULL
  V1         <- NULL
  V3         <- NULL
  V4         <- NULL
  C          <- NULL
  distance   <- NULL
  D          <- NULL
  colour     <- NULL  
  samplename <- NULL
  SUM        <- NULL
  
  
  ############################################################################## for bed...
  out = lapply(bedlist, function(bed){
    
    if(!is.null(chromsToUse)){
      bed = bed[bed[,1] %in% chromsToUse,]
    }
    
    if(!is.numeric(bed[,2])){
      bed[,2] = as.numeric(bed[,2])
    }
    
    if(!is.numeric(bed[,3])){
      bed[,3] = as.numeric(bed[,3])
    }
    
    ############################################################################ for exp...
    exoOut = lapply(explist, function(x){
      
      IDX = split(x$IDX, x$IDX$V1)
      
      chromRange = x$IDX[ , .(first = min(V4)), by = V1]
      chromRange = chromRange[order(chromRange$first),]
      
      ########################################################################## for chrom...
      regionSIG = lapply(unique(bed[,1]), function(CHROM){
        
        bedi = bed[bed[,1] == CHROM,]
        
        # expand bedi to include all bins
        bedi[,2] = floor(bedi[,2]/attr(x, 'resolution'))*attr(x, 'resolution')
        bedi[,3] = ceiling(bedi[,3]/attr(x, 'resolution'))*attr(x, 'resolution')
        
        # do reduce on bed
        
        # get all unique bins
        bins = apply(bedi[,2:3] ,1, FUN = function(bx){
          seq(bx[1], bx[2], by = attr(x, 'resolution')) # so ugly, much wow
        })
        bins = unlist(bins)
        
        idxi = IDX[[CHROM]]
        idxi = idxi$V4[findInterval(bins, idxi$V2)]
      
        mats = x$MAT[list(idxi)]
  
        #now check if V2 is on the same chromosome
        mats = mats[findInterval(mats$V1, chromRange$first) == findInterval(mats$V2, chromRange$first),]
        
        mats
      })
      regionSIG = rbindlist(regionSIG)
      
      ##########################################################################
      regionSIG$D = abs(regionSIG$V2-regionSIG$V1) * attr(x, 'resolution')
      regionSIG$colour = attr(x, 'colour')
      regionSIG$samplename = attr(x, 'samplename')
      regionSIG$SUM = sum(regionSIG$V3)
      
      regionSIG
    })
    exoOut = rbindlist(exoOut)
    
    VALS = NULL
    VALS = exoOut[ , sum(V3, na.rm = T), by = list(D, colour, samplename, SUM)]

    VALS$P = VALS$V1/VALS$SUM
    VALS$V1 = NULL
    VALS$SUM = NULL
    colnames(VALS) = c('distance', 'colour','samplename','P')

    data.table::setkey(VALS, distance)
    
    VALS
  })
  
  if(!is.null(names(bedlist))){
    names(out) = names(bedlist) 
  } else {
    names(out) = paste0('bed',1:length(bedlist))
  }
  
  out = rbindlist(out, idcol = 'region')
  
  data.table::setkey(out, distance)
  out[,c(1:2,5,4,3)]
}

# Utils -------------------------------------------------------------------

# RAW RCP in, lfc out
#' @export
#' @keywords internal
#' @title RCP log2 foldchange
#' 
#'  RAW RCP in, lfc out
#'
#' @param dt a data.table of rcp
#' @param contrast the name of the contrast-sample
#' @param breaks a set of numbers to use as intervals
#'
#' @return a data.table with the log2 fold changes compared to the `contrast`.
#' @export
RCPlfc = function(dt, contrast, breaks){
  
  intervalMid = (diff(breaks)/2)+breaks[-length(breaks)]
  intervalMid[1] = 0
  
  # intersect with cuts
  CT = cut(dt$distance, breaks, labels = F, include.lowest = T)
  
  dt$cut = intervalMid[CT]
  
  SPREAD = dcast(dt, cut ~ samplename, value.var = "P", fun.aggregate = mean)
  
  i = which(colnames(SPREAD) == contrast)
  j = 2:ncol(SPREAD); j = j[j != i]
  
  out = lapply(j, function(J){
    log2(as.matrix(SPREAD[,J, with = F]) / as.matrix(SPREAD[,i, with = F]))
  })
  
  out = as.data.frame(do.call('cbind',out))
  
  out$distance <- SPREAD$cut
  
  setDT(out)
  out <- melt.data.table(out, id.vars = 'distance')
  setDF(out)
  colnames(out)[2:3] = c('samplename','P')
  out
}







# 
# 
# ! do
# 
# quantify.RCP_discovery = function(discovery, ...){
#   
#   
#   DISC = split(discovery, list(discovery$region, discovery$samplename))
#   # calculate the slope
#   tmp = lapply(DISC, function(dat){
#     slopevector = calculate.slope(cbind(dat$distance,dat$P))
#     slopevector$distance = 10**slopevector$pos
#     slopevector$region = unlist(dat[1,1])
#     slopevector$samplename =  unlist(dat[1,4])
#     slopevector$colour =  unlist(dat[1,5])
#     slopevector[,-1]
#     
#   })
#   
#   
# }
# 
# calculate.slope <- function( rcp, min.bin = 4, max.bin = 8, step = 0.2, bin.increase = 0.3 ){
#   add.step <- seq(0, 1, by=step)
#   add.step <- add.step[!add.step == 1] #remove the 1, because it is the same as 0
#   x.pos <- c()
#   slopes <- c()
#   x <- log10(rcp[,1]); y <- log10(rcp[,2]) #move into log-space
#   for( add in add.step ){
#     br <- seq( min.bin + add, max.bin + add, by = bin.increase ) #create the breaks
#     for( i in 2:(length(br)-1) ){ 
#       #select everything between the breaks
#       xs <- x[x > br[i] & x < br[i+1]]; ys <- y[x > br[i] & x < br[i+1]]; 
#       if(length(xs) > 0){
#         #print(lm(ys~xs)$coeff[2])
#         slopes <- c(slopes, lm(ys~xs)$coeff[2] )
#         x.pos <- c(x.pos, median(xs) )
#       }	
#       
#     }
#   }	
#   slope.df <- data.frame(pos=x.pos, slope=slopes)
#   slope.df <- unique(slope.df) #there can be double entries, so remove these
#   slope.df <- slope.df[order(slope.df[,1]),]
#   slope.df
# }
# 

# 
# 
# 
# 
# 
# 
# 
# 
# 
# !man export


