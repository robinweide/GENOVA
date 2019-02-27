ATAlist <- ATAlist_WT

getTADstats <- function(MAT, NAasZero = T, pixWidth = 1, donutWidth = 5, pseudoCount = NULL){
  # store rows that are completely empty: these should not by set to 0 in NAasZero
  NArows = which(rowSums(MAT) == 0)
  if(NAasZero){
    MAT[is.na(MAT)] <- 0
  }
  
  MAT[,NArows]  <- NA
  MAT[NArows,]  <- NA
  
  # count empty rows
  Nempty = length(NArows)
  
  if(is.null(pseudoCount)){
    pseudoCount <- min(MAT[MAT > 0 & !is.na(MAT)] )
  }
  
  MAT[MAT == 0 & !is.na(MAT)] <- pseudoCount
  
  #############################################################################
  # get cornerPixel
  #############################################################################
  
  CPloc = c(1+(nrow(MAT)/4), 
            (nrow(MAT)/4)*3)
  
  CPlocDown = (CPloc[1] - ((pixWidth-1)/2) ) : (CPloc[1] + ((pixWidth-1)/2))
  CPlocUp = (CPloc[2] - ((pixWidth-1)/2) ) : (CPloc[2] + ((pixWidth-1)/2))
  
  cornerPixMean =  mean(MAT[CPlocDown,  CPlocUp], na.rm = T)
  
  #############################################################################
  # get donut: split into in- and outside TAD
  #############################################################################
  donutlocDown = (CPloc[1] -((donutWidth-1)/2) ):
                 (CPloc[1] + ((donutWidth-1)/2))
  donutlocUp = (CPloc[2] -((donutWidth-1)/2) ):
               (CPloc[2] + ((donutWidth-1)/2))
  
  donutlocDownInside = donutlocDown[donutlocDown > max(CPlocDown)]
  donutlocDownOutside = donutlocDown[donutlocDown < min(CPlocDown)]
  donutlocUpInside = donutlocUp[donutlocUp < max(CPlocUp)]
  donutlocUpOutside = donutlocUp[donutlocUp > min(CPlocUp)]
  
  inside = mean(MAT[donutlocDownInside, donutlocUpInside], na.rm = T)
  outside = mean(MAT[donutlocDownOutside, donutlocUpOutside], na.rm = T)
  
  if(sum(complete.cases(as.vector(MAT[donutlocDownInside, donutlocUpInside]))) == 0){
    inside = 0
  }
  
  if(sum(complete.cases(as.vector(MAT[donutlocDownOutside, donutlocUpOutside]))) == 0){
    outside = 0
  }
  
  #############################################################################
  # get outside: pixel shifted halfway
  #############################################################################
  midLoc <- ((nrow(MAT)/2) - ((donutWidth-1)/2)) : 
            ((nrow(MAT)/2) + ((donutWidth-1)/2))
  upLoc <- 1:donutWidth
  downLoc <- nrow(MAT):(nrow(MAT)-donutWidth)
  
  upDonut = MAT[midLoc, upLoc]
  downDonut = MAT[downLoc, midLoc]
  
  outsidePixel = mean(c(upDonut, downDonut), na.rm = T)
  
  #############################################################################
  # get surface FC: the middle disagonals in- and outside
  #############################################################################
  sMAT = MAT
  
  #remove outer TADs
  sMAT[1:CPloc[1], 1:CPloc[1]] <- NA
  sMAT[CPloc[2]:nrow(MAT), CPloc[2]:nrow(MAT)] <- NA
  
  # remove above TAD
  sMAT[1:CPloc[1], CPloc[2]:nrow(MAT)] <- NA
  sMAT[ CPloc[2]:nrow(MAT), 1:CPloc[1]] <- NA
  
  # melt
  sMAT <- reshape2::melt(sMAT)
  sMAT = sMAT[sMAT[,1] <= sMAT[,2], ]
  sMAT <- sMAT[!is.na(sMAT$value),]
  sMAT$d = abs(sMAT$Var1 - sMAT $Var2)
  
  surface = 0
  if(length(which(!is.na(sMAT$value))) == 0){
    # only NA's. 
    surface = 0
  } else {
    sMAT$loc = 'outside'
    sMAT[sMAT$Var1 >= CPloc[1] & sMAT$Var1 <= CPloc[2] & 
           sMAT$Var2 >= CPloc[1] & sMAT$Var2 <= CPloc[2], "loc"] <- 'inside'
    
    # throw away diags with less then 10 observations
    tmpi = table(sMAT[sMAT$loc == 'inside',]$d) > 10
    tmpo = table(sMAT[sMAT$loc == 'outside',]$d) > 10
    
    tmpi = as.numeric(names(tmpi)[tmpi])
    tmpo = as.numeric(names(tmpo)[tmpo])
    sMAT = sMAT[sMAT$d %in% tmpi[tmpi %in% tmpo],]
    
    
    OUTSIDE = sMAT[sMAT$loc == 'outside',]
    INSIDE = sMAT[sMAT$loc == 'inside',]
    
    dat = data.frame(IN = tapply(INSIDE$value, INSIDE$d, mean, na.rm = T),
                     OU = tapply(OUTSIDE$value, OUTSIDE$d, mean, na.rm = T))
    
    dat$TADstrength = log2( dat$IN / dat$OU )
    surface = mean(dat$TADstrength, na.rm = T)
  }
  

  
  #############################################################################
  # get output
  #############################################################################    
  datOut = data.frame(pixel = cornerPixMean,
                      insideDonut = inside,
                      outsideDonut = outside,
                      #pixelStrength = log((cornerPixMean ** 2)/(inside * inside)),
                      pixelStrength = log((cornerPixMean**2)/(inside * outside)),
                      #cornerStrength = log((cornerPixMean ** 2)/(inside * inside)),
                      #outsidePixel = outsidePixel,
                      TADstrength = surface,
                      numberOfEmptyBins = Nempty)

  return(datOut)
}



quantifyATA <- function(ATAlist, NAasZero = T, pixWidth = 1, donutWidth = 5){
  

  outDF <- list() # make a df with a line per loop (add column for color and name)
  
  resOut <- NULL
  for(i in 1:length(ATAlist)){
    ATAout <- ATAlist[[i]]$STACK.list

    avgsPix   <- lapply(ATAout, getTADstats, NAasZero = NAasZero, 
                        pixWidth = pixWidth, donutWidth = donutWidth)
    
    avgsPix   <- data.table::rbindlist(avgsPix)
    avgsPix <- as.data.frame(avgsPix)

    TADs <- ATAlist[[i]]$TADs.bed
    df     <- cbind(data.frame(name = names(ATAlist)[i],
                         TADIDX = 1:length(ATAout)),
                    TADs)
    df <- cbind(df, avgsPix)
    outDF[[names(ATAlist)[i]]] <- df
  }
  
 return(data.table::rbindlist(outDF)) 
}


