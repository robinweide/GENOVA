loadCooler = function(cooler, balancing = T, scale_bp = NULL, scale_cis = F){

  ABS = data.table::as.data.table(rhdf5::h5read(file = cooler,
                             name = "bins"))
  ABS$bin = 1:nrow(ABS)

  if(!'weight' %in% colnames(ABS)){
    balancing = F
    warning('No balancing is done due to missing weights.')
  }

  SIG = data.table::as.data.table(rhdf5::h5read(file = cooler,
                             name = "pixels"))

  rhdf5::h5closeAll()

  SIG[,1] = SIG[,1] + 1
  SIG[,2] = SIG[,2] + 1

  if(balancing){
    SIG = balance_cooler(ABS, SIG)
  }

  ABS$weight = NULL
  colnames(ABS) = paste0('V', 1:4)
  colnames(SIG) = paste0('V', 1:3)

  if (!is.null(scale_bp)) {
    
    if(scale_cis){
      
      chromRange = ABS[ , .(first = min(V4)), by = V1]
      chromRange = chromRange[order(chromRange$first),]
      
      F1 = findInterval(SIG$V1, chromRange$first)
      F2 = findInterval(SIG$V2, chromRange$first)
      
      SIG$V3 <- scale_bp * SIG$V3 / sum(SIG[ifelse(F1 == F2, T, F), 3])
      
    } else {
      SIG$V3 <- scale_bp * SIG$V3 / sum(SIG$V3)
    }

  }

  return(list(SIG, ABS))

}



#
# out = loadCooler("Dixon2012-H1hESC-HindIII-allreps-filtered.50kb.cool")
# SIG = out[[1]]
# ABS = out[[2]]



balance_cooler = function(ABS, SIG){

  WEIGHTS = data.table::data.table(stats::setNames(ABS[, c('bin','weight'),],
                                                   c('bin1_id','bin1_weight')),
                                   key = 'bin1_id')
  WEIGHTS$bin1_weight[is.nan(WEIGHTS$bin1_weight)] = 1

  data.table::setkey(SIG, 'bin1_id')
  SIG <- SIG[WEIGHTS, nomatch = 0]

  colnames(WEIGHTS) = c('bin2_id','bin2_weight')

  data.table::setkey(WEIGHTS, 'bin2_id')
  data.table::setkey(SIG, 'bin2_id')

  SIG <- SIG[WEIGHTS, nomatch = 0]

  SIG$weight = SIG$bin1_weight * SIG$bin2_weight
  SIG$signal = SIG$count*SIG$weight

  return(SIG[, c(1,2,7)])
}

#
#
#
#
# idChr10 =  (ABS[ABS[,1] == 'chr19',"bin"])
#
# data = SIG[SIG$bin1_id %in% idChr10 & SIG$bin2_id %in% idChr10, ]
#
# RAW = setNames(data[,1:3], c('bin1_id', 'bin2_id','signal')); RAW$set = 'count'
# BAL = data[,c(1:2,7)]; BAL$set = 'balanced'
#
# data = rbind(RAW,BAL)
#
# dataSwitch = data
# dataSwitch$bin1_id = data$bin2_id
# dataSwitch$bin2_id = data$bin1_id
#
# data = rbind(data, dataSwitch)
#
# data$set = factor(data$set, levels = c('count','balanced'))
#
# ggplot(data, aes(x = bin1_id, y = bin2_id, fill= log10(signal))) +
#   theme(panel.background = element_rect(fill = 'white'))+
#   geom_raster() +
#   facet_grid(.~ set) +
#   scale_fill_gradientn(colours = rev(c('black','red','orange','yellow', 'white')))
#
