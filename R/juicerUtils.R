loadJuicer = function(juicerPath, resolution, scale_bp = 1e9, scale_cis = F, balancing = T){

  require(strawr)

  # get metadata of juicer-file
  juicer_metadata = get_juicer_metadata(juicerPath)

  # check if resolution is found!
  if(!any(resolution == juicer_metadata[[2]])){
    stop('resolution is not found.\nAvailable resolutions are ',
         paste(rev(juicer_metadata[[2]]), collapse = ', '))
  }
  juicer_metadata[[1]] = juicer_metadata[[1]][juicer_metadata[[1]]$chrom != 'All',]

  expandedChromosomes = as.data.frame(t(combn(juicer_metadata[[1]]$chrom, m = 2)), stringsAsFactors = F)
  expandedChromosomes = rbind(as.data.frame(cbind(juicer_metadata[[1]]$chrom,juicer_metadata[[1]]$chrom)),
                              expandedChromosomes)
  expandedChromosomes = apply(expandedChromosomes, 2, as.character)
  strawNorm = ifelse(balancing, 'KR', "NONE")
  
  juicerList = lapply(seq_len(nrow(expandedChromosomes)), function(eci){
    
    ec = expandedChromosomes[eci,]
    
    juicer_in <- tryCatch(
      {strawr::straw(norm = strawNorm,
                     fname =juicerPath,
                     chr1loc =  ec[2],
                     chr2loc =  ec[1],
                     unit = 'BP',
                     binsize = resolution)
      },error=function(cond) {
        
        return(NULL)
      })
    
    # remove NaN
    if(!is.null(juicer_in)){
      juicer_in = juicer_in[!is.nan(juicer_in[,3]),]
      if(nrow(juicer_in) == 0){
        juicer_in = NULL
      }
    }
    
    # decorate
    if(!is.null(juicer_in)){
      
      juicer_in$chrom_x = ec[1]
      juicer_in$chrom_y = ec[2]
      
      juicer_in = juicer_in[, c(4,1,5,2,3)]
    }
    juicer_in
    
  })


  juicer_data = data.table::rbindlist(juicerList)

  # split into index and signal-files
  SA = splitJuicerData(juicer_data, resolution)

  SIG = SA[[1]]
  ABS = SA[[2]]

  
  
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



get_juicer_metadata = function(juicerPath){
  file2read = file(juicerPath, "rb")
  number_of_integers_in_file = 125

  # check if magicstring is there
  magic_string = readBin(file2read, character(), n = 1, size = 3)
  if(magic_string != 'HIC'){
    stop('magic sting comparison failed. Check the .hic file!')
  }

  hic_version = readBin(file2read, integer(), n = 1, size = 4)

  master_index = readBin(file2read, integer(), n = 1, size = 8)

  genome = readBin(file2read, character())

  n_attributes = readBin(file2read, integer(), n = 1, size = 4)

  n_chrs = readBin(file2read, integer(), n = 1, size = 4)


  chrom_name_length = lapply(1:n_chrs, function(x){
    chrom = readBin(file2read, character(), size = 1)
    len = readBin(file2read, integer(), n = 1)
    data.table::data.table(chrom, len)
  })
  chrom_name_length = data.table::rbindlist(chrom_name_length)

  n_BP_Res  = readBin(file2read, integer(), n = 1, size = 4)


  resolutions = sapply(1:n_BP_Res, function(x){
    readBin(file2read, integer(), n = 1)
  })

  close(file2read)

  return(list(chrom_name_length, resolutions))
}


splitJuicerData = function(juicer_data, resolution){

  # get all unique bins
  bins_x = stats::setNames(unique(juicer_data[, 1:2]), c('V1', 'V2'))
  bins_y = stats::setNames(unique(juicer_data[, 3:4]), c('V1', 'V2'))

  # make abs
  ABS = unique(rbind(bins_x,bins_y))
  ABS$V3 = ABS$V2 + resolution
  data.table::setkeyv(ABS, c('V1','V2'))
  ABS$index = 1:nrow(ABS)

  # make sig
  SIG = stats::setNames(juicer_data, c('V1','V2','Y1','Y2','signal'))
  data.table::setkeyv(SIG, c('V1','V2'))

  # merge first
  SIG$index1 <- SIG[ABS, index, nomatch = 0]
  SIG[,1:2] = NULL
  colnames(SIG)[1:2] = c('V1','V2')
  data.table::setkeyv(SIG, c('V1','V2'))

  # second round
  SIG$index2 <- SIG[ABS, index, nomatch = 0]

  SIG = stats::setNames(SIG[,c(4,5,3)], paste0('V',1:3))

  colnames(ABS) = paste0('V', 1:4)
  return(list(SIG, ABS))
}


