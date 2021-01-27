loadJuicer = function(juicerPath, resolution, scale_bp = 1e9, scale_cis = F, balancing = T){
  
  try_require("strawr", "loadJuicer", "github")
  
  juicerPath <- normalizePath(juicerPath)
  # get metadata of juicer-file
  juicer_metadata = get_juicer_metadata(juicerPath)

  # check if resolution is found!
  if(!any(resolution == juicer_metadata[[2]])){
    stop('resolution is not found.\nAvailable resolutions are ',
         paste(rev(juicer_metadata[[2]]), collapse = ', '))
  }
  validchrom <- !(grepl("all", juicer_metadata[[1]]$chrom, ignore.case = TRUE))
  
  juicer_metadata[[1]] = juicer_metadata[[1]][validchrom,]

  expandedChromosomes = as.data.frame(t(utils::combn(juicer_metadata[[1]]$chrom, m = 2)), stringsAsFactors = F)
  expandedChromosomes = rbind(as.data.frame(cbind(juicer_metadata[[1]]$chrom,juicer_metadata[[1]]$chrom)),
                              expandedChromosomes)
  expandedChromosomes = apply(expandedChromosomes, 2, as.character)
  strawNorm = ifelse(balancing, 'KR', "NONE")
  
  juicerList = lapply(seq_len(nrow(expandedChromosomes)), function(eci){
 
    ec = expandedChromosomes[eci,]
    juicer_in <- tryCatch(
      {
        out <- NULL
        if("matrix" %in% names(as.list( args(strawr::straw) ))){ # check if straw-version has "matrix"-argument
          out <- strawr::straw(matrix = 'observed', 
                        norm = strawNorm,
                        fname =juicerPath,
                        chr1loc =  ec[2],
                        chr2loc =  ec[1],
                        unit = 'BP',
                        binsize = resolution)
        } else {
          out <- strawr::straw(norm = strawNorm,
                        fname =juicerPath,
                        chr1loc =  ec[2],
                        chr2loc =  ec[1],
                        unit = 'BP',
                        binsize = resolution)
        }
        out
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
  rm(juicerList)
  # split into index and signal-files
  SA = splitJuicerData(juicer_data, resolution)
  rm(juicer_data)
  
  SIG = SA[[1]]
  ABS = SA[[2]]
  rm(SA)
  
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
  
  if (n_attributes > 0) {
    block <- 256*4
    found <- 0
    # Find key-value pairs for each attribute by searching 
    # for string terminators (raw == 0)
    while (found < n_attributes * 2) {
      r <- readBin(file2read, "raw", block)
      if (length(w <- head(which(r == 0), 1))) {
        found <- found + 1
        # Rewind to previous terminator
        seek(file2read, -(block - w), origin = "current")
      }
    }
  }
  
  n_chrs <- readBin(file2read, integer(), n = 1, size = 4)
  
  
  chrom_name_length <- lapply(seq_len(n_chrs), function(x){
    chrom <- readBin(file2read, character(), size = 1)
    len <- readBin(file2read, integer(), n = 1)
    data.frame(chrom = chrom, len = len, stringsAsFactors = FALSE)
  })
  chrom_name_length <- do.call(rbind, chrom_name_length)
  # chrom_name_length = data.table::rbindlist(chrom_name_length)
  
  n_BP_Res  = readBin(file2read, integer(), n = 1, size = 4)
  
  
  resolutions = sapply(1:n_BP_Res, function(x){
    readBin(file2read, integer(), n = 1)
  })
  
  close(file2read)
  
  return(list(chrom_name_length, resolutions))
}


splitJuicerData = function(juicer_data, resolution){

  # get all unique bins
  bins_x <- juicer_data[, list(min = min(x), max = max(x)), 
                        by = .(chrom = chrom_x)]
  bins_y <- juicer_data[, list(min = min(y), max = max(y)), 
                        by = .(chrom = chrom_y)]
  
  # Make abs
  abs <- rbind(bins_x, bins_y)
  abs <- abs[, list(min = min(min), max = max(max)), by = "chrom"]
  abs <- abs[, seq(min, max, by = as.integer(resolution)), by = "chrom"]
  setnames(abs, c("V1", "V2"))
  abs[, c("V4") := seq_len(nrow(abs))]
  
  # Make signal
  ## Search for x-index
  SIG <- juicer_data[abs, on = c(chrom_x = "V1", x = "V2")]
  setnames(SIG, 6, "V1")
  ## Search for y-index
  SIG <- SIG[abs, on = c(chrom_y = "V1", y = "V2")]
  setnames(SIG, 7, "V2")
  
  # Clean up signal
  SIG[, c("chrom_x", "x", "chrom_y", "y") := NULL]
  setcolorder(SIG, c(2,3,1))
  setnames(SIG, c("V1", "V2", "V3"))
  
  # Discard NAs (from the join with non-existing bins)
  SIG <- SIG[!is.na(V3)]
  # Ensure upper-triangular format
  SIG[, c("V1", "V2") := list(pmin(V1, V2), pmax(V1, V2))]
  setkeyv(SIG, c("V1", "V2"))

  # Finishing touches on abs
  abs[, V3 := V2 + as.integer(resolution)]
  setcolorder(abs, c("V1", "V2", "V3", "V4"))

  return(list(SIG, abs))
}


