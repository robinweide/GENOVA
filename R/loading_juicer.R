loadJuicer = function(juicerPath, resolution, scale_bp = 1e9, scale_cis = F, balancing = T){

  try_require("strawr", "loadJuicer", "github")
  if (utils::packageVersion("strawr") < "0.0.9") {
    stop("Loading `*.hic` files requires the {strawr} package version 0.0.9+.")
  }
  
  juicerPath <- normalizePath(juicerPath)
  
  # get metadata of juicer-file
  resolutions <- strawr::readHicBpResolutions(juicerPath)
  chromosomes <- strawr::readHicChroms(juicerPath)

  # check if desired resolution is found!
  if(!any(resolution == resolutions)){
    stop('resolution is not found.\nAvailable resolutions are ',
         paste(rev(resolutions), collapse = ', '))
  }
  validchrom <- !(grepl("all", chromosomes$name, ignore.case = TRUE))
  
  chromosomes <- chromosomes[validchrom, , drop = FALSE]
  
  combinations <- matrix(1, nrow = nrow(chromosomes), ncol = nrow(chromosomes))
  combinations[upper.tri(combinations)] <- NA
  combinations <- data.table(
    V1 = chromosomes$name[col(combinations)[!is.na(combinations)]],
    V2 = chromosomes$name[row(combinations)[!is.na(combinations)]]
  )
  
  if (!is.logical(balancing)) {
    balancing <- balancing %in% c("KR", "T", "TRUE")
  }
  
  strawNorm = ifelse(balancing[1], 'KR', "NONE")
  
  juicerList = lapply(seq_len(nrow(combinations)), function(eci){
 
    ec = unlist(unclass(combinations[eci,]))
    juicer_in <- tryCatch(
      {
        out <- NULL
        if("matrix" %in% names(formals(strawr::straw))){ # check if straw-version has "matrix"-argument
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
        warning(c(
          paste0("Failed to retrieve combination of ", ec[1], 
                 " and ", ec[2], ".\n"),
          cond$message
        ), call. = FALSE)
        return(NULL)
      })
    
    if (is.null(juicer_in)) {
      return(NULL)
    }
    
    # remove NaN
    juicer_in = juicer_in[!is.nan(juicer_in[,3]), , drop = FALSE]
    if (nrow(juicer_in) == 0) {
      juicer_in = NULL
    }
    
    # decorate
    juicer_in$chrom_x = ec[1]
    juicer_in$chrom_y = ec[2]
    juicer_in = juicer_in[, c(4,1,5,2,3)]
    juicer_in
  })

  juicer_data = data.table::rbindlist(juicerList)
  if (nrow(juicer_data) < 1) {
    stop(paste0("Failed to read file: ", juicerPath))
  }
  
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


