#' Construct a contacts-object.
#'
#' Make an object which holds the most needed information of a HiC-experiment at a given resolution.
#'
#' @param signal_path Full path to an Hi-C file. This could be \*.matrix (for HiCpro), \*.cooler or \*.hic (for juicer).
#' @param indices_path Full path the HiC-pro-like index-file. Required when `signal_path` is \*.matrix.
#' @param resolution Set the desired resolution of the matrix when using juicer-data.
#' @param sample_name The name of the sample.
#' @param centromeres A data.frame with three columns per chromosome: chromosome name, start-position and end-position of the centromeric region.
#' @param colour colour associated with sample.
#' @param z_norm Normalise the matrices per-chromosome with a Z-score.
#' @param balancing TRUE (default) will perform matrix balancing for .cooler and KR for.hic.
#' @param scale_bp Scale contacts to have agenome-wide sum of `scale_bp` reads (default: 1e9). Set to NULL to skip this.
#' @param scale_cis Only scale with cis-contacts.
#' @param legacy Get a pre-v1 object (mimics the output of construct.experiment.)
#' @param verbose Do you want updates during the construction?
#' @note
#' Some reference genomes have very small "random" or "patch" chromosomes,
#' which can have zero contacts mapped to it (at certain resolutions).
#' `load_contacts` checks this and omits these chromosomes in the resulting
#' experiment-object. The RMCHROM-flag will also be set to TRUE: this will
#' help other GENOVA-functions to deal better with this problem. There is a
#' slight performance-cost during the construction of the object, however.
#' @examples
#' \dontrun{
#' WT_10kb_hicpro <- load_contacts(signal_path = "WT_10kb_iced.matrix", 
#'                                 indices_path = "WT_10kb_abs.bed", 
#'                                 sample_name = "WT", 
#'                                 colour = "black")
#' WT_10kb_cooler <- load_contacts("WT_10kb.cooler", 
#'                                 balancing =T, 
#'                                 sample_name = "WT", 
#'                                 colour = "black")
#' WT_10kb_juicer <- load_contacts("WT_10kb_iced.hic", 
#'                                 resolution = 10e3, 
#'                                 balancing =T, 
#'                                 sample_name = "WT", 
#'                                 colour = "black")
#' }
#' @return An contacts-object, which is a named list of contacts, indices and 
#' attributes for a Hi-C matrix of a given sample at a given resolution.
#' @export
load_contacts = function(signal_path, 
                         indices_path = NULL,
                         resolution = 10e3,
                         sample_name = NULL,
                         centromeres = NULL,
                         colour = NULL,
                         z_norm = FALSE,
                         scale_bp = 1e9,
                         scale_cis = FALSE,
                         balancing = TRUE,
                         legacy = FALSE,
                         verbose = TRUE){
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  if(is.null(sample_name)){
    
    sample_name = basename(signal_path)
    
  }
  doJuicer = F
  doCooler = F
  doHiCpro = F
  balanced = F
  juicerPath = NULL
  coolerPath = NULL
  ##############################################################################
  #################################################### what are we dealing with?
  ##############################################################################
  inputType = switch(tools::file_ext(tolower(signal_path)), 
                     'matrix' = 'hicpro',
                     'cooler' = 'cooler',
                     'cool' = 'cooler',
                     'mcool' = 'cooler',
                     'hic'    = 'juicer')
  
  if(inputType == 'juicer'){juicerPath = signal_path}
  if(inputType == 'cooler'){coolerPath = signal_path}
  
  if(!is.null(juicerPath)){
    ##################################################################### juicer
    
    
    
    if(grepl(juicerPath ,pattern = '^http')){
      stop('signal_path starts with http, but the current version does not support it.')
      # signal_path points to an URL!
      if(unname(RCurl::url.exists(signal_path, .header = T)[7] == "404")){
        stop('signal_path starts with http, but url is not found (404).')
      }
      
    } else if(!file.exists(juicerPath)){
      stop("juicerPath doesn't point to an existing .hic-file.")
    }
    
    # is strawr installed?
    try_require('strawr', "load_contacts", source = 'github')
    
    doJuicer = T
    
    
  } else if(!is.null(coolerPath)){
    ##################################################################### cooler
    
    # a coolerPath is given. Does it exist?
    if(!file.exists(coolerPath)){
      stop("coolerPath doesn't point to an existing file.")
    }
    
    try_require('rhdf5', 'loadCooler', source = 'BIOC')
    
    doCooler = T
    
  } else if(!is.null(signal_path)){
    ##################################################################### hicpro
    
    # a signal_path is given. Does it exist?
    if(!file.exists(signal_path)){
      stop("signal_path doesn't point to an existing file.")
    }
    
    if(!file.exists(indices_path)){
      stop("indices_path doesn't point to an existing file.")
    }
    
    doHiCpro = T
    
  } else {
    stop('please provide either .matrix/indices_path, a .cooler or a .hic.')
  }
  
  ##############################################################################
  #################################################################### load data
  ##############################################################################
  if(verbose){message('Reading data...')}
  if(doJuicer){
    sig_ind = loadJuicer(juicerPath, resolution, scale_bp = scale_bp, scale_cis = scale_cis, balancing = balancing)
    balanced = balancing
  }
  
  if(doHiCpro){
    sig_ind = loadHiCpro(signal_path, indices_path, scale_bp = scale_bp, scale_cis = scale_cis )
    balanced = any(data.table::fread(signal_path, nrows = 1000)[,3] %% 1 != 0)
  }
  
  if(doCooler){
    sig_ind = loadCooler(coolerPath, scale_bp = scale_bp,balancing = balancing, scale_cis = scale_cis, resolution = resolution)
    balanced = balancing
  }
  
  ##############################################################################
  ########################################################## check index and sig
  ##############################################################################
  signal = sig_ind[[1]]
  if(!data.table::is.data.table(signal)){
    stop('signal is no dt')
  }
  if(is.null(data.table::key(signal))){
    data.table::setkey(signal, "V1", "V2")
  }
  
  index = sig_ind[[2]]
  if(!data.table::is.data.table(index)){
    stop('index is no dt')
  }
  if(is.null(data.table::key(index))){
    data.table::setkey(index, "V1", "V2")
  }
  res <- as.numeric(stats::median(index$V3 - index$V2))
  ##############################################################################
  ################################################################### check bins
  ##############################################################################
  
  ICEmaxID <- max(max(signal$V1), max(signal$V2))
  ABSmaxID <- max(index[, 4])
  if (ICEmaxID > ABSmaxID) {
    stop("Undefined HiC-bins.\nWe checked if the highest Hi-C bin in the signal-file could be found in the indices-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  if (abs(ICEmaxID - ABSmaxID) > ICEmaxID * 0.1) {
    warning("Undefined BED-bins.\nWe checked if the highest Hi-C bin in the indices-file could be found in the signal-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  
  ##############################################################################
  ########################################################## remove empty chroms
  ##############################################################################
  chromRLE = rle(index$V1)
  CS = cumsum(c(1,chromRLE$lengths))
  
  tmpsignal = signal
  tmpsignal$C1 = chromRLE$values[findInterval(tmpsignal$V1, CS)]
  tmpsignal$C2 = chromRLE$values[findInterval(tmpsignal$V2, CS)]
  foundChroms = unique(c(unique(tmpsignal$C1), unique(tmpsignal$C2)))
  RMCHROM = length(foundChroms) != length(chromRLE$values)
  
  ##############################################################################
  ###################################################################### Z-score
  ##############################################################################
  if(z_norm){
    if(verbose){message('Running Z-score normalisation...')}
    signal = zscore_hic(signal, index)
  }
  
  ##############################################################################
  ############################################################# set to upper tri
  ##############################################################################
  if(!all(signal$V1 <= signal$V2)){
    signal[, c("V1", "V2") := list(pmin(V1, V2), pmax(V1, V2))]
    setkeyv(signal, c("V1", "V2"))
  }
  
  ##############################################################################
  ############################################################ check centromeres
  ##############################################################################
  
  if (!is.null(centromeres)) {
    centromeres <- clean_centromeres(centromeres, index)
  } else {
    centromeres <- find_centromeres(signal, index)
  }
  
  ##############################################################################
  ###############################################################  Contruct list
  ##############################################################################
  if( is.null(colour) ){
    colour = 'black'
  }
  out = structure(list(
    # Iced HiC-matrix in three-column format (i.e. from HiC-pro)
    MAT = signal,
    
    # HiC-index in four-column format (i.e. from HiC-pro)
    IDX = index,
    
    # Available chromosomes
    CHRS = foundChroms,
    
    # DF with chr, start, stop
    CENTROMERES = centromeres
    
  ),
  # classes
  class = "contacts",
  # attrs
  znorm = z_norm, samplename = sample_name, colour = colour,
  resolution = res, rmChrom = RMCHROM, balanced = balanced, scale_cis = scale_cis,
  # packgs
  package = "GENOVA")
  
  
  if(legacy){
    warning("Legacy is a soft deprecation option, which ~someday~ will be removed.\nUsing the legacy option means you're longing for the good old days.
            Don't: the future has zoomy functions and shiny plots.")
    out = legacy_out(out)
  }
  
  return(out)
}



loadHiCpro = function(signal_path, indices_path, scale_bp, scale_cis){
  
  # init
  .          <- NULL
  V1         <- NULL
  V4         <- NULL
  
  ABS <- data.table::fread(indices_path, header = F, data.table = T)
  
  SIG = NULL
  if(scale_cis){
    SIG <- read.hicpro.matrix(signal_path, scale_bp = NULL)
    
    chromRange = ABS[ , .(first = min(V4)), by = V1]
    chromRange = chromRange[order(chromRange$first),]
    
    F1 = findInterval(SIG$V1, chromRange$first)
    F2 = findInterval(SIG$V2, chromRange$first)
    
    SIG$V3 <- scale_bp * SIG$V3 / sum(SIG[ifelse(F1 == F2, T, F), 3])
    
  } else {
    SIG <- read.hicpro.matrix(signal_path, scale_bp = scale_bp)
  }
  
  return(list(SIG, ABS))
}


zscore_hic = function(SIG, ABS){
  
  # init
  D          <- NULL
  V3         <- NULL
  C1         <- NULL
  C2         <- NULL
  
  chromRLE = rle(ABS$V1)
  
  CS = cumsum(c(1,chromRLE$lengths))
  SIG$C1 = chromRLE$values[findInterval(SIG$V1, CS)]
  SIG$C2 = chromRLE$values[findInterval(SIG$V2, CS)]
  SIG$D = abs(SIG$V1-SIG$V2)
  
  SIG[SIG$C1 != SIG$C2, 'D'] = 0
  
  Z = SIG[ , V3 := scale(V3), by = list(C1, C2, D)]
  Z = Z[,-c(4:6)]
  Z$V3[is.nan(Z$V3)] = 1
  
  data.table::setkey(Z, "V1", "V2")
  
  return(Z)
}

legacy_out = function(out){
  newOut = list(
    # Iced HiC-matrix in three-column format (i.e. from HiC-pro)
    ICE = out$MAT,
    
    # HiC-index in four-column format (i.e. from HiC-pro)
    ABS = data.table::as.data.table(out$IDX),
    
    # Name of sample
    NAME = attr(out, 'samplename'),
    
    # Resolution of sample
    RES = attr(out, 'resolution'),
    
    # Available chromosomes
    CHRS = out$CHRS,
    
    # colour of sample (optional, but recommended for running RCP)
    COL = attr(out, 'colour'),
    
    # comments
    COMM = "",
    
    # Vector of masked bins
    MASK = vector(),
    
    # DF with chr, start, stop
    CENTROMERES = out$CENTROMERES,
    
    # Are there now missing bins?
    RMCHROM =  attr(out, 'rmChrom')
    
  )
  return(newOut)
  
}

#' Read hicpro three column matrix format.
#'
#' This function loads a HiC-pro file as a matrix. It assumes a three-column
#' layout: bin1, bin2 and score. All scors are normalised to contacts per *scale_bp*
#' total contacts.
#'
#' @param file Full path to file.
#' @param scale_bp Normalising factor. Set to NULL to skip scale_bp.
#' @return A data.table with normalised counts.
read.hicpro.matrix <- function(file, scale_bp = 1e9) {
  data <- data.table::fread(file)
  data.table::setkey(data, "V1", "V2")
  if (!is.null(scale_bp)) {
    data$V3 <- scale_bp * data$V3 / sum(data$V3)
  }
  return(data)
}
