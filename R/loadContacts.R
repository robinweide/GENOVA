#' Construct a contacts-object.
#'
#' Make an object which holds the most needed information of a HiC-experiment at a given resolution.
#'
#' @param signalPath Full path to a HiC-pro-like matrix-file
#' @param indicesPath Full path the HiC-pro-like index-file
#' @param coolerPath Full path the .cooler-file
#' @param juicerPath Full path the .hic-file
#' @param juicerResolution Set the desired resolution of the matrix when using juicer-data.
#' @param sampleName The name of the sample.
#' @param centromeres A data.frame with three columns per chromosome: chromosome name, start-position and end-position of the centromeric region.
#' @param color Color associated with sample.
#' @param Znorm Normalise the matrices per-chromosome with a Z-score.
#' @param balancing TRUE (default) will perform matrix balancing for .cooler and KR for.hic.
#' @param BPscaling Scale contacts to have agenome-wide sum of `BPscaling` reads (default: 1e9). Set to NULL to skip this.
#' @param verbose Do you want updates during the construction?
#' @note
#' Some reference genomes have very small "random" or "patch" chromosomes,
#' which can have zero contacts mapped to it (at certain resolutions).
#' `loadContacts` checks this and omits these chromosomes in the resulting
#' experiment-object. The RMCHROM-flag will also be set to TRUE: this will
#' help other GENOVA-functions to deal better with this problem. There is a
#' slight performance-cost during the construction of the object, however.
#' @examples
#' WT_10kb_hicpro <- loadContacts(signalPath = "WT_10kb_iced.matrix", indicesPath = "WT_10kb_abs.bed", sampleName = "WT", color = "black")
#' WT_10kb_cooler <- loadContacts(coolerPath = "WT_10kb.cooler", balancing =T, sampleName = "WT", color = "black")
#' WT_10kb_juicer <- loadContacts(juicerPath = "WT_10kb_iced.hic", juicerResolution = 10e3, balancing =T, sampleName = "WT", color = "black")
#' @return An contacts-object, which is a named list of contacts, indices and attributes for a Hi-C matrix of a given sample at a given resolution.
#' @export
loadContacts = function(signalPath = NULL, indicesPath = NULL,
                        coolerPath = NULL,
                        juicerPath = NULL,
                        juicerResolution = 10e3,
                        sampleName,
                        centromeres = NULL,
                        color = NULL,
                        Znorm = F,
                        bpscaling = 1e9,
                        balancing = T,
                        verbose = T){

  doJuicer = F
  doCooler = F
  doHiCpro = F
  balanced = F
  ##############################################################################
  #################################################### what are we dealing with?
  ##############################################################################
  if(!is.null(juicerPath)){
    ##################################################################### juicer

    # a juicerpath is given. Does it exist?
    if(!file.exists(juicerPath)){
      stop("juicerPath doesn't point to an existing .hic-file.")
    }

    # is strawr installed?
    try_require('strawr', "loadContacts", source = 'github')

    doJuicer = T

    # did you also load in hicpro?
    if(!is.null(signalPath)){
      warning('Both hicpro- and juicebox-files provided. Will use juicebox.')
    }

    if(!is.null(coolerPath)){
      warning('Both cooler- and juicebox-files provided. Will use juicebox.')
    }

  } else if(!is.null(signalPath)){
    ##################################################################### hicpro

    # a signalPath is given. Does it exist?
    if(!file.exists(signalPath)){
      stop("signalPath doesn't point to an existing file.")
    }

    if(!file.exists(indicesPath)){
      stop("indicesPath doesn't point to an existing file.")
    }

    doHiCpro = T

    # did you also load in cooler?
    if(!is.null(coolerPath)){
      warning('Both cooler- and hicpro-files provided. Will use hicpro.')
    }

  } else if(!is.null(coolerPath)){
    ##################################################################### cooler

    # a coolerPath is given. Does it exist?
    if(!file.exists(coolerPath)){
      stop("coolerPath doesn't point to an existing file.")
    }

    try_require('rhdf5', 'loadCooler', source = 'BIOC')

    doCooler = T

  } else {
    stop('please provide either signal/indicesPath, a coolerPath or a juicerPath.')
  }

  ##############################################################################
  #################################################################### load data
  ##############################################################################
  if(verbose){message('Reading data...')}
  if(doJuicer){
    sig_ind = loadJuicer(juicerPath, juicerResolution, norm = bpscaling, balancing = balancing)
    balanced = balancing
  }

  if(doHiCpro){
    sig_ind = loadHiCpro(signalPath, indicesPath, norm = bpscaling)
    balanced = any(data.table::fread(signalPath, nrows = 1000)[,3] %% 1 != 0)
  }

  if(doCooler){
    sig_ind = loadCooler(coolerPath, norm = bpscaling,balancing = balancing)
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
  res <- as.numeric(median(index$V3 - index$V2))
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
  RMCHROM = length(foundChroms) == length(chromRLE$values)

  ##############################################################################
  ############################################################ check centromeres
  ##############################################################################
  if (!is.null(centromeres)) {
    centromeres <- clean_centromeres(centromeres, res)
  }

  ##############################################################################
  ###################################################################### Z-score
  ##############################################################################
  if(Znorm){
    if(verbose){message('Running Z-score normalisation...')}
    signal = zscore_hic(signal, index)
  }



  ##############################################################################
  ###############################################################  Contruct list
  ##############################################################################
  color = if(is.null(color)){color = 'black'}
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
  znorm = Znorm, samplename = sampleName, color = color,
  resolution = res, rmChrom = RMCHROM, balanced = balanced,
  # packgs
  package = "GENOVA")

  return(out)
}



loadHiCpro = function(signalPath, indicesPath, norm){

  SIG <- read.hicpro.matrix(signalPath, norm = norm)
  ABS <- data.table::fread(indicesPath, header = F, data.table = F)

  return(list(SIG, ABS))
}


zscore_hic = function(SIG, ABS){

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

clean_centromeres <- function(centros, resolution) {
  # Essentially does the same as `reduce(a_granges_object, min.gapwidth = resolution)

  centros <- as.data.frame(centros)

  # Drop unused factor levels
  if (is.factor(centros[, 1])) {
    centros[, 1] <- droplevels(centros[, 1])
  }

  # Set column names and order on chromosome and start
  centros <- setNames(centros, paste0("V", 1:3))
  centros <- centros[order(centros[, 1], centros[, 2]), ]

  # Test wether adjacent entries should be merged
  centros$merge <- c(FALSE, vapply(seq_len(nrow(centros) - 1), function(i) {

    # Is the difference between end and next start sub-resolution?
    test <- abs(centros[i, 3] - centros[i + 1, 2]) < resolution

    # Are the entries on the same chromosome?
    test && centros[i, 1] == centros[i + 1, 1]
  }, logical(1)))

  # Merge the TRUE entries
  while (any(centros$merge)) {
    # What is the next entry to be merged?
    i <- which(centros$merge)[1]

    # Replace end of previous entry with end of current entry
    centros[i - 1, 3] <- centros[i, 3]

    # Remove current entry
    centros <- centros[-i, ]
  }

  # Return centromeres
  centros[, 1:3]
}
