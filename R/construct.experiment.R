#' Construct a HiC-experiment.
#'
#' Make an object which holds the most needed information of a HiC-experiment at a given resolution.
#'
#' @param signalPath Full path to a HiC-pro-like matrix-file or .sig from juicerToGenova.py
#' @param indicesPath Full path the HiC-pro-like index-file or .bed from juicerToGenova.py
#' @param name The name of the sample.
#' @param ignore.checks Skip the checks for empty matrices: EXPERT-ONLY!
#' @param color Color associated with sample.
#' @param comments A place to store some comments.
#' @param centromeres A data.frame with three columns per chromosome: chromosome name, start-position and end-position of the centromeric region.
#' @param BPscaling Scale contacts to have agenome-wide sum of [BPscaling] reads (default: 1000M). Set to NULL to skip this.
#' @note
#' Some reference genomes have very small "random" or "patch" chromosomes, which can have zero contacts mapped to it (at certain resolutions). Construct.experiment checks this and omits these chromosomes in the resulting experiment-object. The RMCHROM-flag will also be set to TRUE: this will help other GENOVA-functions to deal better with this problem. There is a slight performance-cost during the construction of the experiment, however. Therefore, experienced users can set ignore.checks to TRUE to skip all of this, keeping in mind that some functions will not work properly/at all is there are these zero-coverage chromosomes in their data.
#' @examples
#' WT_10kb <- construct.experiment(ignore.checks = T, signalPath = 'WT_10kb_iced.matrix', indicesPath = 'WT_10kb_abs.bed', name = "WT", color = "black")
#' @return An experiment-object, which is a named list of contacts, indices and metadata for a Hi-C matrix of a given sample at a given resolution.
#' @export
construct.experiment <- function(signalPath, indicesPath, name, ignore.checks = F, centromeres = NULL,  color = 1, comments = NULL, BPscaling = 1e9){
  # Check if files exist
  if(!file.exists(signalPath)){stop('Signal file not found.')}
  if(!file.exists(indicesPath)){stop('Index file not found.')}

  ICE <- read.hicpro.matrix(signalPath, norm=BPscaling)
  ABS <- data.table::fread(indicesPath, header = F, data.table = F)

  # check for similar binsizes BED and ICE
  ICEmaxID <- max(max(ICE$V1), max(ICE$V2))
  ABSmaxID <- max(ABS[,4])
  if( ICEmaxID > ABSmaxID){
    stop("Undefined HiC-bins.\nWe checked if the highest Hi-C bin in the signal-file could be found in the indices-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  if( abs(ICEmaxID - ABSmaxID) > ICEmaxID*0.1 ) {
    warning("Undefined BED-bins.\nWe checked if the highest Hi-C bin in the indices-file could be found in the signal-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  chromVector <- as.character(unique(ABS$V1))
  res = as.numeric( median(ABS$V3-ABS$V2)  )

  if(!is.null(centromeres)){
    colnames(centromeres)[1:3] <- paste0('V',1:3)
  }

  # check is all chromosomes have actual data
  RMCHROM = F
  if(!ignore.checks){
    IIDX = unique(c(ICE$V1, ICE$V2))
    for(C in chromVector){
      CIDX = ABS[ABS[,1] == C, 4]

      if(!any(IIDX %in% CIDX)   ){
        warning(paste0("No contacts of ", C))
        ABS = ABS[ !ABS[,1] == C,]
      }

    }

  }
  chromVector.bk = chromVector
  chromVector <- as.character(unique(ABS$V1))

  if(!length(chromVector) == length(chromVector.bk)){
    RMCHROM = T
    warning("Some chromosomes have been removed, since these had no data in the signal-matrix.\nTo turn this off, set ignore.checks to TRUE.")
  }
  # Contruct list
  list(
    # Iced HiC-matrix in three-column format (i.e. from HiC-pro)
    ICE = ICE,

    # HiC-index in four-column format (i.e. from HiC-pro)
    ABS = ABS,

    # Name of sample
    NAME = name,

    # Resolution of sample
    RES = res,

    # Available chromosomes
    CHRS = chromVector,

    # color of sample (optional, but recommended for running RCP)
    COL = color,

    # comments
    COMM = comments,

    # Vector of masked bins
    MASK = vector(),

    # DF with chr, start, stop
    CENTROMERES = centromeres,

    # Are there now missing bins?
    RMCHROM = RMCHROM

  )
}
