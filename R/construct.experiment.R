#' Construct a HiC-experiment.
#'
#' Make a structure which holds the most needed information of a HiC-experiment.
#'
#' @param signalPath Full path to a HiC-pro-like matrix-file or .sig from juicerToGenova.py
#' @param indicesPath Full path the HiC-pro-like index-file or .bed from juicerToGenova.py
#' @param name The name of the sample.
#' @param color Color associated with sample.
#' @param COMMENTS A place to store some comments.
#' @return A list.
#' @export
construct.experiment <- function(signalPath, indicesPath, name, color = 1, comments = NULL){
  # Check if files exist
  if(!file.exists(signalPath)){stop('Signal file not found.')}
  if(!file.exists(indicesPath)){stop('Index file not found.')}

  ICE <- read.hicpro.matrix(signalPath)
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

    # Named list of centromeric bins: list(chr1 = c(1,2,3,4,etc), "chr2" = c(3,4,5,6,etc))
    CENTROMERES = list()
  )
}
