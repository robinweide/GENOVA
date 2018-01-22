#' Construct a HiC-experiment.
#'
#' Make a structure which holds the most needed information of a HiC-experiment.
#'
#' @param signalPath Full path to a HiC-pro-like matrix-file or .sig from juicerToGenova.py
#' @param indicesPath Full path the HiC-pro-like index-file or .bed from juicerToGenova.py
#' @param name The name of the sample.
#' @param ignore.checks Skip the checks for empty matrices: EXPERT-ONLY!
#' @param color Color associated with sample.
#' @param COMMENTS A place to store some comments.
#' @param centromeres A data.frame with three columns per chromosome: chromosome name, start-position and end-position of the centromeric region.
#' @param BPscaling Scale contacts to havebgenome-wide sum of [BPscaling] reads.
#' @return A list.
#' @export
construct.experiment <- function(signalPath, indicesPath, name, ignore.checks = F, centromeres = NULL,  color = 1, comments = NULL,BPscaling = 1e9){
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
    for(C in chromVector){
      CIDX = ABS[ABS[,1] == C, 4]
      x = any( ICE$V1 %in%CIDX)
      y = any( ICE$V2 %in% CIDX )

      if(all(x,y)){

      } else {
        # no ICE-data of chrom. warn and delete!
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
