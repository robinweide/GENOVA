#' Construct a HiC-experiment.
#'
#' Make a structure which holds the most needed information of a HiC-experiment.
#'
#' @param ICEDpath Full path to HiC-pro matrix-file.
#' @param BEDpath Full path the HiC-pro index-file
#' @param SAMPLENAME The name of the sample.
#' @param COLOR Color associated with sample.
#' @param COMMENTS A place to store some comments.
#' @return A list.
#' @export
construct.experiment <- function(ICEDpath,BEDpath, SAMPLENAME, COLOR = 1, COMMENTS = NULL){
  # Check if files exist
  if(!file.exists(ICEDpath)){stop('ICE-matrix file not found.')}
  if(!file.exists(BEDpath)){stop('ICE-index file not found.')}
  
  ICE <- read.hicpro.matrix(ICEDpath)
  ABS <- data.table::fread(BEDpath, header = F, data.table = F)
  
  # check for similar binsizes BED and ICE
  ICEmaxID <- max(max(ICE$V1), max(ICE$V2))
  ABSmaxID <- max(ABS[,4])
  if( ICEmaxID > ABSmaxID){
    stop("Undefined HiC-bins.\nWe checked if the highest Hi-C bin in the Matrix-file could be found in the BED-file: it couldn't.\nPlease check if BED and Matrix are from the same reference genome and resolution!")
  }
  if( abs(ICEmaxID - ABSmaxID) > ICEmaxID*0.1 ) {
    warning("Undefined BED-bins.\nWe checked if the highest Hi-C bin in the BED-file could be found in the Matrix-file: it couldn't.\nPlease check if BED and Matrix are from the same reference genome and resolution!")
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
    NAME = SAMPLENAME,
    
    # Resolution of sample
    RES = res,
    
    # Available chromosomes
    CHRS = chromVector,
    
    # Color of sample (optional, but recommended for running RCP)
    COL = COLOR,
    
    # Comments
    COMM = COMMENTS
    
  )
}
