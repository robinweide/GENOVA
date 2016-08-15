#' Load, parse and check BEDPE.
#'
#' This function loads a BEDPE file. Columns 2 and 4 and 5 and 6 will be sorted numerically.
#' If only a BED-file is present, duplicate BED-columns 1,2,3 as: 1,2,2,1,3,3.
#'
#' @param loop.bed Path to bedpe file.
#' @return A data.frame with the bedpe-values.
#' @export
read.bedpe <- function(loop.bed){
  # Load bedpe of loops of interest
  # Chromosome-columns are factors
  LOI <- read.delim(loop.bed, stringsAsFactors=T, h = F)
  # Store extra columns
  LOI.extra <- NA
  if(ncol(LOI) > 6){
    LOI.extra <- LOI[,7:ncol(LOI)]
    LOI <- LOI[,1:6]
  }
  # Split chromosome
  chromosomeAvector <- LOI[,1]
  chromosomeBvector <- LOI[,4]
  # Check if chromosome entries have 'chr' or add.
  if(all(grepl(chromosomeAvector,pattern = "^chr")) == FALSE){
    chromosomeAvector <- gsub(chromosomeAvector, pattern = '^', replacement = 'chr', perl = T)
    chromosomeBvector <- gsub(chromosomeBvector, pattern = '^', replacement = 'chr', perl = T)
  }
  # Sort POS for each BED-entry
  posAdf <- data.frame(t(apply(LOI[,2:3], 1, sort)))
  posBdf <- data.frame(t(apply(LOI[,5:6], 1, sort)))
  # Bind to new DF
  BEDPE.new <- data.frame(V1 = chromosomeAvector,
                          V2 = posAdf[,1],
                          V3 = posAdf[,2],
                          V4 = chromosomeBvector,
                          V5 = posBdf[,1],
                          V6 = posBdf[,2])
  # Add extra columns
  if(!is.na(LOI.extra)){
    BEDPE.new <- cbind(BEDPE.new, LOI.extra)
  }
  # Return bedpe
  return(BEDPE.new)
}
