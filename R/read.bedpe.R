#' Load, parse and check BEDPE.
#'
#' This function loads a BEDPE file. Columns 2 and 4 and 5 and 6 will be sorted numerically.
#' If only a BED-file is present, duplicate BED-columns 1,2,3 as: 1,2,2,1,3,3.
#'
#' @param loop.bed Path to bedpe file.
#' @return A data.frame with the bedpe-values.
read.bedpe <- function(loop.bed, header = T){
  # Load bedpe of loops of interest
  # Chromosome-columns are factors
  LOI <- read.delim(loop.bed, stringsAsFactors=F, h = header)
  # Store extra columns
  LOI.extra <- NA
  if(ncol(LOI) > 6){
    LOI.extra <- LOI[,7:ncol(LOI)]
    LOI <- LOI[,1:6]
  }
  # remove empty rows
  s <- apply(is.na(LOI), 1, sum)
  LOI <- LOI[-which(s > 0),]
  if(is.data.frame(LOI.extra)){
    LOI.extra <- LOI.extra[-which(s > 0),]
  } else {
    LOI.extra <- LOI.extra[-which(s > 0)]
  }
  # Split chromosome
  chromosomeAvector <- as.character(LOI[,1])
  chromosomeBvector <- as.character(LOI[,4])
  # remove whitespaces
  chromosomeAvector <- gsub("[[:space:]]","",chromosomeAvector)
  chromosomeBvector <- gsub("[[:space:]]","",chromosomeBvector)
  # # Check if chromosome entries have 'chr' or add.
  # chromosomeAvector[!grepl(chromosomeAvector,pattern = "^chr")] <- gsub(chromosomeAvector[!grepl(chromosomeAvector,pattern = "^chr")], pattern = '^', replacement = 'chr', perl = T)
  # chromosomeBvector[!grepl(chromosomeBvector,pattern = "^chr")] <- gsub(chromosomeBvector[!grepl(chromosomeBvector,pattern = "^chr")], pattern = '^', replacement = 'chr', perl = T)
  # Sort POS for each BED-entry
  posAdf <- cbind(ifelse(LOI[,2] < LOI[,3], LOI[,2], LOI[,3]), ifelse(LOI[,2] < LOI[,3], LOI[,3], LOI[,2]))
  posBdf <- cbind(ifelse(LOI[,5] < LOI[,6], LOI[,5], LOI[,6]), ifelse(LOI[,5] < LOI[,6], LOI[,6], LOI[,5]))
  # check if posA is indeed upstream of posB
  posAdf.Smaller <- na.omit(posAdf[posAdf[,1] < posBdf[,1],] )
  posBdf.Smaller <- na.omit(posAdf[posAdf[,1] > posBdf[,1],] )
  posAdf.Larger <- na.omit(posAdf[posAdf[,1] > posBdf[,1],] )
  posBdf.Larger <- na.omit(posAdf[posAdf[,1] < posBdf[,1],])
  posAdfsorted <- rbind(posAdf.Larger, posBdf.Larger)
  posBdfsorted <- rbind(posAdf.Smaller, posBdf.Smaller)
  # Bind to new DF
  BEDPE.new <- data.frame(V1 = chromosomeAvector,
                          V2 = posAdfsorted[,1],
                          V3 = posAdfsorted[,2],
                          V4 = chromosomeBvector,
                          V5 = posBdfsorted[,1],
                          V6 = posBdfsorted[,2])
  # Return bedpe
  return(BEDPE.new)
}
