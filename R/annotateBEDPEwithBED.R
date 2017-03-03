#' Annotate loops
#'
#' Using some genomicRanges-magic, annotate bedpe-files with a bed file.
#'
#' @param bedpe Loop-dataframe with at least six columns
#' @param bed Bed-file, containing a column with annotations
#' @param annotCol Which column-number should we take for the annotation?
#' @param gap Maximal allowed gap between bedpe-bed overlaps. Default = 0.
#' @param verbose Should this function be chatty?
#' @return An annotated bedpe-dataframe
#' @import GenomicRanges
#' @export
annotateBEDPEwithBED <- function(bedpe, bed, annotCol, gap = 0,verbose =F){
  # maybe sort anchors, so upstream is always in bedpe5?

  # check for same chomNames: must not have ^chr
  bedpe[,1] <- gsub(bedpe[,1], pattern = "^chr",perl = T, replacement =  "")
  bedpe[,4] <- gsub(bedpe[,4], pattern = "^chr",perl = T, replacement =  "")
  bed[,1] <- gsub(bed[,1], pattern = "^chr",perl = T, replacement =  "")

  # get 5' anchors
  bedpe5 <- bedpe[,c(1,2,3)]
  colnames(bedpe5) <- c("seqnames", "start", "end")
  bedpe5 <- GenomicRanges::makeGRangesFromDataFrame(bedpe5)
  # get 3' anchors
  bedpe3 <- bedpe[,c(4,5,6)]
  colnames(bedpe3) <- c("seqnames", "start", "end")
  bedpe3 <- GenomicRanges::makeGRangesFromDataFrame(bedpe3)

  # get 5' anchors
  bedAnnot <- bed[,c(1,2,3,annotCol)]
  colnames(bedAnnot) <- c("seqnames", "start", "end","annotation")
  bedAnnot <- GenomicRanges::makeGRangesFromDataFrame(bedAnnot,keep.extra.columns = T)

  # Findoverlaps
  FO5 <- GenomicRanges::findOverlaps(bedpe5, bedAnnot, maxgap = gap)
  FO3 <- GenomicRanges::findOverlaps(bedpe3, bedAnnot, maxgap = gap)
  FO5 <- data.table(as.data.frame(FO5)) ; data.table::setkey(FO5, "queryHits")
  FO3 <- data.table(as.data.frame(FO3)) ; data.table::setkey(FO3, "queryHits")

  # loop over all idx'es and store
  # Output Df has first 6 cols of bedpe, followed bij two columns of 5' and 3' annotations.
  # loops with multiple annnotations will have multiple rows.
  outDF <- NULL
  for(i in 1:nrow(bedpe)){
    if(verbose == TRUE){cat("Annotating all the loops:",i," \r")}
    hits5 <- FO5[list(i)]
    hits3 <- FO3[list(i)]
    N5 <- nrow(hits5)
    N3 <- nrow(hits3)
    if(N5 == 0){
      if(N3 == 0){
        thisRow <- cbind(bedpe[i,], NA ,  NA )
        NC <- ncol(thisRow)
        colnames(thisRow)[c((NC-1),NC)] <- c("5primeAnnotation","3primeAnnotation")
        outDF <- rbind(outDF, thisRow )
      } else {
        for(i3 in 1:N3){
          thisRow <- cbind(bedpe[i,], NA , bed[  unlist(hits3[i3,2] ),4] )
          NC <- ncol(thisRow)
          colnames(thisRow)[c((NC-1),NC)] <- c("5primeAnnotation","3primeAnnotation")
          outDF <- rbind(outDF, thisRow )
        }
      }
    } else { # hits5 is not empty
      if(N3 == 0){
        thisRow <- cbind(bedpe[i,], bed[  unlist(hits5[i5,2]) ,4]  ,  NA )
        NC <- ncol(thisRow)
        colnames(thisRow)[c((NC-1),NC)] <- c("5primeAnnotation","3primeAnnotation")
        outDF <- rbind(outDF, thisRow )
      } else {
        for(i5 in 1:N5){
          for(i3 in 1:N3){

            thisRow <- cbind(bedpe[i,], bed[  unlist(hits5[i5,2]) ,4]  , bed[  unlist(hits3[i3,2]) ,4]  )
            NC <- ncol(thisRow)
            colnames(thisRow)[c((NC-1),NC)] <- c("5primeAnnotation","3primeAnnotation")
            outDF <- rbind(outDF, thisRow )

          }
        }
      }
    }
  }
  return(outDF)
}
