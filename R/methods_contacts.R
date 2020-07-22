# Subsetting --------------------------------------------------------------

#' @title subset
#' @param x A \code{contacts}-object
#' @param chrom A vector of chromosomes to subset
#' @param start A single postion in bp.
#' @param end A single posttion in bp.
#' @param ... Currently not in use.
#' @examples 
#' \dontrun{
#' WT_chr1 <- subset(WT_new, chrom = 'chr1')
#' WT_acrocentric <- subset(WT_new, chrom = c('chr13',
#'                                            'chr14',
#'                                            'chr15',
#'                                            'chr21',
#'                                            'chr22'))
#' }
#' @rdname subset
#' @export
subset.contacts <- function(x, chrom = NULL, start = 0, end = Inf, ...){
  
  if(is.null(chrom)){
    
    if(start == 0 & end == Inf){
      return(x)
    }
    chrom = x$CHRS
  }
  
  idx_chrom <- x$IDX[chrom]
  idx <- idx_chrom[V2 >= start & V3 <= end]
  
  mat <- x$MAT[V1 %in% idx$V4]
  mat <- mat[V2 %in% idx$V4]
  
  
  
  out = structure(list(
    # Iced HiC-matrix in three-column format (i.e. from HiC-pro)
    MAT = mat,
    
    # HiC-index in four-column format (i.e. from HiC-pro)
    IDX = idx,
    
    # Available chromosomes
    CHRS = chrom,
    
    # DF with chr, start, stop
    CENTROMERES = x$CENTROMERES[x$CENTROMERES[,1] %in% chrom,]
    
  ))
  
  attributes(out) <- attributes(x)
  
  out
}


is_contacts <- function(x) {
  inherits(x, "contacts")
}
