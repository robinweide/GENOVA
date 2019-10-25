#' Get an average virtual 4C of 2D-regions (e.g. loops).
#'
#' Extract matrices around a defined region a.k.a. the viewpoint
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param viewpoint The viewpoint: will take the middle Hi-C bin if it spans mulitple bins.
#' @param xlim A vector of two with the flanking basepairs up- and downstream 
#' of the viewpoint, resp. 
#' @return A virtual4C_discovery object.
#' @export
virtual_4C <- function(exp, viewpoint, xlim = NULL){
  # ! someday: allow mulitple samples
  vp_idx <- median(bed2idx(exp$IDX, viewpoint))
  
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  signal <- NULL
  if( is.null(xlim) ){
    # run genome-wide ==========================================================
    signal <- exp$MAT[V1 == vp_idx | V2 == vp_idx]
    
    upstream_signal   <- signal[V2 %in% vp_idx][,c(1,3)]
    downstream_signal <- signal[V1 %in% vp_idx][,2:3]
  } else {
    # run for a region =========================================================
    flank <- floor(xlim/attr(exp, 'resolution'))

    range_idx <- unlist(vp_idx - flank[1]):unlist(vp_idx + flank[2])
    
    upstream_signal <- exp$MAT[list(range_idx,vp_idx), nomatch = 0][,c(1,3)]
    downstream_signal <- exp$MAT[list(vp_idx, range_idx), nomatch = 0][,2:3]
  }
  
  signal <- rbind(upstream_signal, downstream_signal, use.names = FALSE)
  colnames(signal) = c('V4','signal')
  setkey(signal, 'V4')
  
  # set basepairs --------------------------------------------------------------
  bed_ranges <- exp$IDX[match(signal$V4, exp$IDX$V4)]
  setkey(bed_ranges, 'V4')
  signal <- signal[bed_ranges]
  signal$mid = rowMeans(signal[,4:5])
  
  if( !is.null(xlim) ){
    signal <- signal[V1 == viewpoint[1,1]]
  }
 
  
  # output ---------------------------------------------------------------------
  signal <- unique(signal[,c(3,6,2)])
  colnames(signal) <- c('chromosome','mid','signal')
  signal <- list(data = signal)
  
  signal <- structure(signal, 
                      class = "virtual4C_discovery",
                      'viewpoint' = viewpoint, 
                      'xlim' = xlim,
                      'sample' = attr(exp, 'sample'),
                      'resolution' = attr(exp, 'resolution'),
                      package = "GENOVA")
  
  signal
}


