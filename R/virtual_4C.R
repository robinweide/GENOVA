#' Get an average virtual 4C of 2D-regions (e.g. loops).
#'
#' Extract matrices around a defined region a.k.a. the viewpoint
#'
#' @param explist Either a single GENOVA \code{contacts} object or list of
#'   GENOVA \code{contacts}.
#' @param viewpoint The viewpoint: will take the middle Hi-C bin if it spans
#'   mulitple bins.
#' @param xlim A vector of two with the flanking basepairs up- and downstream of
#'   the viewpoint, resp.
#' @return A virtual4C_discovery object.
#' @export
virtual_4C <- function(explist, viewpoint, xlim = NULL){
  # ! someday: allow mulitple samples
  explist  <- check_compat_exp(explist)
  expnames <- if (is.null(names(explist))) {
    vapply(explist, attr, character(1L), "samplename")
  } else {
    names(explist)
  }
  IDX <- explist[[1]]$IDX
  
  vp_idx <- median(bed2idx(IDX, viewpoint))
  
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  signal <- NULL
  if( is.null(xlim) ){
    # run genome-wide ==========================================================
    signal <- lapply(seq_along(explist), function(i) {
      explist[[i]]$MAT[V1 == vp_idx | V2 == vp_idx, 
                       list(V1, V2, V3, exp = i)]
    })
    signal <- rbindlist(signal)
    upstream_signal <- signal[V2 %in% vp_idx][, c(1, 3, 4)]
    downstream_signal <- signal[V1 %in% vp_idx][, 2:4]
  } else {
    # run for a region =========================================================
    flank <- floor(xlim/attr(explist[[1]], 'resolution'))

    range_idx <- unlist(vp_idx - flank[1]):unlist(vp_idx + flank[2])
    
    upstream_signal <- lapply(seq_along(explist), function(i) {
      explist[[i]]$MAT[list(range_idx, vp_idx), 
                       nomatch = 0][, list(V1, V3, exp = i)]
    })
    upstream_signal <- rbindlist(upstream_signal)
    downstream_signal <- lapply(seq_along(explist), function(i) {
      explist[[i]]$MAT[list(vp_idx, range_idx),
                       nomatch = 0][, list(V2, V3, exp = i)]
    })
    downstream_signal <- rbindlist(downstream_signal)
  }
  
  signal <- rbind(upstream_signal, downstream_signal, use.names = FALSE)
  colnames(signal) <- c('idx', 'signal', 'experiment')
  signal <- signal[IDX, on = "idx==V4", nomatch = 0]
  
  # set basepairs --------------------------------------------------------------
  signal$mid = rowMeans(signal[, 5:6])
  
  if( !is.null(xlim) ){
    signal <- signal[V1 == viewpoint[1,1]]
  }
 
  
  # output ---------------------------------------------------------------------
  signal <- unique(signal[, c(4, 7, 2, 3)])
  colnames(signal) <- c("chromosome", "mid", "signal", "experiment")
  signal$experiment <- expnames[signal$experiment]
  signal <- list(data = signal)
  
  viewpoint <- data.frame(chrom = viewpoint[1, 1],
                          start = viewpoint[1, 2],
                          end   = viewpoint[1, 3],
                          exp = expnames)
  
  signal <- structure(signal, 
                      class = "virtual4C_discovery",
                      'viewpoint' = viewpoint, 
                      'xlim' = xlim,
                      'sample' = expnames,
                      'resolution' = attr(explist[[1]], 'resolution'),
                      package = "GENOVA")
  
  signal
}


