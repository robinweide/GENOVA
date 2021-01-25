
#' Get the percentage cis-contacts genome-wide or from a region
#'
#' @param exp A contacts-objects or a list of them
#' @param bed A bed-dataframe or a (named) list of them.
#' 
#' @section Resolution recommendation: 500kb-1Mb
#'
#' @return A data.table with the sample, percentage and region.
#' @export
cis_trans = function(exp, bed = NULL){
  
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  if (any(c("MAT", "IDX") %in% names(exp))) {
    exp <- list(exp)
  }
  
  if(!is.null(bed)){
    if(!inherits(bed, what = 'list')) {
      bed <- list(bed) 
      names(bed) = 'bed'
    }
    
    if(is.null(names(bed))){
      names(bed) = paste0('bed ',1:length(bed))
    }
    
    bed = lapply(bed, function(x){x[,1:3]})
    bed = data.table::rbindlist(bed, idcol = 'naam')
    bed = lapply(split(bed[,c(2:4,1)],bed$naam), as.data.frame)
  }
  
  
  
  # lapply over exp
  lapOut = lapply(exp, function(x){
    
    if(attr(x, 'scale_cis')){
      stop('The data is cis-scaled. Run load_contacts with scale_cis = F.')
    }
    
    # lapply over bed
    cis = NULL
    chromRange = x$IDX[ , .(first = min(V4)), by = V1]
    chromRange = chromRange[order(chromRange$first),]
    if(!is.null(bed)){
      thisExpOut = lapply(bed, function(bd){
        
        start <- bed2idx(IDX = x$IDX, bd[,1:3], "start")
        end <- bed2idx(IDX = x$IDX, bd[,1:3], "end")
        seqs <- mapply(seq.int, from = start, to = end, SIMPLIFY = FALSE)
        seqs <- SJ(do.call(c, seqs))
        
        left = x$MAT[seqs, on = "V1"]
        right = x$MAT[seqs, on = c(V2 = "V1")]
        bedMAT = rbind(left, right)
        
        F1 = findInterval(bedMAT$V1, chromRange$first)
        F2 = findInterval(bedMAT$V2, chromRange$first)
        
        cisMat = bedMAT[ , sum(V3, na.rm = TRUE), by = ifelse(F1 == F2, T, F) ]
        cis = data.table::data.table(sample = attr(x, 'samplename'), cis = cisMat[cisMat$ifelse == TRUE, 2] / sum(cisMat$V1))
        cis$region = unique(bd$naam)
        colnames(cis)[2] = 'cis'
        cis
      })
      cis = data.table::rbindlist(thisExpOut)
    } else {
      F1 = findInterval(x$MAT$V1, chromRange$first)
      F2 = findInterval(x$MAT$V2, chromRange$first)
      
      cisMat = x$MAT[ , sum(V3), by = ifelse(F1 == F2, T, F) ]
      cis = data.table::data.table(sample = attr(x, 'samplename'), cis = cisMat[cisMat$ifelse == TRUE, 2] / sum(cisMat$V1))
      cis$region = 'genome-wide'
      colnames(cis)[2] = 'cis'
      cis
    }
    
    cis$cis = cis$cis*100
    cis
    
  })
  lapOut = data.table::rbindlist(lapOut)
  
  return(lapOut)
}
