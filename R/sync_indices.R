#' Synchronising indices across experiments.
#'
#' Particularly when Hi-C data is loaded from different sources, it occurs that
#' the assigned bin indices don't match between datasets. This function
#' re-indexes a list of experiments to a shared set of bin indices.
#'
#' @param explist A \code{list} of GENOVA \code{contacts} objects of the same
#'   resolution.
#'
#' @details Chromosome naming conventions are expected to be equal among input,
#'   e.g. all \code{contacts} objects use \code{"chr1"}, \code{"chr2"} etc. or
#'   all \code{contacts} objects use \code{"1"}, \code{"2"} etc.
#'
#'   The first \code{contacts} object is used as a template. This should only
#'   matter when the different experiments have slightly different end-positions
#'   at e.g. chromosome ends, in which case the ends of the first
#'   \code{contacts} object is used..
#'
#' @return A \code{list} of GENOVA \code{contacts} objects
#' @export
#'
#' @examples
#' \dontrun{
#' # Data loaded from Hi-C Pro
#' exp1 <- load_contacts(signal_path = "exp1_10kb_iced.matrix",
#'                       indices_path = "exp1_10kb_abs.bed",
#'                       sample_name = "exp1")
#' # Data loaded from Juicer
#' exp2 <- load_contacts("exp2_10kb.cooler",
#'                       balancing = TRUE,
#'                       sample_name = "exp2")
#' synched <- sync_indices(list(exp1, exp2))
#' }
sync_indices <- function(explist) {
  
  nexp <- length(explist)
  if (nexp < 1L) {
    return(NULL)
  }
  valid_exps <- vapply(explist, inherits, logical(1), "contacts")
  if (!all(valid_exps)) {
    stop("Can only synchronise indices of GENOVA contacts objects.")
  }
  if (nexp == 1L) {
    message(paste("Single experiment provided to 'sync_indices()'.",
                  "No indices to synchronise. Returning input."))
    return(explist)
  }
  
  res <- resolution(explist)
  if (!all(tail(res, -1) == res[[1]])) {
    stop("Can only synchronise indices of a single resolution.")
  }
  
  idxs <- copy(lapply(explist, `[[`, "IDX"))
  
  # Initialise template from first exp
  template <- idxs[[1]]
  setnames(template, 4, "exp1")
  
  # For-loop because we update template with subsequent exps
  for (i in 2:nexp) {
    jointhis <- idxs[[i]]
    setnames(jointhis, 4, paste0("exp", i))
    # Join on chrom/start
    template <- merge.data.table(template, idxs[[i]], 
                                 on = c("V1", "V2"), all = TRUE)
    # Adopt missing ends from next exp
    template$V3 <- ifelse(is.na(template$V3.x), template$V3.y, template$V3.x)
    template[, V3.x := NULL]
    template[, V3.y := NULL]
  }
  template[, newbin := seq_len(NROW(template))]
  final <- template[, list(V1, V2, V3, V4 = newbin)]
  
  target <- template[["newbin"]]
  expcols <- paste0("exp", seq_len(nexp))
  new_explist <- lapply(seq_len(nexp), function(i) {
    source <- template[, eval(as.symbol(expcols[i]))]
    keep <- !is.na(source)
    
    # Make match vectors
    match <- vector("integer", max(source, na.rm = TRUE))
    match[source[keep]] <- target[keep]
    
    exp <- copy(explist[[i]])
    
    # Set index
    exp$IDX <- final
    
    # Lift over old indices to new
    exp$MAT[, c("V1", "V2") := list(match[V1], match[V2])]
    setkeyv(exp$MAT, c("V1", "V2"))
    
    # Also for centromere entries
    exp$CENTROMERES[, c("start", "end") := list(match[start], match[end])]
    setkeyv(exp$CENTROMERES, "chrom")
    
    return(exp)
  })
  names(new_explist) <- names(explist)
  return(new_explist)
}
