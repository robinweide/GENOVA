#' Read hicpro three column matrix format.
#'
#' This function loads a HiC-pro file as a matrix. It assumes a three-column
#' layout: bin1, bin2 and score. All scors are normalised to contacts per *scale_bp*
#' total contacts.
#'
#' @param file Full path to file.
#' @param scale_bp Normalising factor. Set to NULL to skip scale_bp.
#' @return A data.table with normalised counts.
read.hicpro.matrix <- function(file, scale_bp = 1e9) {
  data <- data.table::fread(file)
  data.table::setkey(data, "V1", "V2")
  if (!is.null(scale_bp)) {
    data$V3 <- scale_bp * data$V3 / sum(data$V3)
  }
  return(data)
}
