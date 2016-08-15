#' Read hicpro three column matrix format.
#'
#' This function loads a HiC-pro file as a matrix. It assumes a three-column
#' layout: bin1, bin2 and score. All scors are normalised to contacts per *norm*
#' total contacts.
#'
#' @param file Full path to file.
#' @param norm Normalising factor.
#' @return A data.table with normalised counts.
#' @export
read.hicpro.matrix <- function(file, norm=100e6){
  data <- data.table::fread(file)
  data.table::setkey(data, "V1", "V2")
  data$V3 <- norm*data$V3/sum(data$V3)
  return(data)
}
