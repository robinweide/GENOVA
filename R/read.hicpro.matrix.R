#' Read hicpro three column matrix format.
#'
#' This function loads a HiC-pro file as a matrix. It assumes a three-column
#' layout: bin1, bin2 and score. All scors are normalised to contacts per *norm*
#' total contacts.
#'
#' @param file Full path to file.
#' @param norm Normalising factor.
#' @return A data.table with normalised counts.
read.hicpro.matrix <- function(file, norm=1e9){
  data <- data.table::fread(file)
  data$V3[!is.finite(data$V3)] <- 0 # Juicer-output will have entries with nan as value. Set these to 0
  data.table::setkey(data, "V1", "V2")
  data$V3 <- norm*data$V3/sum(data$V3)
  return(data)}
