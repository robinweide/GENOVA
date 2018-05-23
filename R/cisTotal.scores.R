#' cisTotal.scores
#'
#' Calculate the genome-wide ratio of the intrachromosomal interactions over
#' the interchromosomal interactions
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp The Hi-C experiment object of a sample:
#' produced by construct.experiment().
#' @return The genome-wide cis/total ratio
#' @import data.table
#' @examples
#' # calculate the ratios for a 1Mb-resolution experiment
#' cis_out <- cisTotal.scores(WT_1MB)
#' @export
cisTotal.scores <- function(exp){
  #sort the TADs otherwise the findInterval function will not work
  cis <- 0
  #loop over chromosomes
  d <- exp$ICE #select the ICE matrix from the HiC data structure
  #select the chromosome names in a specific order
  exp$ABS[,1] <- as.character(exp$ABS[,1])
  chrom <- c(exp$ABS[1,1],
             exp$ABS[which(head(exp$ABS[,1],-1) != tail(exp$ABS[,1],-1))+1,1])
  for( chr in chrom ){
    #select the chromosome ids from the HiC data structure
    chr.id <- exp$ABS[exp$ABS[,1]==chr,]
    #select the chromosomal indexes for chr
    chrom.min <- min(chr.id[,4]); chrom.max <- max(chr.id[,4])
    #remove the trans contacts if only cis interactions are scored
    d.chrom <- d[(d$V1 >= chrom.min & d$V1 <= chrom.max),]
    d.chrom <- d.chrom[d.chrom$V2 >= chrom.min & d.chrom$V2 <= chrom.max,]
    cis <- cis + sum(d.chrom$V3)
  }
  total <- sum(d$V3)
  cis/total
}
