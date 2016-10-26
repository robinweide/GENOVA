#' cisTotal.perChrom
#'
#' Calculate the ratio of the intrachromosomal interactions over the interchromosomal interactions
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @return The per-chromosome cis/total ratio
#' @note The per chromsosome score is on average lower than the whole genome score, because in the per chromosome analysis every ditag is counted twice (for the two interacting chromosomes)
#' @import data.table
cisTotal.perChrom <- function( exp ){
  #sort the TADs otherwise the findInterval function will not work
  cis <- c()
  #loop over chromosomes
  d <- exp$ICE #select the ICE matrix from the HiC data structure
  #select the chromosome names in a specific order
  exp$ABS[,1] <- as.character(exp$ABS[,1])
  chrom <- c(exp$ABS[1,1],exp$ABS[which(head(exp$ABS[,1],-1) != tail(exp$ABS[,1],-1))+1,1])
  for( chr in chrom ){
    chr.id <- exp$ABS[exp$ABS[,1]==chr,] #select the chromosome ids from the HiC data structure
    #select the minimal and maximal chromosomal indexes for chr
    chrom.min <- min(chr.id[,4]); chrom.max <- max(chr.id[,4])
    #remove the trans contacts if only cis interactions are scored
    d.cis <- d[(d$V1 >= chrom.min & d$V1 <= chrom.max) & (d$V2 >= chrom.min & d$V2 <= chrom.max),]
    d.total <- d[(d$V1 >= chrom.min & d$V1 <= chrom.max) | (d$V2 >= chrom.min & d$V2 <= chrom.max),]
    cis <- c(cis, sum(d.cis$V3)/sum(d.total$V3))
  }
  cis
}
