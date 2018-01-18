#' cisTotal.perChrom
#'
#' Calculate the ratio of the intrachromosomal interactions over the interchromosomal interactions
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chromsToUse Only get data of a small set of chromosomes. Genome-wide score will not we given!
#' @return The per-chromosome cis/total ratio
#' @note The per chromsosome score is on average lower than the whole genome score, because in the per chromosome analysis every ditag is counted twice (for the two interacting chromosomes)
#' @import data.table
#' @export
cisTotal.perChrom <- function( exp , chromsToUse = NULL, ...){
  #sort the TADs otherwise the findInterval function will not work
  cis <- c()
  cisSum = 0
  #d <- exp$ICE #select the ICE matrix from the HiC data structure
  #select the chromosome names in a specific order
  exp$ABS[,1] <- as.character(exp$ABS[,1])
  chrom <- NULL
  if(is.null(chromsToUse)){
    chrom <- exp$CHRS
  } else {
    chrom <- chromsToUse
  }

  for( chr in chrom ){
    chr.id <- exp$ABS[exp$ABS[,1]==chr,] #select the chromosome ids from the HiC data structure
    #select the minimal and maximal chromosomal indexes for chr
    chrom.min <- min(chr.id[,4]); chrom.max <- max(chr.id[,4])

    # get singnal-rows with indexes
    allRows = exp$ICE[V1 %in% chrom.min:chrom.max | V2 %in% chrom.min:chrom.max]

    # get cis
    cisRows = sum(allRows[V1 %in% chrom.min:chrom.max & V2 %in% chrom.min:chrom.max,3])
    cisSum = cisSum + cisRows
    # get trans
    transRows = sum(allRows[xor(V1 %in% chrom.min:chrom.max , V2 %in% chrom.min:chrom.max),3])

    # calc perc: cis*2 because we count trans double: chr1:chr2 contacts will be used in chr1 and chr2
    # also, cis-matrices are not symetrical
    percCIS = (cisRows*2)/sum(cisRows*2,transRows)


    cis <- c(cis, percCIS )
  }

  # genomewide
  total <- sum(exp$ICE$V3)
  GW = 100*(cisSum/total)

  # plot!
  df =  data.frame(chromosomes = chrom, `percentage cis` = as.numeric(cis)*100)
  df$chromosomes = factor(df$chromosomes, levels = chrom)
  bxpdat <- boxplot(df[,2],ylim = c(0,100), ylab = 'Percentage Cis', pch = 20, plot = F)
  boxplot(df[,2],ylim = c(0,100), ylab = 'Percentage Cis', pch = 20)
  if(is.null(chromsToUse)){
    abline(h = GW, col = 'red', lty = 3)
    text(.5,GW, labels = round(GW,2) , adj = c(0,-.5), col = 'red')
  }
  text(bxpdat$group,
       bxpdat$out,
       df[which(df[,2] %in% bxpdat$out),1],
       pos = c(2,4)[rep(c(1,2), length(chrom))[order(bxpdat$out)]])

  if(is.null(chromsToUse)){
    invisible(list(perChrom = df, genomeWide = GW))
  } else {
    invisible(list(perChrom = df, genomeWide = NA))
  }
}


