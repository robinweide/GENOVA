#' cisTotal.perChrom
#'
#' Calculate the ratio of the intrachromosomal interactions over the
#' interchromosomal interactions.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param exp The Hi-C experiment object of a sample:
#' produced by construct.experiment().
#' @param chromsToUse Only get data of a small set of chromosomes.
#' Genome-wide score will not we given!
#' @return A boxplot with red line showing genome-wide %cis and a list:
#' @return \item{perChrom}{dataframe with chromosome and percentage cis}
#' @return \item{genomeWide}{the genome-wide percentage cis}
#' @note The per chromsosome score is on average lower than the
#' whole genome score, because in the per chromosome analysis every ditag is
#' counted twice (for the two interacting chromosomes)
#' @import data.table
#' @examples
#' # calculate the ratios for a 1Mb-resolution experiment
#' cisChrom_out <- cisTotal.perChrom(WT_1MB)
#'
#' # using the underlying data, make a barplot
#' plot( cisChrom_out$perChrom, las=2 )
#' abline( h = cisChrom_out$genomeWide, col = 'red' )
#' @export
cisTotal.perChrom <- function(exp , chromsToUse = NULL, ...){
  if(is.null(chromsToUse)){
    chromsToUse = exp$CHRS
  }
  #sort the TADs otherwise the findInterval function will not work
  cis <- c()
  cisSum = 0
  #d <- exp$ICE #select the ICE matrix from the HiC data structure
  #select the chromosome names in a specific order
  exp$ABS[,1] <- as.character(exp$ABS[,1])

  # loop over all chromosomes.
  # store cis and trans in df with chromname

  for( chr in exp$CHRS ){
    #select the chromosome ids from the HiC data structure
    chr.id <- exp$ABS[exp$ABS[,1]==chr,]
    #select the minimal and maximal chromosomal indexes for chr
    chrom.min <- min(chr.id[,4]); chrom.max <- max(chr.id[,4])

    # get singnal-rows with indexes
    allRows = exp$ICE[V1 %in% chrom.min:chrom.max | V2 %in% chrom.min:chrom.max]

    # get cis
    tmp <- allRows[V1 %in% chrom.min:chrom.max & V2 %in% chrom.min:chrom.max,3]
    cisRows <- 0
    if(!nrow(tmp) == 0){
      cisRows = sum(tmp)
      cisSum = cisSum + cisRows
    }
    
    # get trans
    transRows = sum(allRows[xor(V1 %in% chrom.min:chrom.max,
                                V2 %in% chrom.min:chrom.max),3])

    # calc perc: cis*2 because we count trans double: chr1:chr2 contacts will
    # be used in chr1 and chr2
    # also, cis-matrices are not symetrical
    percCIS = (cisRows*2)/sum(cisRows*2,transRows)

    names(percCIS) = chr
    cis <- c(cis, percCIS )
  }

  cis = cis[names(cis) %in% chromsToUse]

  # genomewide
  total <- sum(exp$ICE$V3)
  GW = 100*(cisSum/total)


  # plot!
  df =  data.frame(chromosomes = names(cis),
                   `percentage cis` = as.numeric(cis)*100)
  df$chromosomes = factor(df$chromosomes, levels = exp$CHRS)
  bxpdat <- boxplot(df[,2],
                    ylim = c(0,100),
                    ylab = 'Percentage Cis',
                    pch = 20,
                    plot = F)
  boxplot(df[,2],ylim = c(0,100), ylab = 'Percentage Cis', pch = 20)

  abline(h = GW, col = 'red', lty = 3)
  text(.5,GW, labels = round(GW,2) , adj = c(0,-.5), col = 'red')

  if(length(df[which(df[,2] %in% bxpdat$out),1]) > 0){
    text(bxpdat$group,
         bxpdat$out,
         df[which(df[,2] %in% bxpdat$out),1],
         pos = c(2,4)[rep(c(1,2), length(names(cis)))[order(bxpdat$out)]])
  }

  invisible(list(perChrom = df, genomeWide = GW))

}
