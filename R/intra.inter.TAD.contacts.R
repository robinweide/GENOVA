#' intra.inter.TAD.contacts 
#'
#' Combine a bed file with TAD positions and the Hi-C matrix file to calculate the average coverage over a TAD and between two TADs.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param TAD A bed-dataframe.
#' @param max.neighbor How many surrounding TADs should be taken into account?
#' @return A list containing a dataframe with scores between two TADs and an updated version of the TAD-input file, which now also includes TAD-IDs.
#' @import data.table
#' @export
#'
intra.inter.TAD.contacts <- function ( exp, TAD, max.neighbor = 10 ){
  #error handling
  if( any(TAD[,3] < TAD[,2])){ stop("Incorrectly structured TAD data: end column smaller than start") }

  cat("Calculating the intra- and inter-TAD contact frequency with", max.neighbor, "neighboring TADs. This can take several minutes depending on the size of the dataset\n")
  #sort the TADs otherwise the findInterval function will not work
  TAD <- TAD[order(TAD[,1],TAD[,2]),1:3]
  TAD[,4] <- NA
  chrom <- unique(TAD[,1])
  max.tad <- 1
  if(exists("intra.inter")){ rm(intra.inter) }
  #loop over chromosomes
  for( chr in chrom ){
    cat("Now analyzing", chr, "\r" )
    chr.id <- exp$ABS[exp$ABS[,1]==chr,] #select the chromosome ids from the HiC data structure
    chr.id[,5] <- findInterval(chr.id[,2], TAD[TAD[,1]==chr,2]) #determine which ids overlap with which TAD/domain
    chr.id[,6] <- findInterval(chr.id[,2], TAD[TAD[,1]==chr,3]) #determine which ids overlap with which TAD/domain
    chr.id <- chr.id[chr.id[,5] != chr.id[,6],] #select only those sites that are within a TAD that is given by the user
    d <- exp$ICE #select the ICE matrix from the HiC data structure
    #select the intrachromosomal contact sites for chr
    chrom.min <- min(chr.id[,4]); chrom.max <- max(chr.id[,4])
    d.chrom <- d[J(chrom.min:chrom.max),]
    d.chrom <- d.chrom[d.chrom$V2 >= chrom.min & d.chrom$V2 <= chrom.max,]

    #associate hic-pro indices with TAD positions
    tad.x <- chr.id[match( d.chrom$V1, chr.id[,4]),5]
    tad.y <- chr.id[match( d.chrom$V2, chr.id[,4]),5]
    #select tad(combination)s that are not more than max.neighbor apart
    sel <- tad.y - tad.x <= max.neighbor & tad.x > 0 & tad.y > 0 & !is.na(tad.x) & !is.na(tad.y)
    tad.x <- tad.x[sel]
    tad.y <- tad.y[sel]
    d.chrom <- d.chrom$V3[sel]
    #calculate the sum of the TAD-TAD interactions
    if(exists("intra.inter") ){
      max.tad <- max(intra.inter[,1:2])
      #sum over d.chrom, on the two tad positions
      intra.inter.chr <- ave(d.chrom, tad.x, tad.y, FUN=sum)
      intra.inter.chr <- data.frame(x=tad.x, y=tad.y, score=intra.inter.chr) #create a new data.frame with the results
      intra.inter.chr <- unique(intra.inter.chr)
      intra.inter.chr[,1] <- intra.inter.chr[,1] + max.tad
      intra.inter.chr[,2] <- intra.inter.chr[,2] + max.tad
      intra.inter <- rbind(intra.inter, intra.inter.chr)
      TAD[TAD[,1]==chr,4] <- (max.tad+1):(max.tad+nrow(TAD[TAD[,1]==chr,]))
    }else{
      intra.inter <- ave(d.chrom, tad.x, tad.y, FUN=sum)
      intra.inter <- data.frame(x=tad.x, y=tad.y, score=intra.inter)
      intra.inter <- unique(intra.inter)
      TAD[TAD[,1]==chr,4] <- max.tad:(max.tad+nrow(TAD[TAD[,1]==chr,])-1)
    }
    #give every TAD an ID corresponding to the IDs used in intra.inter
  }
  #combine the TAD positions and 
  list(hic=intra.inter, tad=TAD)
}
