#' chromosomeMatrix
#'
#' Plot a genome-wide overview-matrix: values are log10(observed/expected).
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp The Hi-C experiment object of a sample: produced by
#' construct.experiment().
#' @param color.fun Optional color-function
#' @param cut.off Adjust the maximum value of the color-scale
#' @param chromsToUse Only get data of a small set of chromosomes.
#' Genome-wide score will not we given!
#' @param remove A vector of (partial) names of chromosomes to leave out.
#' @param plot A \code{logical}: plot the resulting matrix?
#' @note Please use only low-resolution matrices with this (e.g. 1Mb resolution)
#' @return A plot of the chromosome matrix and an assignable list:
#' @return \item{rawCounts}{the observed chromosome-chromosome matrix}
#' @return \item{normMat}{the expected chromosome-chromosome matrix}
#' @examples
#' \dontrun{
#' # use chromosomeMatrix without chromosome M and Y.
#' chromosomeMatrix(Hap1_WT_1MB, remove = c("chrM", "chrY"), cut.off = 2)
#' }
#' @export
chromosomeMatrix <- function(exp, color.fun = NULL, cut.off = NULL,
                             chromsToUse = NULL, remove = NULL, plot = TRUE) {
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  exp$IDX[["V1"]] <- as.character(exp$IDX[["V1"]])
  chrom <- c(
    exp$IDX[1, V1],
    exp$IDX[which(head(exp$IDX[["V1"]], -1) != tail(exp$IDX[["V1"]], -1)) + 1, V1]
  )

  chrA <- factor(exp$IDX[exp$MAT[["V1"]], V1], levels = chrom)
  chrB <- factor(exp$IDX[exp$MAT[["V2"]], V1], levels = chrom)

  # calculate the number of interactions per/between chromosome(s)
  chrom.count <- aggregate(V3 ~ chrA + chrB, data = exp$MAT, sum)
  chrom.count.inv <- chrom.count[chrom.count[, 1] != chrom.count[, 2], c(2, 1, 3)]
  names(chrom.count.inv) <- names(chrom.count)
  chrom.count <- rbind(chrom.count, chrom.count.inv)

  # make sure all the chromosomes are covered
  all.combinations <- cbind(
    rep(chrom, each = length(chrom)),
    rep(chrom, length(chrom))
  )
  chrom.count <- merge(all.combinations, chrom.count,
    by = c(1, 2),
    all.x = T, sort = F
  )
  chrom.count[is.na(chrom.count[, 3]), 3] <- 0

  # for some reason two merge operation are necessary to get
  # the correct order
  chrom.count <- merge(all.combinations, chrom.count,
    by = c(1, 2),
    all.x = T, sort = F
  )

  # create a normalization matrix
  chrom.n <- table(factor(exp$IDX[["V1"]], levels = chrom))
  norm.mat <- outer(chrom.n, chrom.n, "*")
  chrom.mat <- matrix(chrom.count[, 3], ncol = length(chrom), byrow = T)

  if (!is.null(chromsToUse)) {
    keep.vec <- which(rownames(norm.mat) %in% chromsToUse)

    chrom.mat <- chrom.mat[keep.vec, keep.vec]
    norm.mat <- norm.mat[keep.vec, keep.vec]
  }

  if (!is.null(remove)) {
    if (!any(remove %in% exp$CHRS)) {
      warning("Chromosome names in remove do not match exp!")
    } else {
      remove.vec <- which(rownames(norm.mat) %in% remove)

      chrom.mat <- chrom.mat[-remove.vec, -remove.vec]
      norm.mat <- norm.mat[-remove.vec, -remove.vec]
    }
  }


  mat <- list(rawCounts = chrom.mat, normMat = norm.mat)
  if (plot) {
    chrom.comparison.plot(mat, color.fun, cut.off)
  }
  invisible(mat)
}

