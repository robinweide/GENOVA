#' Chromosome matrix
#'
#' Computes sums of chromosome interactions genome-wide along with expected 
#' values, making it easier to identify enriched trans-contacts that may 
#' indicate translocations.
#'
#' @param explist Either a single GENOVA \code{contacts} object or a list of
#'   GENOVA \code{contacts} objects.
#' @param expected A \code{character} of length one, one of \code{"bins"},
#'   \code{"sums"}, \code{"trans"} or \code{"cis"}, determining how the expected
#'   component is calculated. See details.
#' @param include_chr A \code{character} containing the chromosome names to
#'   include. Setting this to \code{"all"} will refrain form dropping
#'   chromosomes based on inclusion criteria, but will still apply the exclusion
#'   criteria..
#' @param exclude_chr A \code{character} containing the chromosome names to
#'   exclude. When a chromosome is both in the \code{include_chr} argument and
#'   in the \code{exclude_chr} argument, exlusion criteria overrule inclusion
#'   criteria.
#' @param sort_chr A \code{logical} of length one deciding wether the chromosome
#'   names should by sorted by their name.
#'
#' @details The \code{expected} part can be calculated various ways, and is returned as
#'   proportion of the total. Setting \code{expected = "bins"} (the default)
#'   will calculate the expected part as the product of a chromosome pair's
#'   number of bins. Setting \code{expected = "sums"} calculates the expected as
#'   the chi squared null hypothesis, namely that the proportions are
#'   conditional on the sums over individual chromosomes. Likewise,
#'   \code{expected = "trans"} and \code{expected = "cis"} will calculate the
#'   same as \code{expected = "sum"}, but omit cis- and trans-combinations 
#'   from the calculation respectively.
#'   
#' @return A \code{chrommat_discovery} object with 2 elements:
#' @return \itemize{
#'   \item{\strong{\code{obs}}, a \code{n} x \code{n} x \code{length(explist)} 
#'   \code{array}, wherein \code{n} is the number of included chromosomes, 
#'   containing sums of contacts for chromosome combinations}
#'   \item{\strong{\code{exp}}, an \code{array} of the same dimensions of 
#'   \code{obs}, but contains the expected values}
#'   }
#' @export
#'
#' @section Resolution recommendation: 500kb-1Mb+
#'
#' @examples
#' \dontrun{
#' cm <- chromosome_matrix(list(WT = WT_1Mb, KO = KO_1Mb))
#' visualise(cm)
#' }
chromosome_matrix <- function(
  explist, 
  expected = c("bins", "sums", "trans", "cis", "regress"), 
  include_chr = "all",
  exclude_chr = c("chrM", "chrY", "M", "Y"),
  sort_chr = TRUE
) {
  expected <- match.arg(expected)
  
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  # Set experiment parameters
  explist <- check_compat_exp(explist)
  expnames <- expnames(explist)
  n_exp <- length(explist)
  
  # Curate chromosomes
  idx <- copy(explist[[1]]$IDX)
  chrom <- character(max(idx[["V4"]]))
  chrom[idx[["V4"]]] <- idx[["V1"]]
  if (length(include_chr) == 1 && include_chr != "all") {
    chrom[!(chrom %in% include_chr)] <- NA_character_
  }
  chrom[chrom %in% exclude_chr] <- NA_character_
  valid_chrom <- unique(chrom)
  valid_chrom <- valid_chrom[nzchar(valid_chrom) & !is.na(valid_chrom)]
  
  # Count bins per chromosome
  chrom_bins <- idx[, list(N = .N), by = "V1"][, setNames(N, V1)]
  chrom_bins <- chrom_bins[names(chrom_bins) %in% valid_chrom]
  chrom_bins <- outer(chrom_bins, chrom_bins)
  # Express as proportion of total bins
  chrom_bins <- chrom_bins / sum(chrom_bins)
  
  # Sort if asked
  if (sort_chr) {
    chrsort <- gsub("M", "Z", gsub("chr", "", dimnames(chrom_bins)[[1]]))
    chrsort <- paste0(ifelse(grepl("[^0-9]", chrsort), "9", "0"), chrsort)
    chrsort <- substr(chrsort, nchar(chrsort) - 1, nchar(chrsort))
    chrom_bins <- chrom_bins[order(chrsort), order(chrsort)]
  }
  
  # Dummy template to fill later
  template <- chrom_bins
  template[] <- 0
  
  # Calculate sums per chromosome combination
  counts <- vapply(explist, function(xp) {
    dat <- xp$MAT
    dat <- dat[, list(V1 = chrom[V1], V2 = chrom[V2], V3 = V3)]
    dat <- dat[!is.na(V1) & !is.na(V2)]
    dat <- dat[, list(V3 = sum(V3)), by = c("V1", "V2")]
    with(dat, `[<-`(template, i = cbind(c(V1, V2), c(V2, V1)), c(V3, V3)))
  }, template)
  
  if (expected == "sums") {
    # Like the chisq.test expected
    E <- vapply(seq_len(n_exp), function(i) {
      m <- counts[, , i, drop = TRUE]
      m <- outer(colSums(m), rowSums(m), "*")
      m / sum(m)
    }, template)
  } else if (expected == "trans") {
    # Like the chisq.test expected, but setting cis to 0
    E <- vapply(seq_len(n_exp), function(i) {
      m <- counts[, , i, drop = TRUE]
      diag(m) <- 0
      m <- outer(colSums(m), rowSums(m), "*")
      m / sum(m)
    }, template)
  } else if (expected == "cis") {
    E <- vapply(seq_len(n_exp), function(i) {
      m <- counts[, , i, drop = TRUE]
      m <- diag(m)
      m <- outer(m, m)
      m / sum(m)
    }, template)
  } else if (expected == "regress"){
    E <- vapply(seq_len(n_exp), function(i) {
      m <- counts[, , i, drop = TRUE]
      long <- data.frame(
        x = as.factor(row(m)),
        y = as.factor(col(m)),
        value = log10(as.vector(m)),
        bins = as.vector(chrom_bins)
      )
      keep <- with(long, x != y)
      fit <- lm(value ~ 0 + x + y, long, keep)
      m[] <- predict(fit, newdata = long)
      10^m / sum(10^m, na.rm = TRUE)
    }, template)
  } else {
    # Just match up the proportional bins with the observed
    E <- vapply(seq_len(n_exp), function(i){chrom_bins}, template)
  }
  
  # Output formatting stuff
  dimnames(counts)[[3]] <- expnames
  dimnames(E)[[3]] <- expnames
  
  structure(list(obs = counts, exp = E),
            resolution = resolution(explist)[1],
            package = "GENOVA",
            mode = expected,
            class = c("chrommat_discovery", "discovery"))
}


