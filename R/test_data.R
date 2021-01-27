
#' Test data
#'
#' Retrieves a test \code{contacts} object from a local cache at either 40k or
#' 150k resolution. The test data only contains human chromosomes 21 and 22.
#'
#' @param res The resolution of test data to use. Either \code{"40k"} or
#'   \code{"150k"}.
#' @param download A \code{logical(1)}. If \code{FALSE} (default), will not
#'   download the data if it cannot be found in the cache. If \code{TRUE}, data
#'   is downloaded when not found in the cache.
#'
#' @return A GENOVA \code{contacts} at 40k or 150k resolution.
#' @export
#'
#' @details The data originates from wildtype Hap1 cells from the Haarhuis
#'   \emph{et al}. (2017) publication. The full data is available in raw fastq
#'   files or valid pair files at GEO accession GSE95014. The data were mapped
#'   against the 'hg19' reference genome using HiC-Pro v2.7.7.
#'
#'   The test data only needs to be downloaded only once as it is stored in a
#'   local cache. The local cache can be cleared using the
#'   \code{\link[GENOVA]{erase_GENOVA_cache}()} function.
#'
#' @examples
#' \dontrun{
#' # Fails the first time
#' exp <- get_test_data("40k")
#'
#' # Should succees if `download = TRUE`
#' exp <- get_test_data("40k", download = TRUE)
#'
#' # From now on, it should succeed
#' exp <- get_test_data("40k")
#'
#' # Additional calls with `download = TRUE` do not re-download the data.
#' exp <- get_test_data("40k", download = TRUE)
#' }
get_test_data <- function(res, download = FALSE) {
  res <- switch(
    as.character(res),
    "40k" = "data/test_40k.rds",
    "150k" = "data/test_150k.rds",
    error("No such optional data exists.", .call = FALSE)
  )
  
  path <- get_GENOVA_data_filepath(res, mustWork = FALSE)
  if (!file.exists(path)) {
    if (literalTRUE(download)) {
      check <- download_GENOVA_data()
      path <- get_GENOVA_data_filepath(res, mustWork = TRUE)
    } else {
      stop("Data is not in cache and `download` is FALSE")
    }
  }
  exp <- readRDS(path)
  setDT(exp$MAT)
  setDT(exp$IDX)
  setDT(exp$CENTROMERES)
  setkeyv(exp$MAT, c("V1", "V2"))
  setkeyv(exp$IDX, c("V1", "V2", "V3"))
  setkeyv(exp$CENTROMERES, c("chrom", "start", "end"))
  return(exp)
}

#' Erase cache.
#'
#' GENOVA can download test data with the
#' \code{\link[GENOVA]{get_test_data}()} function. This function undoes that
#' by deleting the files in the local cache. Test data can then be re-downloaded
#' using the \code{get_test_data(..., download = TRUE)} function.
#'
#' @return \code{0} for success and \code{1} for failure.
#' @export
#'
#' @examples
#' NULL
erase_GENOVA_cache <- function() {
  try_require("pkgfilecache", "download_GENOVA_data", "CRAN")
  pkg_info <- pkgfilecache::get_pkg_info("GENOVA")
  erased <- pkgfilecache::erase_file_cache(pkg_info)
  return(erased)
}


download_GENOVA_data <- function() {
  try_require("pkgfilecache", "download_GENOVA_data", "CRAN")
  
  pkg_info <- pkgfilecache::get_pkg_info("GENOVA")
  
  local_filenames <- list(
    c("data", "test_150k.rds"),
    c("data", "test_40k.rds")
  )
  
  urls <- c(
    "https://github.com/teunbrand/GENOVA/raw/testdat_branch/data/test_150k.rds", 
    "https://github.com/teunbrand/GENOVA/raw/testdat_branch/data/test_40k.rds"
  )
  
  md5sums <- c(
    "8b89fa2344a544ec01f7218d8dcbde86",
    "6d011f0fe5632f24bc97497d26df791b"
  )
  
  res <- pkgfilecache::ensure_files_available(
    pkg_info, local_filenames, urls, md5sums = md5sums
  )
  res$file_status <- NULL
  return(res)
}

list_GENOVA_data <- function() {
  try_require("pkgfilecache", "download_GENOVA_data", "CRAN")
  pkg_info <- pkgfilecache::get_pkg_info("GENOVA")
  available <- pkgfilecache::list_available(pkg_info)
  return(available)
}

get_GENOVA_data_filepath <- function(filename, mustWork = TRUE) {
  try_require("pkgfilecache", "download_GENOVA_data", "CRAN")
  pkg_info <- pkgfilecache::get_pkg_info("GENOVA")
  path <- pkgfilecache::get_filepath(pkg_info, filename, mustWork = mustWork)
  return(path)
}


