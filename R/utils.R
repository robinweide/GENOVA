#' Match bed-like entries to Hi-C bin indices
#'
#' @param IDX The IDX-slot of a \code{contacts} object
#' @param bed A 3-column data.frame
#' @param mode A \code{character} of length 1 indicating what position of the
#'   \code{bed} argument to match with the indices. Possible values:
#'   \code{"center"}, \code{"start"} or \code{"end"}.
#'   
#' @return An \code{integer} vector of length \code{nrow(bed)} and parallel to
#'   \code{bed} with indices to the Hi-C matrix.
#'
#' @details Out of bounds values are matched to nearest bin.
#' @export

bed2idx <- function(IDX, bed, mode = c("centre", "start", "end")) {
  if (!inherits(bed, "data.frame") | is.data.table(bed)) {
    bed <- as.data.frame(bed)
  }
  if (anyNA(bed[1:3])) {
    stop("Cannot match `NA`s to indices.", call. = FALSE)
  }

  # American/British spelling
  mode <- gsub("center", "centre", mode)
  mode <- match.arg(mode, c("centre", "start", "end"))

  # Reformat bed depending on mode
  bed <- cbind.data.frame(
    V1 = bed[, 1],
    V2 = switch(mode,
      "centre" = (bed[, 2] + bed[, 3]) / 2,
      "start" = bed[, 2],
      "end" = bed[, 3]
    )
  )
  
  class(IDX) <- "data.frame"

  # Assign entries to shared chromosomes
  chroms <- intersect(IDX[, 1], bed[, 1])
  bed_group <- match(bed[, 1], chroms)
  IDX_group <- match(IDX[, 1], chroms)
  
  # Split by chromosome
  bed_chrom <- split(bed[, 2], bed_group)
  IDX_chrom <- split(IDX[, c(2, 4)], IDX_group)
  
  # Match bed entry to idx
  out <- mapply(function(i, j) {
    j[pmax(findInterval(i, j[, 1]), 1), 2]
  }, i = bed_chrom, j = IDX_chrom, SIMPLIFY = FALSE)
  unsplit(out, bed_group)
}

#' Upper triangle sparse symmetric triplet matrix to dense matrix
#'
#' Convenience function for returning square, dense symmetric matrices from the
#' upper triangle of triplet format symmetric matrices.
#'
#' @param x A \code{numeric} vector containing the data for matrix elements.
#' @param i A \code{integer} vector parallel to \code{x} containing row
#'   positions.
#' @param j A \code{integer} vector parallel to \code{x} containing column
#'   positions.
#' @param dim A \code{integer} of length 1 specifying the dimensions of the
#'   output matrix.
#' @param offset A \code{integer} of length 1 noting a potential offset in the
#'   \code{i} and \code{j} arguments.
#'
#'   Particularly useful when lookup up square regions around the diagonal in
#'   Hi-C data.
#'
#' @note For reasons of speed, no checks are performed wether the input is
#'   compatible with sensible output.
#'
#' @return A \code{dim * dim} sized \code{matrix}.
#'
#' @keywords internal
dt_matrix <- function(x, i, j, dim, offset) {
  m <- matrix(0, dim, dim)
  m[matrix(c(i, j, j, i) - offset, 2 * length(i), 2)] <- c(x, x)
  m
}

#' Get a matrix from a BED-like entry
#'
#' Extracts a square symmetric matrix around the diagonal from a \code{contacts}
#' object.
#'
#' @param exp The \code{contacts} objects produced by
#'   \code{construct.experiment()}
#' @param chrom A \code{character} of length 1: the chromosome.
#' @param start An \code{integer} of length 1 noting the start position in bp.
#' @param end An \code{integer} of length 1 noting the end position in bp.
#'
#' @return A list with the \code{X} and \code{Y} coordinates and a \code{matrix
#'   Z} containing the contacts at these coordinates.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get the TP53-locus in an experiment mapped to hg19
#' mat <- select_subset(WT, "chr17", 75e5, 76e5)
#'
#' # Plot the region
#' image(mat)
#' }
select_subset <- function(exp, chrom, start, end) {
  # Restrict data.table core usage
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)

  idx <- which(exp$IDX[, 1] == chrom & exp$IDX[, 2] >= start & exp$IDX[, 2] <= end)
  pos <- rowMeans(exp$IDX[idx, 2:3])
  i <- exp$IDX[idx, V4]
  min <- i[1] - 1
  len <- length(i)
  structure(list(
    x = pos,
    y = pos,
    z = exp$MAT[CJ(V1 = i, V2 = i),
                dt_matrix(V3, V1, V2, len, min),
                nomatch = NULL]
  ), class = c("contacts_matrix", "list"), 
  chrom = chrom, resolution = resolution(exp))
}

#' @export
#' @noRd
as.matrix.contacts_matrix <- function(x, ...) {
  out <- x$z
  dimnames(out) <- list(x$x, x$y)
  out
}

# taken from ggplot
try_require <- function(package, fun, source = NULL) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }

  if (source == 'BIOC') {
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from Bioconductor and try again.", call. = FALSE)
  } else   if (source == 'github') {

    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from github and try again.", call. = FALSE)
  } else {
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install and try again.", call. = FALSE)
  }

}

#' Check compatability of a list of GENOVA experiments
#'
#' Checks if the indices (ABS slot) across experiments are identical.
#'
#' @inheritParams APA
#'
#' @return A \code{list} of GENOVA experiment(s).
#'
#' @keywords internal
check_compat_exp <- function(explist) {
  if (!is.list(explist)) {
    stop("Supply either a GENOVA experiment or list of GENOVA experiments",
         call. = FALSE
    )
  }
  
  # Re-list of only one experiment was given
  if (any(c("MAT", "IDX") %in% names(explist))) {
    explist <- list(explist)
  }
  
  # Test equality of experiments in list
  if (length(explist) > 1) {
    equal <- vapply(seq_along(explist)[-1], function(i) {
      literalTRUE(all.equal(explist[[1]]$IDX, explist[[i]]$IDX))
    }, logical(1))
    
    if (any(!equal)) {
      stop(paste(
        "Indices of experiment(s)",
        paste(which(!equal) + 1, collapse = " & "),
        "are not equal to indices of experiment 1"
      ), call. = FALSE)
    }
  }
  
  return(explist)
}



# Equivalent to isTRUE fron R>3.5
literalTRUE <- function(x) is.logical(x) && length(x) == 1L && !is.na(x) && x

GENOVA_THEME <- function() {
  p = ggplot2::theme(panel.background = ggplot2::element_blank(),
                     legend.key =  ggplot2::element_rect(fill = 'white'),
                     strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                     panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                     text = ggplot2::element_text(color = 'black'),
                     axis.text = ggplot2::element_text(colour = 'black'),
                     strip.text = ggplot2::element_text(colour = 'black') )
  return(p)
}

cache_chroms <- function(exp) {
  first <- exp$IDX[, list(V4 = min(V4)), by = V1][order(V4)]
  
  chrom <- findInterval(exp$MAT$V1, first$V4)
  cis <- findInterval(exp$MAT$V2, first$V4) == chrom
  rle <- rle(paste0(chrom, "-", cis))
  x <- as.data.table(tstrsplit(rle$values, "-"))
  x <- x[, list(chrom = first$V1[as.integer(V1)],
                cis = as.logical(V2),
                lengths = rle$lengths)]
  x[, ends := cumsum(x$lengths)]
  x[, starts := x$ends - x$lengths + 1]
}

check_input <- function(x) {
  cl <- sys.call(-1)
  today <- format(Sys.Date(), "%d-%m") == "01-04"
  if (today & !("antoni" %in% names(cl))) {
    image(GENOVA::antoni, col = grDevices::grey.colors(255))
  }
  return(invisible(!today))
}

#' Standardising genomic location inputs
#'
#' The main raison d'etre for this function is to standardise genomic locations
#' from user input to a data.frame with the names 'chrom', 'start', 'end'.
#'
#' @param chrom One of the following:
#'  * A character with a chromosome name (e.g. "chr1")
#'  * A data.frame with 'chrom', 'start', 'end' equivalent columns.
#'  * A character with a complete locus (e.g. "chr1:1000-2000"). Handled by the 
#'    interpret_location_function.
#' @param start An integer
#' @param end An integer
#' @param check A logical of length 1. Should we check wether 'chrom' is
#'   character, start is numeric and end is numeric and non contain NAs?
#' @param singular A logical of length 1. Is the input expected to describe one 
#'   and only one locus?
#'
#' @return A data.frame with the columns 'chrom', 'start' and 'end'
#' @noRd
#' @examples \dontrun{
#' # Should handle single loci
#' standardise_location("chr1:1000-2000")
#' standardise_location(data.frame(x = "chr1", y = 1000, z = 2000))
#' standardise_location("chr1", 1000, 2000)
#' 
#' # Should handle vectorised as well
#' standardise_location(c("chr1:1000-2000", "chr2:3000-4000"))
#' standardise_location(data.frame(x = c("chr1", "chr2"), 
#'                                 y = c(1000, 3000), 
#'                                z = c(2000, 4000)))
#' standardise_location(c("chr1", "chr2"), c(1000, 3000), c(2000, 4000))
#' 
#' # If a single locus is expected, throw a warning when multiple are given
#' standardise_location(c("chr1:1000-2000", "chr2:3000-4000"), singular = TRUE)
#' }
standardise_location <- function(chrom, start = NULL, end = NULL, 
                                 check = TRUE, singular = FALSE) {
  if (is.atomic(chrom)) {
    chrom <- as.character(chrom)
    test <- strsplit(chrom, "\\:|\\-")
    if (all(lengths(test) == 1)) {
      # chrom-start-end branch
      out <- data.frame(chrom = chrom, start = start, end = end)
    } else {
      # Likely "chr1:100-200" format branch
      if (!all(lengths(test) == 3)) {
        stop("Cannot interpret character location.", call. = FALSE)
      }
      out <- interpret_location_string(chrom, feature_id = FALSE)
    }
  } else {
    # bed branch
    if (!inherits(chrom, "data.frame") | is.data.table(chrom)) {
      out <- as.data.frame(chrom)[1:3]
    } else {
      out <- as.data.frame(chrom[1:3]) # To convert tibbles, see #188
    }
    if (!missing(start) | !missing(end)) {
      message("The input in `chrom` has been prioritised over ",
              "`start` and `end`.")
    }
    setnames(out, 1:3, c("chrom", "start", "end"))
  }
  
  if (check) {
    # Check for NAs
    if (anyNA(out)) {
      stop("Location cannot contain `NA`s.", call. = FALSE)
    }
    # Check chromosome
    if (!is.character(out[[1]])) {
      if (is.factor(out[[1]])) {
        out[[1]] <- as.character(out[[1]])
      } else {
        stop("The chromosome name must be a character.", call. = FALSE)
      }
    }
    # Check start-end
    # Shouldn't check for integers specifically, since -Inf and Inf are
    # special values part of the doubles universe, and mean take the entire
    # chromosome.
    if (!is.numeric(out[[2]]) | !is.numeric(out[[3]])) {
      stop("The starts and ends of the location must be numeric.")
    }
  }
  
  if (singular) {
    if (NROW(out) > 1L) {
      stop("A single location is expected, while multiple were detected.", 
           call. = FALSE)
    }
  }
  return(out)
}

#' Interpreter for string locations
#' 
#' This tries to interpret locations given as strings
#'
#' @param location One of the following:
#' * A character of locations, e.g. c("chr1:100-200", "chr2:200-300")
#' * A character of paired locations, e.g. c("chr1:100-200;chr2:200-300")
#' * A character of bin indices, e.g. "1111"
#' * A character of paired bin indices, e.g. "1111,2222"
#' @param IDX The IDX data.table in a contacts object for coupling bin indices
#'   to locations. Should not apply for location characters.
#' @param feature_id A logical of length 1. Should we attach IDs for unique 
#'   locations?
#'
#' @return A data.frame with 'chrom', 'start' or 'end' columns. If paired 
#' locations were given, these columns occur twice with "_up" and "_down" appended.
#' @noRd
#' @examples \dontrun{
#' interpret_location_string(c("chr1:100-200", "chr2:200-300"))
#' interpret_location_string(c("chr1:100-200;chr2:200-300"))
#' interpret_location_string(c("1111", "2222"), IDX = exp$IDX)
#' interpret_location_string(c("1111,2222"), IDX = exp$IDX)
#' }
interpret_location_string <- function(location, IDX = NULL, feature_id = TRUE) {
  location <- as.character(location)
  if (feature_id) {
    fid <- match(location, unique(location))
  }
  location <- gsub(",|\\.", "", location)
  location <- tstrsplit(location, ";")
  location <- lapply(location, function(str) {
    str <- tstrsplit(str, "\\:|\\-")
    if (length(str) == 3) {
      # Assume bed-entry
      df <- data.frame(
        chrom = str[[1]],
        start = as.integer(str[[2]]),
        end   = as.integer(str[[3]])
      )
      return(df)
    } else if (length(str) == 1) {
      # Assume IDX index
      if (is.null(IDX)) {
        warning("No IDX found. Returning location as bins.", 
                call. = FALSE)
        df <- data.frame(location = as.integer(str[[1]]))
        return(df)
      } else {
        str <- as.integer(str[[1]])
        df <- as.data.frame(IDX[.(str), on = "V4"][, V4 := NULL])
        colnames(df) <- c("chrom", "start", "end")
      }
      return(df)
    } else {
      warning("Cannot interpret location. Returning location as-is.", 
              call. = FALSE)
      df <- as.data.frame(str)
      colnames(df) <- paste0("V", seq_len(ncol(df)))
      return(df)
    }
  })
  
  if (length(location) == 2) {
    # Assume paired-end locations
    colnames(location[[1]]) <- paste0(colnames(location[[1]]), "_up")
    colnames(location[[2]]) <- paste0(colnames(location[[2]]), "_down")
  }
  
  location <- do.call(cbind, location)
  if (feature_id) {
    location <- cbind(location, feature_id = fid)
  }
  return(location)
}

import_bigwig <- function(file, chr, start, end) {
  try_require("rtracklayer", "import_bigwig", "Bioconductor")
  loc <- IRanges::IRangesList(IRanges::IRanges(start, end))
  names(loc) <- chr
  im <- rtracklayer::import(file, selection = loc,
                            as = "RleList")[[chr]]
  im <- data.table(len  = im@lengths, val = im@values)
  im[, end := cumsum(len)]
  im[, start := end - len + 1]
  x <- integer(nrow(im) * 2)
  x[c(TRUE, FALSE)] <- im$start
  x[c(FALSE, TRUE)] <- im$end
  x <- pmax(pmin(x, end), start)
  data.table(x = x, y = rep(im$val, each = 2))
}
