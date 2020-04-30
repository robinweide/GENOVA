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
  list(x = pos,
       y = pos,
       z = exp$MAT[CJ(V1 = i, V2 = i),
                   dt_matrix(V3, V1, V2, len, min),
                   nomatch = NULL])
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

GENOVA_THEME = function(){
  p = ggplot2::theme(panel.background = ggplot2::element_blank(),
                     legend.key =  ggplot2::element_rect(fill = 'white'),
                     strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                     panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                     text = ggplot2::element_text(color = 'black'),
                     axis.text = ggplot2::element_text(colour = 'black'),
                     strip.text = ggplot2::element_text(colour = 'black') )
  return(p)
}



# Resolution methods ------------------------------------------------------

#' Get the resolution of GENOVA objects
#'
#' @param x An object from the GENOVA package
#'
#' @return A \code{integer} with the resolution
#' @export
#'
#' @examples
#' \dontrun{
#' resolution(contacts_object)
#' resolution(discovery_object)
#' }
resolution <- function(x) {
  UseMethod("resolution")
}

#' @export
#' @method resolution default
resolution.default <- function(x) {
  attr(x, "resolution", exact = TRUE)
}

#' @export
#' @method resolution list
resolution.list <- function(x) {
  ans <- lapply(x, attr, which = "resolution", exact = TRUE)
  ans <- lapply(ans, function(y) {
    if(is.null(y)) {
      return(NA_integer_)
    } else {
      as.numeric(y)
    }
  })
  ans <- unlist(ans)
  ans
}

# Expnames methods --------------------------------------------------------

# This is basically needed for all the analysis function to have a consistent
# way to prevent naming inconsistencies between functions such as in Github
# issue #175 and #153

# Getters


#' @name expnames
#' @title Sample names for an object
#'
#' @description Looks in the appropriate spot for sample names and returns or sets these.
#'
#' @param x The object for which to retrieve sample names.
#' @param simplify A \code{logical} of length 1: if the \code{x} argument is a
#'   list, should the sample names be returned as a flat character?
#' @param value A non-\code{NA} \code{character} vector of appropriate length.
#'
#' @return A \code{character} or \code{list} with sample names
#' @export
#'
#' @examples
#' \dontrun{
#' expnames(discovery_object)
#' expnames(contacts_object)
#' }
expnames <- function(x, simplify = TRUE) {
  UseMethod("expnames")
}

#' @export
#' @method expnames default
expnames.default <- function(x, simplify = TRUE) {
  NULL
}

#' @export
#' @method expnames contacts
expnames.contacts <- function(x, simplify = TRUE) {
  attr(x, "samplename", exact = TRUE)
}

#' @export
#' @method expnames ARMLA_discovery
expnames.ARMLA_discovery <- function(x, simplify = TRUE) {
  tail(dimnames(x[[1]]), 1)[[1]]
}

#' @export
#' @method expnames list
expnames.list <- function(x, simplify = TRUE) {
  ans <- lapply(x, expnames)
  ans <- unlist(ans)
  ans
}

# Setters

#' @export
#' @rdname expnames
`expnames<-` <- function(x, value) {
  if (any(is.na(value))) {
    stop("No new name can be NA.",
         call. = FALSE)
  }
  if (any(!is.character(value))) {
    stop("New names should be of type `character`.",
         call. = FALSE)
  }
  UseMethod("expnames<-")
}

#' @export
#' @method `expnames<-` default
`expnames<-.default` <- function(x, value) {
  x
}

#' @export
#' @method `expnames<-` contacts
`expnames<-.contacts` <- function(x, value) {
  if (length(value) != 1L) {
    stop("The new expname should be length 1.",
         call. = FALSE)
  }
  attr(x, "samplename") <- value
  x
}

#' @export
#' @method `expnames<-` ARMLA_discovery
`expnames<-.ARMLA_discovery` <- function(x, value) {
  dim <- dim(x[["signal"]])
  if (length(value) != tail(dim, 1)) {
    stop("The new expnames should be of the same length", 
         " as the existing expnames", call. = FALSE)
  }
  x[] <- lapply(x, function(y) {
    if (is.array(y) || is.matrix(y)) {
      dimnames(y)[[length(dim(y))]] <- value
    }
    return(y)
  })
  x
}