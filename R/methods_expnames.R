# Expnames methods --------------------------------------------------------

# This is basically needed for all the analysis function to have a consistent
# way to prevent naming inconsistencies between functions such as in Github
# issue #175 and #153

# Getters -----------------------------------------------------------------

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

#' @export
#' @method expnames genomescore_discovery
expnames.genomescore_discovery <- function(x, simplify = TRUE) {
  x <- names(x[[1]])
  setdiff(x, c("window", "position", "mid", "bin",
               "chrom", "start", "end", "chromosome"))
}

#' @export
#' @method expnames saddle_discovery
expnames.saddle_discovery <- function(x, simplify = TRUE) {
  unique(x$saddle$exp)
}

#' @export
#' @method expnames RCP_discovery
expnames.RCP_discovery <- function(x, simplify = TRUE) {
  levels(x$raw$samplename)
}

#' @export
#' @method expnames IIT_discovery
expnames.IIT_discovery <- function(x, simplify = TRUE) {
  tail(colnames(x$results), -2)
}

#' @export
#' @method expnames chrommat_discovery
expnames.chrommat_discovery <- function(x, simplify = TRUE) {
  dimnames(x$obs)[[3]]
}

# Setters -----------------------------------------------------------------

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

#' @export
#' @method `expnames<-` genomescore_discovery
`expnames<-.genomescore_discovery` <- function(x, value) {
  i <- names(x[[1]])
  i <- which(!(i %in% c("window", "position", "mid", "bin",
                        "chrom", "start", "end", "chromosome")))
  if (length(value) != length(i)) {
    stop("The new expnames should be of the same length", 
         " as the existing expnames", call. = FALSE)
  }
  if (inherits(x, "data.frame")) {
    names(x)[i] <- value
  } else {
    names(x[[1]])[i] <- value
  }
  return(x)
}

#' @export
#' @method `expnames<-` saddle_discovery
`expnames<-.saddle_discovery` <- function(x, value) {
  oldnames <- expnames(x)
  if (length(value) != length(oldnames)) {
    stop("The new expnames should be of the same length",
         " as the existing expnames", call. = FALSE)
  }
  value <- setNames(value, oldnames)
  x$saddle$exp <- value[x$saddle$exp]
  return(x)
}

#' @export
#' @method `expnames<-` RCP_discovery
`expnames<-.RCP_discovery` <- function(x, value) {
  oldnames <- expnames(x)
  if (length(value) != length(oldnames)) {
    stop("The new expnames should be of the same length",
         " as the existing expnames", call. = FALSE)
  }
  levels(x$raw$samplename) <- value
  if ("smooth" %in% names(x)) {
    levels(x$smooth$samplename) <- value
  }
  return(x)
}

#' @export
#' @method `expnames<-` chrommat_discovery
`expnames<-.chrommat_discovery` <- function(x, value) {
  oldnames <- expnames(x)
  if (length(value) != length(oldnames)) {
    stop("The new expnames should be of the same length",
         " as the existing expnames", call. = FALSE)
  }
  dimnames(x$obs)[[3]] <- value
  dimnames(x$exp)[[3]] <- value
  return(x)
}
