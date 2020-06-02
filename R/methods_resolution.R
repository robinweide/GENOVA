# Resolution methods ------------------------------------------------------

#' Get the resolution of GENOVA objects
#'
#' @param x An object from the GENOVA package
#' @param ... Optional arguments (not used).
#'
#' @return A \code{integer} with the resolution
#' @export
#' 
#' @note If the \code{ggplot2} namespace is loaded after GENOVA, use 
#' \code{GENOVA::resolution()} on GENOVA objects.
#'
#' @examples
#' \dontrun{
#' resolution(contacts_object)
#' resolution(discovery_object)
#' }
resolution <- function(x, ...) {
  UseMethod("resolution")
}

#' @export
#' @method resolution default
resolution.default <- function(x, ...) {
  ggplot2::resolution(x, ...)
}

#' @export
#' @method resolution discovery
resolution.discovery <- function(x, ...) {
  attr(x, "resolution", exact = TRUE)
}

#' @export
#' @method resolution contacts
resolution.contacts <- function(x, ...) {
  attr(x, "resolution", exact = TRUE)
}

#' @export
#' @method resolution list
resolution.list <- function(x, ...) {
  ans <- lapply(x, resolution)
  ans <- lapply(ans, function(y) {
    if(is.null(y) || !is.finite(y)) {
      return(NA_integer_)
    } else {
      as.integer(y)
    }
  })
  ans <- unlist(ans)
  ans
}

# Contact matrix is a list, so we bypass the list resolution method
#' @export
#' @method resolution contact_matrix
resolution.contact_matrix <- function(x, ...) {
  attr(x, "resolution", exact = TRUE)
}

# What follows is protection against ggplot2::resolution(GENOVA_object)
# It's not an ideal solution but better than freezing your R session
# It'll output the correct resolution anyway, but with warnings.

#' @export
#' @noRd
#' @method as.double contacts
as.double.contacts <- function(x, ...) {
  warning("Cannot coerce a 'contacts' object to numeric.", call. = FALSE)
  return(c(0, GENOVA::resolution(x)))
}

#' @export
#' @noRd
#' @method range contacts
range.contacts <- function(x, ...) {
  warning("Cannot compute range on a 'contacts' object.", call. = FALSE)
  return(c(0, GENOVA::resolution(x)))
}

#' @export
#' @noRd
#' @method as.double discovery
as.double.discovery <- function(x, ...) {
  warning("Cannot coerce a 'discovery' object to numeric.", call. = FALSE)
  return(c(0, GENOVA::resolution(x)))
}

#' @export
#' @noRd
#' @method range discovery
range.discovery <- function(x, ...) {
  warning("Cannot compute range on a 'discovery' object.", call. = FALSE)
  return(c(0, GENOVA::resolution(x)))
}