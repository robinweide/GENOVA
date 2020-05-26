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
  ggplot2::resolution(x)
}

#' @export
#' @method resolution discovery
resolution.discovery <- function(x) {
  attr(x, "resolution", exact = TRUE)
}

#' @export
#' @method resolution contacts
resolution.contacts <- function(x) {
  attr(x, "resolution", exact = TRUE)
}

#' @export
#' @method resolution list
resolution.list <- function(x) {
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
resolution.contact_matrix <- function(x) {
  attr(x, "resolution", exact = TRUE)
}