
# Documentation -----------------------------------------------------------

#' @title Quantification of results
#' @name quantify
#' @description Here be a placholder description

# Functions ---------------------------------------------------------------

# If you fail it is best to fail graciously

#' @export
#' @rdname quantify
#' @usage NULL
quantify.default <- function(discovery, ...) {
  stop("No quantify method for class '", class(discovery),
       "' has been implemented.", call. = FALSE)
}
