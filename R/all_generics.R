#' @export
#' @rdname visualise
visualise <- function(discovery, ...) {
  UseMethod("visualise", discovery)
}

#' @export
quantify <- function(discovery, ...) {
  UseMethod("quantify", discovery)
}

#' @export
#' @rdname bundle
bundle <- function(..., collapse = "_") {
  UseMethod("bundle", list(...)[[1]])
}

#' @export
#' @rdname unbundle
unbundle <- function(discovery, ...) {
  UseMethod("unbundle", discovery)
}
