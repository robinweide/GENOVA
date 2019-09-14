#' @export
#' @rdname visualise
visualise <- function(discovery, ...) {
  UseMethod("visualise", discovery)
}

#' @export
quantify <- function(discovery, ...) {
  UseMethod("quantify", discovery)
}
