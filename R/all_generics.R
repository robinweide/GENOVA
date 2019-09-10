#' @export
#' @rdname APA_PESCAn_visualisation
visualise <- function(discovery, ...) {
  UseMethod("visualise", discovery)
}
