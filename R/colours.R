# Robin's colours ---------------------------------------------------------

bezier_corrected_hot <- c("#ffffff","#ffd4aa","#f4a86e","#db8047",
                              "#bb5c2e",
                              "#943e1f","#682616","#3a160e","#000000")
bezier_corrected_divergent <- c('#0065b2', '#2399de', '#67cdfe', 
                        '#f5f5f5', 
                        '#fcb09d', '#ed6855', '#be2a21')
bezier_corrected_whiteRed <- c('#ffffff', '#ffeae2', '#ffd5c6', '#ffbea9', 
                         '#ffa78c', 
                         '#ff8e6f', '#ff7251', '#ff4f30', '#ff0000')
bezier_corrected_greenPink <- c('#38793a','#589b59','#7abd79','#a8dfa6',
                                '#f5f5f5',
                                '#fbbffb','#f28bf2','#d364d5','#b23fb6')

# ggplot scale functions --------------------------------------------------

#' @name GENOVA_colour_scales
#' @title GENOVA gradient colour scales
#'
#' @description These scale functions create colour gradients that have been
#'   created to be perceptually linear for the sequential scales and
#'   perceptually triangular for the divergent scales.
#'
#' @inheritParams ggplot2::continuous_scale
#' @inheritDotParams ggplot2::continuous_scale
#' @param palette A \code{character} vector of length 1 for palette options.
#'   Either \code{"hot"} or \code{"whitered"} for sequential scales, or
#'   \code{"divergent"} or \code{"greenpink"} for divergent scales.
#' @param midpoint A \code{numeric} indicating where the midpoint in the
#'   divergent scale should occur.
#'
#' @return A \code{ScaleContinuous} gg-object.
#'
#' @examples
#' NULL
NULL

#' @rdname GENOVA_colour_scales
#' @export
scale_colour_GENOVA <- 
  function(..., palette = "hot", na.value = "grey50", guide = "colourbar",
           aesthetics = "colour") 
  {
    palette <- match.arg(palette, c("hot", "whitered"))
    if (palette == "hot") {
      palette <- scales::gradient_n_pal(bezier_corrected_hot)
    } else {
      palette <- scales::gradient_n_pal(bezier_corrected_whiteRed)
    }
    ggplot2::continuous_scale(
      aesthetics, "colour", palette,
      na.value = na.value, guide = guide, ...
    )
  }

#' @rdname GENOVA_colour_scales
#' @export
scale_fill_GENOVA <-
  function(..., aesthetics = "fill") {
    scale_colour_GENOVA(aesthetics = aesthetics, ...)
  }

#' @rdname GENOVA_colour_scales
#' @export
scale_colour_GENOVA_div <- 
  function(..., palette = "divergent", midpoint = NA,
           na.value = "grey50", guide = "colourbar", aesthetics = "colour") 
    
  {
    palette <- match.arg(palette, c("divergent", "greenpink"))
    if (palette == "divergent") {
      palette <- scales::gradient_n_pal(bezier_corrected_divergent)
    } else {
      palette <- scales::gradient_n_pal(bezier_corrected_greenPink)
    }
    
    if (is.na(midpoint)) {
      rescaler <- scales::rescale
    } else {
      rescaler <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
        scales::rescale_mid(x, to, from, midpoint)
      }
    }
    ggplot2::continuous_scale(
      aesthetics, "colour_div", palette, na.value = na.value,
      guide = guide, ..., rescaler = rescaler
    )
  }

#' @rdname GENOVA_colour_scales
#' @export
scale_fill_GENOVA_div <- function(..., aesthetics = "fill") {
  scale_colour_GENOVA_div(..., aesthetics = aesthetics)
}
