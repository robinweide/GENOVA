# Robin's colours ---------------------------------------------------------

#' @name GENOVA_colours
#' @title Colours in GENOVA
#' 
#' @description GENOVA comes with a few built-in colour palettes. The gradients
#' in GENOVA's palettes have been carefully chosen to be perceptually linear,
#' which is great for Hi-C visualisations.
#' 
#' @section Sequential palettes:
#' 
#' The sequential palettes are used for displaying absolute values, such as the
#' (normalised) contacts values in Hi-C matrices, or the average of aggregates.
#' The default palette is the 'Hot' palette. The default sequential palette can 
#' be changed by setting the global options.
#' 
#' \preformatted{
#' 
#' options("GENOVA.colour.palette" = "whitered")
#' 
#' options("GENOVA.colour.palette" = "hot")
#' 
#' }
#' 
#' \subsection{Hot}{
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{colour-hot.png}{options: style="width:255px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{colour-hot.png}\out{\end{center}}
#'   }
#' }
#' 
#' \subsection{White-Red}{
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{colour-whitered.png}{options: style="width:255px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{colour-whitered.png}\out{\end{center}}
#'   }
#' }
#' 
#' @section Diverent palettes:
#' 
#' The divergent palettes are used for displaying relative values, such as
#' Z-score normalised contacts, differences and fold changes.
#' 
#' \subsection{Divergent}{
#' 
#' The default divergent palette is called \code{"divergent"} and goes from a blue
#' at low values, to a light grey at the midpoint to red at high values.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{colour-divergent.png}{options: style="width:255px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{colour-divergent.png}\out{\end{center}}
#'   }
#' }
#' 
#' \subsection{Green-Pink}{
#' 
#' A secondary divergent palette, \code{"greenpink"} is only used when two objects in
#' a plot require divergent palettes, but need to be discriminated from 
#' oneanother. 
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{colour-greenpink.png}{options: style="width:255px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{colour-greenpink.png}\out{\end{center}}
#'   }
#' }
#' 
#' @section Details: 
#' The colours in ggplot based visualisations, such as 
#' \code{\link[GENOVA]{visualise}()} and \code{\link[GENOVA]{pyramid}()}, are 
#' based on the 
#' \code{\link[=GENOVA_colour_scales]{scale_(colour|fill)_GENOVA}()} and 
#' \code{\link[=GENOVA_colour_scales]{scale_(colour|fill)_GENOVA_div}()}
#' functions.
#' 
#' For visualisations in base R, such as \code{\link[GENOVA]{hic_matrixplot}()} 
#' and \code{\link[GENOVA:GENOVA_image]{image()}},
#' the \code{\link[grDevices]{colorRampPalette}()} 
#' function is used.
NULL

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

# Colour choices now controlled by global options, e.g.
#
# options("GENOVA.colour.palette" = "hot")
# options("GENOVA.colour.palette" = "whitered")
#
# Should propagate to downstream functions
.choose_palette <- function(x = NULL) {
  if (is.null(x)) {
    x <- getOption("GENOVA.colour.palette", default = "hot")
  }
  if (length(x) == 1L && is.character(x)) {
    x <- match.arg(x, c("hot", "whitered"))
    x <- switch(x,
                hot = bezier_corrected_hot,
                whitered = bezier_corrected_whiteRed
                # Reserved room for future palettes
    )
  }
  x
}

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
#'   \code{"divergent"} or \code{"greenpink"} for divergent scales. In 
#'   sequential scales, if \code{palette = NULL} an attempt will be made to find
#'   \code{"GENOVA.colour.palette"} in the global options.
#' @param midpoint A \code{numeric} indicating where the midpoint in the
#'   divergent scale should occur.
#'   
#' @seealso The \code{\link[GENOVA:GENOVA_colours]{GENOVA colours}} description.
#'
#' @return A \code{ScaleContinuous} gg-object.
#'
#' @examples
#' \dontrun{
#' require(ggplot2)
#' p <- ggplot(faithfuld, 
#'             aes(waiting, eruptions, fill = density)) +
#'   geom_tile()
#' p + scale_fill_GENOVA()
#' p + scale_fill_GENOVA_div()
#' }
NULL

#' @rdname GENOVA_colour_scales
#' @export
scale_colour_GENOVA <- 
  function(..., palette = NULL, na.value = "grey50", guide = "colourbar",
           aesthetics = "colour") 
  {
    if (!is.function(palette)) {
      palette <- .choose_palette(palette)
      palette <- scales::gradient_n_pal(palette)
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
