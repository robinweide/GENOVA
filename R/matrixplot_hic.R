getEmptyCols <- function(mat, Nmad = 2, slop = 0.9, lowDiag = 2, highDiag = 10, zeroNA = T) {
  if (zeroNA) {
    mat$z[mat$z == 0] <- NA
  }
  
  delta <- row(mat$z) - col(mat$z)
  mat$z[abs(delta) <= lowDiag | abs(delta) >= highDiag] <- NA
  
  tres <- median(mat$z[!is.na(mat$z)]) - (Nmad * mad(mat$z[!is.na(mat$z)]))
  
  possibles <- table(which(mat$z < tres, arr.ind = TRUE)[, 2])
  possiblespossibles <- possibles[possibles >= ceiling(max(colSums(!is.na(mat$z))) * slop)]
  
  foundBins <- as.numeric(names(possiblespossibles))
  
  return(foundBins)
}

# https://stackoverflow.com/questions/29105175/find-neighbouring-elements-of-a-matrix-in-r
gregor <- function(mat) {
  n <- nrow(mat)
  mat.pad <- rbind(NA, cbind(NA, mat, NA), NA)
  ind <- 2:(n + 1) # row/column indices of the "middle"
  neigh <- rbind(
    SE = as.vector(mat.pad[ind + 1, ind + 1]),
    NW = as.vector(mat.pad[ind - 1, ind - 1])
  )
  return(neigh)
}


fillNAs <- function(mat, MADtreshold = 2) {
  
  # get all NA-locations
  eCols <- getEmptyCols(mat, Nmad = MADtreshold)
  
  # make a smoothed version
  completeSmooth <- gregor(mat$z)
  completeSmooth[completeSmooth == min(completeSmooth, na.rm = T)] <- NA
  CSmean <- apply(completeSmooth, 2, mean, na.rm = T)
  CSmat <- matrix(CSmean, nrow = nrow(mat$z), ncol = ncol(mat$z))
  
  # fill in the blanks
  mat$z[eCols, ] <- CSmat[eCols, ]
  mat$z[, eCols] <- CSmat[eCols, ]
  
  return(mat)
}

draw_exon <- function(genes, chrom, y.pos, width, rotate = F) {
  genes.chr <- genes[genes$chrom == chrom, ]
  
  es <- as.numeric(unlist(strsplit(as.character(genes.chr$exonStart), ",")))
  ee <- as.numeric(unlist(strsplit(as.character(genes.chr$exonEnd), ",")))
  
  strands <- rep(genes.chr$strand, unlist(lapply(strsplit(as.character(genes.chr$exonStart), ","), length)))
  
  add <- ifelse(strands == "+", 0.025 * width + width, -0.025 * width - width)
  add2 <- ifelse(genes.chr$strand == "+", 0.025 * width + width, -0.025 * width - width)
  if (rotate) {
    rect((y.pos - width / 2) + add, es, (y.pos + width / 2) + add, es, col = "black")
    segments(y.pos + add2, genes.chr$txStart, y.pos + add2, genes.chr$txEnd)
  } else {
    rect(es, (y.pos - width / 2) + add, ee, (y.pos + width / 2) + add, col = "black")
    segments(genes.chr$txStart, y.pos + add2, genes.chr$txEnd, y.pos + add2)
  }
}

plot.triangle <- function(bed, chrom, y1, y2, start, end, rotate = F) {
  y.scale <- abs(y1 - y2) / 0.15
  x.wid <- (end - start) * 0.012 * y.scale
  
  
  # first plot positive
  sel.bed <- bed[bed[, 1] == chrom &
                   bed[, 3] < end &
                   bed[, 2] > start &
                   bed[, 6] == "+", ]
  if (nrow(sel.bed) > 0) {
    col <- "red"
    add <- ifelse(sel.bed[, 6] == "+", x.wid, -x.wid)
    add.list <- unique(add)
    if (rotate) {
      triangle <- lapply(sel.bed[, 2], function(x) list(
        yy = c(x, x, x + add.list),
        xx = c(y1, y2, (y1 + y2) / 2)
      ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col = col, border = col))
    } else {
      # with polygon
      triangle <- lapply(sel.bed[, 2], function(x) list(
        xx = c(x, x, x + add.list),
        yy = c(y1, y2, (y1 + y2) / 2)
      ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col = col, border = col))
    }
  }
  
  # then repeat for negative
  sel.bed <- bed[bed[, 1] == chrom &
                   bed[, 3] < end &
                   bed[, 2] > start &
                   bed[, 6] == "-", ]
  if (nrow(sel.bed) > 0) {
    col <- "blue"
    add <- ifelse(sel.bed[, 6] == "+", x.wid, -x.wid)
    add.list <- unique(add)
    if (rotate) {
      triangle <- lapply(sel.bed[, 2], function(x) list(
        yy = c(x, x, x + add.list),
        xx = c(y1, y2, (y1 + y2) / 2)
      ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col = col, border = col))
    } else {
      # with polygon
      triangle <- lapply(sel.bed[, 2], function(x) list(
        xx = c(x, x, x + add.list),
        yy = c(y1, y2, (y1 + y2) / 2)
      ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col = col, border = col))
    }
  }
  
  # check if there are any without orientation and plot these as rectangles
  sel.bed <- bed[bed[, 1] == chrom &
                   bed[, 3] < end &
                   bed[, 2] > start &
                   bed[, 6] != "-" &
                   bed[, 6] != "+", ]
  if (nrow(sel.bed) > 0) {
    col <- "black"
    add <- ifelse(sel.bed[, 6] != "+", 0, 0)
    add.list <- unique(add)
    if (rotate) {
      rect(y1, sel.bed[, 2], y2, sel.bed[, 3], col = col, border = col)
    } else {
      rect(sel.bed[, 2], y1, sel.bed[, 3], y2, col = col, border = col)
    }
  }
}

plot.rectangle <- function(bed, chrom, y1, y2, start, end, col, rotate = F) {
  sel.bed <- bed[bed[, 1] == chrom & bed[, 2] < end & bed[, 3] > start, ]
  if (nrow(sel.bed) > 0) {
    if (rotate) {
      rect(y1, sel.bed[, 2], y2, sel.bed[, 3], col = col, border = col)
    } else {
      rect(sel.bed[, 2], y1, sel.bed[, 3], y2, col = col, border = col)
    }
  }
}

plot.genes <- function(genes, chrom, start, end, y.pos, rotate = F) {
  draw_exon(genes, chrom = chrom, y.pos = y.pos, width = 0.05, rotate = rotate)
}



features <- function(mat1, chrom, yMax = NULL, genes = NULL, chip1 = NULL,
                     chip2 = NULL, autoCHIP = T, rotate = F,
                     type = c("triangle", "triangle"), col = NULL) {
  
  # - set bounds for three tracks
  if (is.null(col)) {
    col <- "black"
  }
  if (!is.null(genes)) {
    # BW in 1,  bed/null in 2
    if (typeof(chip1) == "character" && typeof(chip2) != typeof(chip1)) {
      gene.pos <- 1.05
      chip1.y1 <- 0
      chip1.y2 <- 0.5
      chip2.y1 <- 0.6
      chip2.y2 <- 0.8
      # BW in 2,  bed/null in 1
    } else if (typeof(chip2) == "character" && typeof(chip2) != typeof(chip1)) {
      gene.pos <- 1.05
      chip1.y1 <- 0
      chip1.y2 <- 0.2
      chip2.y1 <- 0.3
      chip2.y2 <- 0.8
    } else {
      gene.pos <- 1.05
      chip1.y1 <- 0
      chip1.y2 <- 0.375
      chip2.y1 <- 0.425
      chip2.y2 <- 0.8
    }
  } else {
    # BW in 1,  bed/null in 2
    if (typeof(chip1) == "character" && typeof(chip2) != typeof(chip1)) {
      gene.pos <- 0
      chip1.y1 <- 0
      chip1.y2 <- 0.675
      chip2.y1 <- 0.725
      chip2.y2 <- 1.1
      # BW in 2,  bed/null in 1
    } else if (typeof(chip2) == "character" && typeof(chip2) != typeof(chip1)) {
      gene.pos <- 0
      chip1.y1 <- 0
      chip1.y2 <- 0.375
      chip2.y1 <- 0.425
      chip2.y2 <- 1.1
    } else {
      gene.pos <- 0
      chip1.y1 <- 0
      chip1.y2 <- 0.5
      chip2.y1 <- 0.7
      chip2.y2 <- 1.1
    }
  }
  y.values <- list(
    gene.pos = gene.pos, chip1.y1 = chip1.y1, chip1.y2 = chip1.y2,
    chip2.y1 = chip2.y1, chip2.y2 = chip2.y2
  )
  
  
  # make chip smaller
  
  if (max(mat1$x) - min(mat1$x) > 2e6) {
    if (typeof(chip1) == "list") {
      newSize <- diff(unlist(y.values[2:3])) * 0.75
      D <- diff(unlist(y.values[2:3])) - newSize
      # y.values[[2]] = y.values[[2]] + D/2
      y.values[[3]] <- y.values[[3]] - D
    }
    
    if (typeof(chip2) == "list") {
      newSize <- diff(unlist(y.values[4:5])) * 0.75
      D <- diff(unlist(y.values[4:5])) - newSize
      y.values[[4]] <- y.values[[4]] + D
      # y.values[[5]] = y.values[[5]] - D/2
    }
  }
  
  
  
  # - plot an empty canvas c(0,1.3)
  
  if (rotate) {
    plot(0,
         type = "n", ylim = rev(range(mat1$x)), xlim = c(1.3, 0), axes = F,
         xlab = "", ylab = ""
    )
  } else {
    plot(0,
         type = "n", xlim = range(mat1$x), ylim = c(0, 1.3), axes = F,
         xlab = "", ylab = ""
    )
  }
  
  
  # - plot inner track c(0.1,0.4) -> chipInner
  
  if (typeof(chip1) == "list") { # BED!
    
    if (ncol(chip1) < 6) {
      type[1] <- "rectangle"
      # there is no orientation in column 6
    } else if (autoCHIP & !any(c("-", "+") %in% chip1[, 6])) {
      type[1] <- "rectangle"
    }
    
    if (type[1] == "triangle") {
      # +
      tmp_step <- (y.values$chip1.y2 - y.values$chip1.y1) / 2
      tmp_step <- tmp_step * 0.9
      blabla <- plot.triangle(chip1[chip1[, 6] == "+", ],
                              chrom = chrom,
                              y1 = y.values$chip1.y1,
                              y2 = y.values$chip1.y1 + tmp_step,
                              start = min(mat1$x), end = max(mat1$x), rotate = rotate
      )
      # -
      blabla <- plot.triangle(chip1[chip1[, 6] == "-", ],
                              chrom = chrom,
                              y1 = y.values$chip1.y2 - tmp_step,
                              y2 = y.values$chip1.y2, start = min(mat1$x),
                              end = max(mat1$x), rotate = rotate
      )
      # # other
      blabla <- plot.rectangle(chip1[chip1[, 6] != "-" & chip1[, 6] != "+", ],
                               chrom = chrom, y1 = y.values$chip1.y1,
                               y2 = y.values$chip1.y2, start = min(mat1$x),
                               end = max(mat1$x), col = col[1], rotate = rotate
      )
    } else if (type[1] == "rectangle") {
      blabla <- plot.rectangle(chip1,
                               chrom = chrom, y1 = y.values$chip1.y1,
                               y2 = y.values$chip1.y2, start = min(mat1$x),
                               end = max(mat1$x), col = col[1], rotate = rotate
      )
    }
  } else if (typeof(chip1) == "character") { # BW!
    
    blabla <- suppressWarnings(plot.bw(chip1, chrom, min(mat1$x),
                                       max(mat1$x), y.values$chip1.y1,
                                       y.values$chip1.y2,
                                       col = col[1],
                                       rotate = rotate, yMax = yMax[1]
    ))
  }
  
  # - plot middle track c(0.5,0.8) -> chipOuter
  
  if (typeof(chip2) == "list") { # BED!
    
    if (ncol(chip2) < 6) {
      type[2] <- "rectangle"
      # there is no orientatio in column 6
    } else if (autoCHIP & !any(c("-", "+") %in% chip2[, 6])) {
      type[2] <- "rectangle"
    }
    
    if (type[2] == "triangle") {
      # +
      tmp_step <- (y.values$chip2.y2 - y.values$chip2.y1) / 2
      tmp_step <- tmp_step * 0.9
      blabla <- plot.triangle(chip2[chip2[, 6] == "+", ],
                              chrom = chrom,
                              y1 = y.values$chip2.y1,
                              y2 = y.values$chip2.y1 + tmp_step,
                              start = min(mat1$x), end = max(mat1$x), rotate = rotate
      )
      # -
      blabla <- plot.triangle(chip2[chip2[, 6] == "-", ],
                              chrom = chrom,
                              y1 = y.values$chip2.y2 - tmp_step,
                              y2 = y.values$chip2.y2, start = min(mat1$x),
                              end = max(mat1$x), rotate = rotate
      )
      # # other
      blabla <- plot.rectangle(chip2[chip2[, 6] != "-" & chip2[, 6] != "+", ],
                               chrom = chrom, y1 = y.values$chip2.y1,
                               y2 = y.values$chip2.y2, start = min(mat1$x),
                               end = max(mat1$x), col = col[2], rotate = rotate
      )
    } else if (type[2] == "rectangle") {
      blabla <- plot.rectangle(chip2,
                               chrom = chrom, y1 = y.values$chip2.y1,
                               y2 = y.values$chip2.y2, start = min(mat1$x),
                               end = max(mat1$x), col = col[2], rotate = rotate
      )
    }
  } else if (typeof(chip2) == "character") { # BW!
    
    blabla <- suppressWarnings(plot.bw(chip2, chrom, min(mat1$x), max(mat1$x),
                                       y.values$chip2.y1, y.values$chip2.y2,
                                       col = col[2], rotate = rotate,
                                       yMax = yMax[2]
    ))
  }
  
  
  # - plot outer track c(0.9,1.2) -> genes
  
  
  if (!is.null(genes)) {
    blabla <- plot.genes(genes, chrom, min(mat1$x), max(mat1$x),
                         y.pos = y.values$gene.pos, rotate = rotate
    )
  }
  
  
  # - if rotate, plot an empty canvas for quadrant 4
  
  if (rotate) {
    plot(0,
         type = "n", xlim = range(mat1$x), ylim = c(0, 1.3), axes = F,
         xlab = "", ylab = ""
    )
  }
}





plot.bw <- function(file, chrom, start, end, y1, y2, col,
                    yMax = NULL, rotate = F) {

  d <- import_bigwig(file, chrom, start, end)
  
  max.val <- max(d$y)
  
  if (!is.null(yMax)) {
    d[, y := pmin(y, yMax)]
    max.val <- yMax
  } else {
    message("chip.yMax not given for a .bw-track: yMax is ", max.val)
  }

  y.range <- abs(y1 - y2)
  y.val <- y.range * d$y / max.val
  d[, y := (y / max.val) * y.range]
  if (rotate) {
    polygon(d$y, d$x, col = col)
  } else {
    polygon(d$x, d$y, col = col)
  }
}

# overlay TAD positions with the hi-c data
draw.tads <- function(tads, chrom, tads.type = "lower", tads.colour = "#006837",
                      lwd = 2) {
  
  # check if tads is a list of dataframes or a data.frame
  if (inherits(tads, "list")) { # this seems to be a tad!
    if (all(unlist(lapply(tads, inherits, "data.frame")))) { # all DFs in list!
      tadList <- T
    } else {
      stop("tads is not a (list of) data.frame!")
    }
  } else if (inherits(tads, "data.frame")) {
    tmpList <- list() # cool. there is one df of tads. Put in a list
    tmpList[[1]] <- tads
    tads <- tmpList
    rm(tmpList)
  } else {
    stop("tads is not a (list of) data.frame!")
  }
  # if you survived this, you will now have a list of dfs!
  
  # get the tads colour, resize and type
  if (length(tads) != length(tads.type)) {
    tads.type <- rep(tads.type[1], length(tads))
  }
  
  if (length(tads) != length(tads.colour)) {
    tads.colour <- rep(tads.colour[1], length(tads))
  }
  
  
  for (listIDX in 1:length(tads)) {
    tad <- tads[[listIDX]]
    
    # tads has min. 6 cols : 1 and 4 are the chrom. 2 is outer-edge 5' anchor.
    # 6 is outer edge 3' anchor
    tad <- tad[tad[, 1] == chrom, ]
    
    if (tads.type[listIDX] == "both") {
      rect(tad[, 2], tad[, 2], tad[, 3], tad[, 3],
           border = tads.colour[listIDX],
           lwd = lwd
      )
    } else if (tads.type[listIDX] == "lower") {
      segments(tad[, 2], tad[, 2], tad[, 2], tad[, 3],
               col = tads.colour[listIDX],
               lwd = lwd
      )
      segments(tad[, 2], tad[, 3], tad[, 3], tad[, 3],
               col = tads.colour[listIDX],
               lwd = lwd
      )
    } else if (tads.type[listIDX] == "upper") {
      segments(tad[, 2], tad[, 2], tad[, 3], tad[, 2],
               col = tads.colour[listIDX],
               lwd = lwd
      )
      segments(tad[, 3], tad[, 2], tad[, 3], tad[, 3],
               col = tads.colour[listIDX],
               lwd = lwd
      )
    } else {
      stop("Wrong option for TAD plot type: upper, lower and both are allowed")
    }
  }
}

draw1loop <- function(radius, x.midpoint, y.midpoint, lty = 1, col = "black", lwd = 1) {
  x <- seq(x.midpoint - radius, x.midpoint + radius, 1)
  y <- seq(y.midpoint - radius, y.midpoint + radius, 1)
  
  curve((1 * (radius^2 - (x - x.midpoint)^2)^0.5 + y.midpoint),
        add = TRUE,
        from = (x.midpoint - radius), to = (x.midpoint + radius), lty = lty, col = col, lwd = lwd
  )
  curve((-1 * (radius^2 - (x - x.midpoint)^2)^0.5 + y.midpoint),
        add = TRUE,
        from = (x.midpoint - radius), to = (x.midpoint + radius), lty = lty, col = col, lwd = lwd
  )
}


draw.loops <- function(loops, chrom, start, end, radius = 1e5, col = "black", lwd = 1, lty = 1, type = "upper") {
  subLoop <- loops[loops[, 1] == chrom & loops[, 2] >= start & loops[, 6] <= end, 1:6]
  if (nrow(subLoop) == 0) {
    return()
  }
  subLoop$xmid <- apply(subLoop[, 2:3], 1, mean)
  subLoop$ymid <- apply(subLoop[, 5:6], 1, mean)
  
  # all to upper
  subLoop[subLoop$xmid < subLoop$ymid, ] <- setNames(subLoop[subLoop$xmid < subLoop$ymid, c(1:6, 8, 7)], colnames(subLoop))
  
  if (type == "lower") {
    subLoop <- setNames(subLoop[, c(1:6, 8, 7)], colnames(subLoop))
  } else if (type == "both") {
    subLoopLower <- setNames(subLoop[, c(1:6, 8, 7)], colnames(subLoop))
    subLoop <- rbind(subLoop, subLoopLower)
  }
  
  for (i in 1:nrow(subLoop)) {
    xmid <- subLoop[i, "xmid"]
    
    ymid <- subLoop[i, "ymid"]
    
    draw1loop(
      radius = radius, x.midpoint = xmid[1],
      y.midpoint = ymid[1], lty = lty, col = col, lwd = lwd
    )
  }
}




#' hic.matrixplot
#'
#' Plot a matrix (or two) for a region of interest with annotations
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp1 The Hi-C experiment object: produced by construct.experiment().
#' @param exp2 Optional: a Hi-C experiment object.
#' @param chrom Chromosome of the region of interest. Alternatively one of the 
#'   following: \itemize{
#'     \item{A 3-column, 1-row \code{data.frame} in BED-format.}
#'     \item{A single \code{character} describing a locus in UCSC notation,
#'        e.g. \code{"chr1:30,000,000-35,000,000"}}
#'   }
#'   The latter two options automatically provide the \code{start} and 
#'   \code{end} arguments too.
#' @param start,end Start and end position of the region of interest. Ignored if
#'   the \code{chrom} argument is a BED-format data.frame.
#' @param colour_lim,cut.off A \code{numeric} of length 2 for the lower and 
#'   upper colour limits of the hic-matrix plot. \code{cut.off}: older syntax.
#' @param chip A list of feature tracks, can be bed structure
#' (i.e. data frames) or a path to bigwig file (i.e. character variable),
#' maximum length of 4. Placement is inner-top, outer-top, outer-left,
#' inner-left.
#' @param inferno White/Red/black or White/Red coloscale?
#' @param cexTicks Change size of numbers on the axis
#' @param chip.colour A vector of four, parallel to [chip], of the to use colour.
#' @param chip.yMax A vector of four, parallel to [chip], of the maximum 
#' height of the biwigs. If only one value is given, all be set on this value.
#' @param type Should a rectangle or a triangle be drawn?
#' Note that for a triangle a 6th strand column should be included
#' @param guessType If an element in the chip is a dataframe,
#' infer the type from column 6? If true, column six should hold "+" and "-".
#' If a row has other characters, we view this entry as no-oriention and
#' thus plot a rectangle.
#' @param coplot When drawing together two experiments: [dual] is bottom
#' triangle exp2, top triangle exp1; [diff] plots a substraction of exp2-exp1
#' @param genes Structure with gene information,
#' will only be combined with bed structure
#' @param tads BED-like dataframe or a list of these data.frames
#' @param tads.type How to show TADS: upper, lower and or both
#' @param tads.colour Which colour do you want the TADs to have?
#' @param loops BED-like dataframe or a list of these data.frames
#' @param loops.type How to show loops: upper, lower and or both
#' @param loops.radius Set the size of the loop-circle to X bp. Can help 
#' visibility.
#' @param loops.colour Which colour do you want the loops to have?
#' @param skipAnn Do not plot outside-annotation. Can be used to plot other
#' things next to the matrix.
#' @param symmAnn Put features 1&2 also on verical (ignore chip-entries 3&4)
#' @param check.genome Check if reference genomes in exp1 and exp2 are the same
#' @param smoothNA Set to TRUE to perform a Nadaraya/Watson normalization. 
#' This will try to eliminate white stripes: this is only cosmetic and has no 
#' effect on the compartment-scores.
#' @param fillNAtreshold Set the amount strength of outlier correction for 
#' fillNA.
#' @param rasterise Set to true to use a bitmap raster instead of polygons.
#' @param colour_bar A \code{logical} of length 1, indicating whether a
#'   colour-bar legend should be drawn at the right.
#' @param addnames When the \code{coplot} argument is \code{"dual"}, display
#'   names of the samples? (default: \code{TRUE})
#' @param antoni Logical: plot an explorer of the microscopic world
#' @section Resolution recommendation: A resolution in the ballpark of 
#'   \code{(end - start) / 500}.
#' @note
#' To plot genes, a gene-model data.frame must be made. This can be done via a
#' multitude of ways (e.g. biomart, UCSC table browser). The resulting
#' data.frame must have one exon per row with the following columns:
#' geneID | transcriptID | chrom | txStart | txEnd | exonStart | exonEnd |
#' strand.
#' Alternatively, if you download the knowngene table from UCSC, you can
#' directly use this table (with exons combined per row), by renaming
#' exonStarts and exonEnds to exonStart and exonEnd.
#' @family matrix plots
#' @examples
#' \dontrun{
#' # plot two matrices of Hap-1 cells, including their respective loop-calls
#' hic.matrixplot(
#'   exp1 = Hap1_Haarhuis2017_10kb,
#'   exp2 = Hap1_Sanborn2015_10kb,
#'   chrom = "chr7",
#'   start = 25e6,
#'   end = 30e6,
#'   loops = list(Haarhuis2017_Loops, sanborn2015_Loops),
#'   loops.colour = c("blue", "green"),
#'   loops.type = c("upper", "lower"),
#'   loops.resize = c(20e3, 20e3), # expand for visibility
#'   colour_lim = c(0, 25)
#' ) # upper limit of contacts
#' }
#' @return A matrix-plot
#' @export
hic_matrixplot <- function(exp1, exp2 = NULL, chrom, start, end, 
                           colour_lim = NULL,
                           chip = list(NULL, NULL, NULL, NULL), inferno = NULL,
                           cexTicks = 1, chip.colour = "black", chip.yMax = NULL,
                           type = rep("triangle", 4), guessType = T,
                           coplot = "dual", genes = NULL,
                           tads = NULL, tads.type = "lower", loops = NULL,
                           loops.type = "lower", tads.colour = "#1faee3", # direct complementary colour of mid-orange inferno
                           loops.radius = NULL, loops.colour = "#1faee3",
                           skipAnn = F, symmAnn = F,
                           check.genome = T, smoothNA = F, fillNAtreshold = 2, 
                           rasterise = FALSE, addnames = TRUE, cut.off = NULL,
                           colour_bar = FALSE,
                           antoni = F) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if (is.null(loops.radius)) {
    loops.radius <- attr(exp1, "res") * 5
  }
  
  if (length(chip) < 3) {
    symmAnn <- T
  }
  
  if (!is.null(colour_lim)) {
    if (length(colour_lim) != 2 || !is.numeric(colour_lim)) {
      stop("The `colour_lim` argument should be a numeric of length 2.")
    }
  }
  
  loc <- standardise_location(chrom = chrom, start = start, end = end, 
                              singular = TRUE)
  
  if (!check_input()) {
    return(invisible())
  }
  
  # some error handling
  if (!is.null(exp2)) {
    if (any(attr(exp1, "rmChrom"), attr(exp2, "rmChrom"))) {
      check.genome <- F
    }
    # make sure the resolutions are the same
    if (attr(exp1, "res") != attr(exp2, "res")) {
      stop("The Hi-C matrices should have the same resolution")
    }
    
    if (check.genome) {
      if (!all(exp1$IDX[["V4"]] == exp2$IDX[["V4"]])) {
        msg <- paste0(
          "Not all ICE indexes are the same.\nAre you these ",
          "experiments were mapped to the same genome (-build)?"
        )
        stop(msg)
      }
    }
  }
  
  # if only one colour is given use it for all feature tracks
  if (length(chip.colour) == 1) {
    chip.colour <- rep(chip.colour, 4)
  }
  
  # fill up empty yMax-elements
  if (length(chip.yMax) < 4) {
    for (i in (length(chip.yMax) + 1):4) {
      chip.yMax[i] <- chip.yMax[1]
    }
  }
  
  # if only one colour is given use it for all feature tracks
  if (length(type) == 1) {
    type <- rep(type, 4)
  }
  
  # create a plotting layout
  w <- 6
  lay <- matrix(4, nrow = w, ncol = w)
  lay[2:w, 2:w] <- 1
  lay[1, ] <- 2
  lay[, 1] <- 3
  lay[1, 1] <- 4
  if (colour_bar) {
    lay <- cbind(lay, c(5, rep(6, nrow(lay) - 1)))
  }
  layout(lay)
  par(mar = rep(1, 4), xaxs = "i", yaxs = "i")
  # layout
  
  # get a matrix from the experiment
  mat1 <- select_subset(exp1, loc$chrom, loc$start, loc$end)
  if (!is.null(exp2)) {
    mat2 <- select_subset(exp2, loc$chrom, loc$start, loc$end)
  }
  
  
  
  if (smoothNA) {
    mat1 <- fillNAs(mat1, MADtreshold = fillNAtreshold)
    
    if (!is.null(exp2)) {
      mat2 <- fillNAs(mat2, MADtreshold = fillNAtreshold)
    }
  }
  
  
  ZnormScale <- attr(exp1, "znorm")
  if (!is.null(exp2)) {
    if (attr(exp1, "znorm") != attr(exp2, "znorm")) {
      stop("One experiment is Z-normalised and the other is not.")
    }
    
    ZnormScale <- attr(exp1, "znorm") | attr(exp2, "znorm")
  }
  
  ##############################################################################
  # plot AVL for fun
  ##############################################################################
  if (antoni) {
    load(url("https://github.com/robinweide/GENOVA/raw/dev/R/antoni.Rdat"))
    CR <- colorRampPalette(c(
      "#ffffff", "#f0f0f0", "#d9d9d9",
      "#bdbdbd", "#969696", "#737373", "#525252",
      "#252525", "#000000"
    ))
    
    image(antoni[, nrow(antoni):1 ],
          col = rev(CR(10)),
          axes = F, ylim = rev(range(antoni))
    )
    
    ##############################################################################
    # plot only exp2
    ##############################################################################
  } else if (is.null(exp2)) {
    
    # set cutoffs
    if (is.null(colour_lim)) {
      base <- if (ZnormScale) c(-1, 1) else c(0, 1)
      if (!is.null(cut.off)) {
        colour_lim <- base * cut.off
      } else {
        colour_lim <- base * max(quantile(abs(mat1$z), 0.99))
        message(
          "No colour limits were given: using 99% percentile: ",
          round(colour_lim[2]), "."
        )
      }
    }
    mat1$z[] <- pmax(pmin(mat1$z, colour_lim[2]), colour_lim[1])
    
    wr <- if (ZnormScale) {
      bezier_corrected_divergent
    } else if (is.null(inferno)) {
      .choose_palette()
    } else if (inferno) {
      bezier_corrected_hot
    } else {
      bezier_corrected_whiteRed
    }
    wr <- colorRampPalette(wr)
    image.default(mat1, col = wr(1e4), axes = FALSE, ylim = rev(range(mat1$x)),
                  zlim = colour_lim,
                  useRaster = literalTRUE(rasterise))
  }
  ##############################################################################
  # plot exp1 and exp2
  ##############################################################################
  else {
    if (coplot == "dual") {
      mat1$z[upper.tri(mat1$z)] <- mat2$z[upper.tri(mat2$z)]
      
      # set cutoffs
      if (is.null(colour_lim)) {
        base <- if (ZnormScale) c(-1, 1) else c(0, 1)
        if (!is.null(cut.off)) {
          colour_lim <- base * cut.off
        } else {
          colour_lim <- base * max(quantile(abs(mat1$z), 0.99))
          message(
            "No colour limits were given: using 99% percentile: ",
            round(colour_lim[2]), "."
          )
        }
      }
      mat1$z[] <- pmax(pmin(mat1$z, colour_lim[2]), colour_lim[1])
      
      wr <- if (ZnormScale) {
        bezier_corrected_divergent
      } else if (is.null(inferno)) {
        .choose_palette()
      } else if (inferno) {
        bezier_corrected_hot
      } else {
        bezier_corrected_whiteRed
      }
      wr <- colorRampPalette(wr)
      image.default(mat1, col = wr(1e4), axes = FALSE, ylim = rev(range(mat1$x)), 
            zlim = colour_lim,
            useRaster = literalTRUE(rasterise))
      
      expnames <- c(expnames(exp1), expnames(exp2))
      if (length(expnames) > 1) {
        txt <- seq(min(mat1$x), max(mat1$x), length.out = 20)[c(2, 19)]
        text(x = rev(txt)[1], y = txt[1], labels = expnames[1], adj = c(1, 1))
        text(x = rev(txt)[2], y = txt[2], labels = expnames[2], adj = c(0, 0))
      }
    } else {
      mat1$z <- mat2$z - mat1$z
      # set cutoffs
      if (is.null(colour_lim)) {
        base <- c(-1, 1)
        if (!is.null(cut.off)) {
          colour_lim <- base * cut.off
        } else {
          colour_lim <- base * max(quantile(abs(mat1$z), 0.99))
          message(
            "No colour limits were given: using 99% percentile: ",
            round(colour_lim[2]), "."
          )
        }
      }
      mat1$z[] <- pmax(pmin(mat1$z, colour_lim[2]), colour_lim[1])
      
      wr <- colorRampPalette(bezier_corrected_divergent)
      image.default(mat1, col = wr(500), axes = F, ylim = rev(range(mat1$x)), 
                    zlim = colour_lim)
    }
  }
  
  # draw pretty axes and boxes
  box(lwd = 2)
  size.region <- diff(range(mat1$x))
  if (size.region > 40e6) {
    axis(2,
         at = seq(0, 3e9, by = 10e6), labels = seq(0, 3e9, by = 10e6) / 1e6,
         lwd = 2, cex.axis = cexTicks
    )
    axis(3,
         at = seq(0, 3e9, by = 10e6), labels = seq(0, 3e9, by = 10e6) / 1e6,
         lwd = 2, cex.axis = cexTicks
    )
  } else if (size.region > 2e6) {
    axis(2,
         at = seq(0, 3e9, by = 1e6), labels = seq(0, 3e9, by = 1e6) / 1e6,
         lwd = 2, cex.axis = cexTicks
    )
    axis(3,
         at = seq(0, 3e9, by = 1e6), labels = seq(0, 3e9, by = 1e6) / 1e6,
         lwd = 2, cex.axis = cexTicks
    )
  } else {
    lab <- seq(0, 3e9, by = 500e3) / 1e6
    lab <- sprintf("%.1f", lab)
    axis(2, at = seq(0, 3e9, by = 500e3), labels = lab, lwd = 2, cex.axis = cexTicks)
    axis(3, at = seq(0, 3e9, by = 500e3), labels = lab, lwd = 2, cex.axis = cexTicks)
  }
  
  # draw tads on the image plot
  if (!is.null(tads)) {
    draw.tads(tads, loc$chrom, tads.type = tads.type, tads.colour = tads.colour)
  }
  
  # draw loops on the image plot
  if (!is.null(loops)) {
    draw.loops(loops, chrom = loc$chrom, start = loc$start, end = loc$end, 
               type = loops.type, radius = loops.radius, col = loops.colour, lwd = 2)
  }
  
  # fill up empty elements
  if (length(chip) < 4) {
    for (i in (length(chip) + 1):4) {
      chip[i] <- list(NULL)
    }
  }
  if (!skipAnn) {
    
    # plot the features horizontal
    features(mat1, loc$chrom, genes,
             chip1 = chip[[1]], chip2 = chip[[2]],
             autoCHIP = guessType, yMax = chip.yMax[1:2], col = chip.colour[1:2], type = type[1:2]
    )
    
    # clone 1 and two to feature entries 3 and 4
    if (symmAnn) {
      if (!is.null(chip[[1]])) {
        chip[[3]] <- chip[[1]]
      }
      if (!is.null(chip[[2]])) {
        chip[[4]] <- chip[[2]]
      }
      
      chip.yMax[3:4] <- chip.yMax[1:2]
      chip.colour[3:4] <- chip.colour[1:2]
    }
    
    features(mat1, loc$chrom, genes, chip[[3]], chip[[4]],
             autoCHIP = guessType, yMax = chip.yMax[3:4], col = chip.colour[3:4],
             type = type[3:4], rotate = T
    )
  }
  if (colour_bar) {
    plot.new()
    par(plt = c(0, 0.2, 0.03, 0.97))
    # If difference is extremely small, let it be at least something to 
    # prevent zero-range bug
    if (diff(colour_lim) < 1000 * .Machine$double.eps) {
      colour_lim <- c(-10, 10) * .Machine$double.eps
    }
    m <- t(as.matrix(seq(colour_lim[1], colour_lim[2], 
                         length.out = 255)))
    image(1, m[1,], m, col = wr(255), xaxt = "n", yaxt = "n",
          xlab = "", ylab = "")
    title <- "Contacts"
    title <- if (!is.null(exp2) && coplot != "dual") "Difference" else title
    title <- if (ZnormScale) "Z-score" else title
    axis(side = 4, lwd = 0, lwd.ticks = 1, lend = 1)
    mtext(title, side = 4, line = 2.5)
  }
  
}
