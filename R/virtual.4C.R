#' Get an average virtual 4C of 2D-regions (e.g. loops).
#'
#' Extract matrices around a defined set of pixels, like possible loops.
#' Matrices are rescaled, so that all viewpoints and anchors are overlapping.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param loop.bed Bedpe file containing the loop positions: produced by read.bedpe().
#' @param smallTreshold The minimal size of loops. Too small loops lead to messy plots.
#' @param verbose Produces a progress-indication.
#' @return A vector containing the y-axis scores per rescaled bin.
virtual.4C <- function (experiment, loop.bed, smallTreshold = 225000, verbose = F)
{
    hicdata <- experiment$ICE
    bed <- experiment$ABS
    resolution <- experiment$RES
    if (is.null(data.table::key(hicdata))) {
        data.table::setkey(hicdata, V1, V2)
    }
    bed.p <- paste0(bed[, 1], ":", bed[, 2])
    loop.bed <- loop.bed[abs(loop.bed[, 6] - loop.bed[, 2]) >=
        smallTreshold, ]
    for (i in 1:length(loop.bed[, 1])) {
        if (loop.bed[i, 2] > loop.bed[i, 5]) {
            loop.bed <- loop.bed[i, c(4, 5, 6, 1, 2, 3)]
        }
    }
    loop.bed <- na.exclude(loop.bed[1:6])
    loop.4c <- loop.bed[, c(1, 2, 6)]
    loop.4c.length <- nrow(loop.4c)
    levels(loop.4c$V1) <- levels(bed$V1)
    results.vector <- rep(0, 100)
    SL <- 0
    for (i in 1:loop.4c.length) {
        if (verbose) {
            cat(paste0(i, " of ", loop.4c.length, " loops."),
                "\r")
        }
        loopSize <- loop.4c[i, 3] - loop.4c[i, 2]
        if (loopSize < smallTreshold) {
            next
        }
        halfLoopSize <- floor((loopSize/resolution)/1.5) * resolution
        viewpoint <- paste0(loop.4c[i, 1], ":", as.integer(resolution *
            floor(loop.4c[i, 2]/resolution)))
        start <- paste0(loop.4c[i, 1], ":", as.integer(resolution *
            floor((loop.4c[i, 2] - halfLoopSize)/resolution)))
        end <- paste0(loop.4c[i, 1], ":", as.integer(resolution *
            floor((loop.4c[i, 3] + halfLoopSize)/resolution)))
        viewpoint <- gsub(" ", "", viewpoint, fixed = TRUE)
        start <- gsub(" ", "", start, fixed = TRUE)
        end <- gsub(" ", "", end, fixed = TRUE)
        viewpoint.pos <- bed[match(viewpoint, bed.p), 4]
        start.pos <- bed[match(start, bed.p), 4]
        end.pos <- bed[match(end, bed.p), 4]
        if (!all(is.finite(viewpoint.pos), is.finite(start.pos),
            is.finite(end.pos))) {
            next
        }
        x <- rep(start.pos:end.pos, each = length(viewpoint.pos:viewpoint.pos))
        y <- rep(viewpoint.pos:viewpoint.pos, length(viewpoint.pos:viewpoint.pos))
        sel3 <- na.omit(hicdata[list(x, y)])
        sel5 <- na.omit(hicdata[list(y, x)])
        sel5 <- data.frame(sel5)[, c(2, 1, 3)]
        colnames(sel5) <- colnames(sel3)
        sel <- rbind(sel3, sel5)
        if (nrow(sel) < 10) {
            next
        }
        missingX <- data.frame(start.pos:end.pos)[!data.frame(start.pos:end.pos)[,
            1] %in% sel$V1, ]
        if (!is.null(dim(missingX))) {
            missingX.DF <- data.frame(V1 = missingX, V2 = rep(viewpoint.pos,
                length(missingX)), V3 = 0)
            sel <- rbind(sel, missingX.DF)
        }
        sel <- sel[order(sel$V1), ]$V3
        sel.resized <- approx(seq_along(sel), sel, n = 100)$y
        results.vector <- results.vector + sel.resized
        SL <- SL + 1
    }
    return(results.vector/SL)
}
