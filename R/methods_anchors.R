# Documentation -----------------------------------------------------------

#' @name anchors
#' @title Anchors for Hi-C
#'
#' @description Anchors are indices to the matrix in a \code{contacts} object.
#'   Anchor functions translate genomic coordinates, typically BED-formatted
#'   \code{data.frame}s, into indices corresponding to locations in the Hi-C
#'   matrix. Anchors are used in aggregate analysis functions to indicate where
#'   parts of the matrix should be looked up. Anchors come in different types,
#'   depending on the aggregate analysis functions they are used for.
#'
#' @param IDX The indices slot of a GENOVA \code{contacts} object.
#' @param res The resolution attribute of a GENOVA \code{contacts} object.\emph{
#'   Not for ATA}.
#' @param bed,bedlist A BED-formatted \code{data.frame} with the following 3
#'   columns: \enumerate{ \item A \code{character} giving the chromosome names.
#'   \item An \code{integer} with start positions. \item An \code{integer} with
#'   end positions. } For the \code{bedlist} variant, a named \code{list} of the 
#'   above. (\emph{CSCAn only}).
#' @param bedpe A BEDPE-formatted \code{data.frame} with the following 6
#'   columns. \emph{APA only}: \enumerate{ \item A \code{character} giving the
#'   chromosome names of the first coordinate. \item An \code{integer} giving
#'   the start positions of the first coordinate. \item An \code{integer} giving
#'   the end positions of the first coordinate. \item A \code{character} giving
#'   the chromosome names of the second coordinate. \item An \code{integer}
#'   giving the start positions of the second coordinate. \item An
#'   \code{integer} giving the end positions of the second coordinate. }
#' @param mode A \code{character} vector of length 1 indicating which
#'   interactions to retain. Possible values: \code{"cis"}, \code{"trans"} or
#'   \code{"both"}. \emph{PE-SCAn, C-SCAn and APA only}.
#' @param dist_thres An \code{integer} vector of length 2 indicating the minimum
#'   and maximum distances in basepairs between anchorpoints. For ATA-type
#'   anchors, the minimum and maximum sizes of TADs.
#' @param min_compare An \code{integer} vector of length 1 indicating the
#'   minimum number of pairwise interactions on a chromosome to consider.
#'   \emph{PE-SCAn and C-SCAn only}.
#' @param padding A \code{numeric} of length 1 to determine the padding around
#'   TADs, expressed in TAD widths. \emph{ATA only}.
#' @param group_direction A \code{logical} of length 1 which when \code{TRUE} 
#'   will mirror groups for anchors where the left anchor location is larger 
#'   than the right anchor location. Left and right refer to bedlist elements 
#'   generating combinations. \emph{CSCAn only}.
#' @param strand A \code{character} of the length \code{nrow(bed)}. Overrules
#'   an attempt to infer strand from \code{start > end} information. 
#'   \emph{ARA only}.
#'
#' @return A \code{anchors} object with two colums in \code{matrix} format.
#'
#' @details Anchors are calculated within aggregate repeated matrix lookup
#'   analysis, but can also be provided as the '\code{anchors}' argument for
#'   these functions.
#'
#'   Anchors are specific for a resolution of a \code{contacts} object and
#'   cannot be interchanged freely between resolutions.
#'
#'   The '\code{mode}' argument determines what pairwise interactions are
#'   reported for \code{APA} and \code{PE-SCAn}. \code{"cis"} returns pairwise
#'   interactions within a chromosome; \code{"trans"} gives these between
#'   chromosomes and \code{"both"} returns both these types of interactions.
#'
#' @section Anchor types:
#'
#'   \subsection{PE-SCAn anchors}{ \code{anchors_PESCAn()} takes all pairwise
#'   interactions of genomic coordinates from a BED-like \code{data.frame} and
#'   maps these to indices of the Hi-C matrix. It is used within the
#'   \code{\link[GENOVA]{PESCAn}} function. Wether these pairwise interactions
#'   are allowed to cross chromosome boundaries is determined by the
#'   '\code{mode}' argument, which default to \code{"cis"} to only take pairwise
#'   interactions on the same chromosome.}
#'
#'   \subsection{C-SCAn anchors}{ \code{anchors_CSCan()}, like
#'   \code{anchors_PESCAn}, takes all pairwise interactions of genomic
#'   coordinates, but crosswise between unique combinations of BED-like
#'   \code{data.frame}s in the \code{bedlist} argument. It is used within the
#'   \code{\link[GENOVA]{CSCAn}} function. Has a \code{group} attribute to
#'   keep track from which combination of BED-like \code{data.frame} the anchor
#'   originated.}
#'
#'   \subsection{APA anchors}{ \code{anchors_APA()} takes a BEDPE-formatted
#'   \code{data.frame} and translates the coordinates in the first 3 and last 3
#'   columns to indices of the Hi-C matrix. It is used within the
#'   \code{\link[GENOVA]{APA}} function. The '\code{mode}' argument defaults to
#'   \code{"both"} but optionally allows for either cis- or trans-interactions
#'   too. }
#'
#'   \subsection{Extended loops anchors}{ \code{anchors_extendedloops()} takes
#'   the same input as \code{anchors_APA()}, but transforms these coordinates to
#'   combinations of 5' and 3' anchors outside existing loops to get 'extended'
#'   loops. Based on the extended loops algorithm described in Haarhuis \emph{et
#'   al.} (2017).}
#'
#'   \subsection{ATA anchors}{ \code{anchors_ATA()} takes the genomic
#'   coordinates of TADs in a BED-formatted \code{data.frame} and translates to
#'   indices of the Hi-C matrix. It is used within the \code{\link[GENOVA]{ATA}}
#'   function. In contrast to the PE-SCAn anchors and ATA anchors, ATA anchors
#'   are positions on the matrix's diagonal. The '\code{padding}' argument
#'   controls how large the region around a TAD should be expanded. Since TADs
#'   have variable sizes, ATA anchors can be calculated without resolution.}
#'
#'   \subsection{ARA anchors}{ \code{anchors_ARA()} takes a BED-formatted
#'   \code{data.frame} and translates these to Hi-C matrix indices on the
#'   diagonal. It is used within the \code{\link[GENOVA]{ARA} function}. In
#'   contrast to other anchors, ARA anchors can take on a directionality. If the
#'   start positions are larger than the end positions, the anchor is assigned a
#'   \emph{reverse} direction. Else they are given the \emph{forward}
#'   direction.}
#'
#' @seealso \code{\link[GENOVA]{bed2idx}} for general genomic coordinates to
#'   Hi-C index conversion.
#'
#' @examples
#' \dontrun{
#' # PE-SCAn
#' anch <- anchors_PESCAn(WT_20kb$IDX, attr(WT_20kb, "resolution"),
#'                        super_enhancers)
#' PESCAn(list(WT_20kb, KO_20kb), anchors = anch)
#'
#' # APA
#' anch <- anchors_APA(WT_20kb$IDX, attr(WT_20kb, "resolution"), loops)
#' APA(list(WT_20kb, KO_20kb), anchors = anch)
#'
#' # APA with extended loops
#' ex_anch <- anchors_extendedloops(WT_20kb$IDX, attr(WT_20kb, "resolution"),
#'                                  loops)
#' APA(list(WT_20kb, KO_20kb), anchors = ex_anch)
#'
#' # ATA
#' anch <- anchors_ATA(WT_10kb$IDX, tads)
#' ATA(list(WT_10kb, KO_10kb), anchors = anch)
#'
#' # ARA
#' anch <- anchors_ARA(WT_20kb$IDX, ctcf_sites)
#' ARA(list(WT_20kb, KO_20kb), anchors = anch)
#' }
NULL

# Types --------------------------------------------------------------

#' @rdname anchors
#' @export
anchors_PESCAn <- function(IDX, res, bed,
                           dist_thres = c(5e6L, Inf),
                           min_compare = 10L,
                           mode = c("cis", "trans", "both")) {
  mode <- match.arg(mode)

  # Match and clean indices
  newbed <- cbind(bed, idx = bed2idx(IDX, bed, mode = "centre"))
  newbed <- newbed[!is.na(newbed$idx), ]
  newbed <- newbed[order(newbed$idx), ]

  proxy <- seq_len(nrow(newbed))
  if (mode == "cis") {

    # Generate cis combinations
    # idx <- split(newbed$idx, newbed[, 1])
    idx <- split(proxy, newbed[, 1])
    idx <- idx[lengths(idx) >= max(min_compare, 2)]
    if (length(idx) == 0) {
      stop("Too few comparable positions per chromosome.", call. = FALSE)
    }
    idx <- lapply(idx, function(x) t(utils::combn(x, 2)))
    idx <- do.call(rbind, idx)
    is_cis <- TRUE
  } else {

    # Generate all combinations
    # idx <- t(utils::combn(newbed$idx, 2))
    idx <- t(utils::combn(proxy, 2))
    is_cis <- newbed[c(idx), 1]
    # is_cis <- newbed[match(idx, newbed$idx), 1]
    is_cis <- matrix(is_cis, ncol = 2)
    is_cis <- is_cis[, 1] == is_cis[, 2]

    if (mode == "trans") {

      # Exclude cis when trans
      idx <- idx[!is_cis, ]
      is_cis <- FALSE
    }
  }
  
  faux_bedpe <- cbind(newbed[idx[, 1], 1:3], newbed[idx[, 2], 1:3])
  rownames(idx) <- names_from_bedpe(faux_bedpe)
  # Replace proxy by true indices
  idx[] <- newbed[c(idx), "idx"]
  

  # Filtering distances only makes sense in cis
  if (mode != "trans") {
    # Covert dists to resolution space
    dist_thres <- dist_thres / res
    dist_thres <- sort(dist_thres)

    # Exclude cis combinations based on min/max_dist
    cis <- idx[is_cis, , drop = FALSE]
    dist <- abs(cis[, 1] - cis[, 2])
    keep <- dist >= dist_thres[1] & dist <= dist_thres[2]
    if (sum(keep) < 1) {
      stop("No pairwise interactions are within the distance threshold",
           call. = FALSE)
    }
    cis <- cis[keep, , drop = FALSE]

    # Recombine and order
    idx <- rbind(idx[!is_cis, , drop = FALSE], cis)
    idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  }

  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "PESCAn"
  idx
}

#' @export
#' @rdname anchors
anchors_CSCAn <- function(IDX, res, bedlist,
                          dist_thres = c(50e3, 2e6),
                          min_compare = 10L,
                          mode = c("cis", "trans", "both"),
                          group_direction = FALSE) {
  mode <- match.arg(mode)
  if (length(bedlist) < 2 || !inherits(bedlist, "list")) {
    stop("Less than two 'bedlist' elements found. For self-interaction of a",
         " single BED-like data.frame, see '?anchors_PESCAn'.",
         call. = FALSE)
  }
  
  chroms <- lapply(lapply(bedlist, `[[`, 1), unique)
  chroms <- table(unlist(chroms))
  chroms <- names(chroms)[chroms >= 2]

  if (length(chroms) == 0) {
    stop("No common chromosomes found between 'bedlist' argument elements",
         call. = FALSE)
  }
  
  # Match and clean indices
  if (is.null(names(bedlist))) {
    names(bedlist) <- LETTERS[seq_along(bedlist)]
  }
  if (length(unique(names(bedlist))) < length(bedlist)) {
    stop("Please provide unique names for every element in the 'bedlist' argument",
         call. = FALSE)
  }
  
  beds <- lapply(bedlist, function(x) {
    x <- x[!duplicated(x),]
    newbed <- cbind.data.frame(chrom = x[, 1], 
                               idx = bed2idx(IDX, x, mode = "centre"))
    rownames(newbed) <- names_from_bed(x)
    newbed <- newbed[!is.na(newbed$idx), ]
    newbed <- newbed[order(newbed$idx),]
  })
  
  cbn <- as.data.frame(combn(names(beds), 2), stringsAsFactors = FALSE)
  if (mode == "cis") {
    # Generate cis combinations
    splitbeds <- lapply(beds, function(x) {split(seq_len(nrow(x)), x$chrom)[chroms]})
    lens <- lapply(splitbeds, lengths)

    idx <- lapply(cbn, function(i) {
      left <- i[[1]]; right <- i[[2]]
      keep <- (lens[[left]] * lens[[right]]) >= max(min_compare, 2)
      proto_idx <- rbindlist(
        mapply(CJ, 
               splitbeds[[left]][keep], 
               splitbeds[[right]][keep], 
               SIMPLIFY = FALSE)
      )[, combi := paste0(i[1], "-", i[2])]
      proto_idx[, rnames := paste0(
        rownames(beds[[left]])[V1], ";",
        rownames(beds[[right]])[V2]
      )]
      # Substitute proxy for real idx
      proto_idx[, V1 := beds[[left]][V1, "idx"]]
      proto_idx[, V2 := beds[[right]][V2, "idx"]]
      return(proto_idx)
    })
    idx <- rbindlist(idx)
    is_cis <- TRUE
  } else {
    ind <- lapply(lapply(beds, nrow), seq_len)
    idx <- rbindlist(lapply(cbn, function(i) {
      left <- i[[1]]; right <- i[[2]]
      # Generate all combinations
      out <- CJ(ind[[left]], ind[[right]])
      # Check cis/trans
      out[, is_cis := beds[[left]][[1]][V1] == beds[[right]][[1]][V2]]
      # Attach rownames
      out[, rnames := paste0(rownames(beds[[left]])[V1], ";",
                             rownames(beds[[right]])[V2])]
      
      # Substitute actual indices, add combination column
      out[, c("V1", "V2", "combi") := list(beds[[i[1]]][[2]][V1],
                                           beds[[i[2]]][[2]][V2],
                                           paste0(i[1], "-", i[2]))]
    }))
    if (mode == "trans") {
      # Exclude cis when trans
      idx <- idx[is_cis == FALSE, ]
      is_cis <- FALSE
    } else {
      is_cis <- idx[["is_cis"]]
    }
    idx[, is_cis := NULL]
    # idx <- as.data.frame(idx)
  }
  
  if (group_direction) {
    str <- lapply(strsplit(idx$combi, ""), rev)
    str <- as.data.frame(do.call(rbind, str))
    str <- do.call(paste0, str)
    flip <- idx[, V2 < V1]
    idx[flip, combi := str[flip]]
  }
  
  # Sort start-end
  idx[, c("V1", "V2") := list(pmin(V1, V2), pmax(V1, V2))]
  
  # Filtering distances only makes sense in cis
  if (mode != "trans") {
    # Covert dists to resolution space
    dist_thres <- dist_thres / res
    dist_thres <- sort(dist_thres)
    
    # Exclude cis combinations based on min/max_dist
    cis <- idx[is_cis, , drop = FALSE]
    dist <- cis[, abs(V1 - V2)]
    keep <- dist >= dist_thres[1] & dist <= dist_thres[2]
    if (sum(keep) < 1) {
      stop("No pairwise interactions are within the distance threshold",
           call. = FALSE)
    }
    cis <- cis[keep]
    
    # Recombine and order
    idx <- rbind(idx[!is_cis], cis)
  }
  rnames <- idx$rnames
  combi <- rle(idx[["combi"]])
  idx[, combi := NULL]
  idx[, rnames := NULL]
  idx <- as.matrix(idx)
  dimnames(idx) <- NULL
  rownames(idx) <- rnames
  
  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "CSCAn"
  attr(idx, "group") <- combi
  idx
}

#' @rdname anchors
#' @export
anchors_APA <- function(IDX, res, bedpe,
                        dist_thres = c(0, Inf),
                        mode = c("both", "cis", "trans")) {
  mode <- match.arg(mode)

  posnames <- names_from_bedpe(bedpe)
  
  # Convert bedpe to idx matrix
  newbed <- cbind(bedpe[, 1:6],
    idx1 = bed2idx(IDX, bedpe[, 1:3], mode = "centre"),
    idx2 = bed2idx(IDX, bedpe[, 4:6], mode = "centre")
  )
  rownames(newbed) <- posnames
  newbed <- stats::na.exclude(newbed)
  newbed <- newbed[!duplicated(newbed[, 7:8]), ]

  # Shortcut when no additional filtering is needed
  if (identical(dist_thres, c(0, Inf)) && mode == "both") {
    # Ordering
    idx <- cbind(
      pmin(newbed$idx1, newbed$idx2),
      pmax(newbed$idx1, newbed$idx2)
    )
    rownames(idx) <- rownames(newbed)
    idx <- idx[order(idx[, 1]), , drop = FALSE]
    class(idx) <- c("anchors", "matrix")
    attr(idx, "type") <- "APA"
    return(idx)
  }

  newbed$is_cis <- newbed[, 1] == newbed[, 4]

  if (mode == "trans") {
    newbed <- newbed[!newbed$is_cis, ]
  } else {
    # Convert distances to bins
    dist_thres <- dist_thres / res
    dist_thres <- sort(dist_thres)

    # Filter cis on distances
    cis <- newbed[newbed$is_cis, ]
    dist <- abs(cis$idx1 - cis$idx2)
    cis <- cis[dist >= dist_thres[1] & dist <= dist_thres[2], ]
    if (mode == "cis") {
      newbed <- cis
    } else {
      newbed <- rbind(newbed[!newbed$is_cis, ], cis)
    }
  }

  # Ordering
  idx <- cbind(
    pmin(newbed$idx1, newbed$idx2),
    pmax(newbed$idx1, newbed$idx2)
  )
  rownames(idx) <- rownames(newbed)
  idx <- idx[order(idx[, 1]), , drop = FALSE]
  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "APA"
  return(idx)
}

#' @rdname anchors
#' @export
anchors_ATA <- function(IDX, bed,
                        dist_thres = c(225000, Inf),
                        padding = 1) {
  if (!inherits(bed, "data.frame")) {
    bed <- as.data.frame(bed)[, 1:3]
  }

  # Setup parameters
  width <- abs((bed[, 3] - bed[, 2]))
  mid <- round((bed[, 2] + bed[, 3]) / 2)
  keep <- width  > dist_thres[1] & width < dist_thres[2]
  if (sum(keep) < 1) {
    stop("There are no TADs large enough to pass the distance thresholds.",
         call. = FALSE)
  }
  rnames <- names_from_bed(bed[keep,])

  # Resize regions
  bed <- data.frame(bed[,1],
                    mid - width * padding,
                    mid + width * padding)[keep, ]

  # Translate to Hi-C indices
  idx <- cbind(
    bed2idx(IDX, bed, mode = "start"),
    bed2idx(IDX, bed, mode = "end")
  )
  # Sort
  idx <- cbind(
    pmin(idx[, 1], idx[, 2]),
    pmax(idx[, 1], idx[, 2])
  )
  rownames(idx) <- rnames
  idx <- idx[order(idx[, 1]), , drop = FALSE]
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  
  # Attribute to let matrix lookup methods know it is performing ATA
  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "TADs"
  attr(idx, "padding") <- padding
  return(idx)
}

#' @rdname anchors
#' @export
anchors_ARA <- function(IDX, bed, strand = NULL) {
  if (!inherits(bed, "data.frame")) {
    bed <- as.data.frame(bed)[, 1:3]
  }
  # Translate to indices
  idx <- bed2idx(IDX, bed)
  is_dup <- duplicated(idx)
  idx <- idx[!is_dup]
  idx <- unname(cbind(idx, idx))
  rownames(idx) <- names_from_bed(bed[!is_dup,])

  # Attach direction if necessary
  if (is.null(strand)) {
    f <- rle(ifelse(bed[!is_dup, 2] < bed[!is_dup, 3], "+", "-"))
  } else {
    if (length(strand) != nrow(bed)) {
      stop(paste0("Attempting to calculate stranded ARA anchors, but strand ",
                  "information is not of the same length as BED rows."))
    }
    f <- rle(ifelse(strand[!is_dup] == "-", "-", "+"))
  }

  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "ARA"
  attr(idx, "dir") <- f
  return(idx)
}

#' @rdname anchors
#' @export
anchors_extendedloops <- function(IDX, res, bedpe,
                                  dist_thres = c(3e4, 3e6)) {
  # Convert distances to bins
  dist_thres <- dist_thres / res
  dist_thres <- sort(dist_thres)
  
  # Filter for cis
  bedpe <- bedpe[bedpe[, 1] == bedpe[, 4], ]
  
  posnames <- paste0(bedpe[, 1], ":", 
                     format(bedpe[, 2], scientific = FALSE, trim = TRUE), "-", 
                     format(bedpe[, 3], scientific = FALSE, trim = TRUE), ";",
                     bedpe[, 4], ":", 
                     format(bedpe[, 5], scientific = FALSE, trim = TRUE), "-", 
                     format(bedpe[, 6], scientific = FALSE, trim = TRUE))
  
  # Convert bedpe to idx matrix
  newbed <- data.table(chrom = bedpe[, 1],
                       idx1  = bed2idx(IDX, bedpe[, 1:3], mode = "centre"),
                       idx2  = bed2idx(IDX, bedpe[, 4:6], mode = "centre")
  )
  
  # Order stuff
  newbed <- stats::na.exclude(newbed)
  newbed[, c("idx1", "idx2") := list(pmin(idx1, idx2), pmax(idx1, idx2))]
  newbed <- newbed[!duplicated(newbed), ]
  newbed <- newbed[order(idx1, idx2),]
  
  # Filter very large loops
  newbed <- newbed[abs(idx1 - idx2) < dist_thres[2]]
  
  # Find unique 5' anchors within minimum distance
  newbed[, uni5p := c(1, diff(idx1) > dist_thres[1]), by = chrom]
  newbed[, uni5p := cumsum(uni5p)]
  newbed[, row := seq_len(nrow(newbed))]

  # Determine cumulative maximum of 3' anchors
  cmax <- newbed[, list(max = max(idx2)), by = c("chrom", "uni5p")][["max"]]
  
  # Make all combinations of loop anchors
  combi <- newbed[, CJ(i = row[!duplicated(uni5p)], j = row), by = chrom]
  
  # Filter out unwanted loops from the same 5' anchors
  combi <- combi[newbed[i, uni5p] < newbed[j, uni5p], list(i, j)]
  
  # Translate back to actual indices
  combi <- combi[, list(idx1 = newbed[i, idx1], idx2 = newbed[j, idx2], i)]
  
  # Filter out loops larger than maximum distance
  combi <- combi[idx1 + dist_thres[2] > idx2]
  
  # Filter out loops smaller than cumulative maximum of 3' anchors
  combi <- combi[idx2 > cmax[newbed[i, uni5p]] + dist_thres[1]]
  
  # Exclude duplicates
  idx <- combi[!duplicated(combi), list(idx1, idx2)]
  idx <- as.matrix(idx)
  dimnames(idx) <- NULL
  class(idx) <- c("anchors", "matrix")
  attr(idx, "type") <- "APA"
  return(idx)
}

# Utilities ---------------------------------------------------------------

#' Finish anchors for repeated matrix lookup
#'
#' @inheritParams anchors_shift
#'
#' @return A \code{anchors} object of the same type.
#' @examples
#' \dontrun{
#' anch <- anchors_APA(WT_20kb$ABS, WT_20kb$RES, loops)
#' anchors_finish(WT_20kb$ABS, anch, -10:10, 0)
#' }
anchors_finish <- function(IDX, anchors, rel_pos, shift = 0) {
  anchors <- anchors_filter_oob(IDX, anchors, rel_pos)
  anch_id <- parse(text = paste0("seq.int(", nrow(anchors), ")"),
                   keep.source = FALSE)

  # Quick workaround for no shift
  if (shift == 0) {
    shft_id <- expression(0)
    attr(anchors, "anch_id") <- anch_id
    attr(anchors, "shft_id") <- shft_id
    return(anchors)
  }

  # Shift anchors as necessary
  shift <- anchors_shift(IDX, anchors, rel_pos, shift)
  shft_id <- parse(text = paste0("seq.int(", nrow(shift), ") + ",
                                 nrow(anchors)), keep.source = FALSE)

  # Take care of direction attribute
  attris <- attributes(anchors)
  if ("dir" %in% names(attris)) {
    attris$dir$lengths <- c(attris$dir$lengths, attr(shift, "dir")$lengths)
    attris$dir$values  <- c(attris$dir$values,  attr(shift, "dir")$values)
  }
  if ("group" %in% names(attris)) {
    attris$group$lengths <- c(attris$group$lengths, 
                              attr(shift, "group")$lengths)
    attris$group$values <- c(attris$group$values,
                             attr(shift, "group")$values)
  }

  # Make final anchors
  fin <- rbind(anchors, shift)

  # Re-attach attributes
  extra <- setdiff(names(attris), names(attributes(fin)))
  attributes(fin) <- c(attributes(fin), attris[extra])
  attr(fin, "anch_id") <- anch_id
  attr(fin, "shft_id") <- shft_id
  fin
}

#' Shift anchors
#'
#' Shifts anchors upstream by a specified amount, unless they become out of
#' bounds. In that case, shift the anchors downstream.
#'
#' @inheritParams rep_mat_lookup
#' @inheritParams anchors_PESCAn
#' @param shift An \code{integer} of length 1 indicating how many bins the
#'   anchors should be shifted.
#'
#' @return A \code{anchors} object with two columns in \code{matrix} format.
#'
#' @details The resulting matrix contains indices to the Hi-C matrix slot in the
#'   GENOVA experiment.
#'
#'   An index is considered out of bounds when that index plus the shift size
#'   and maximum relative position would belong to a different chromosome.
#'
#' @seealso \code{\link[GENOVA]{PESCAn}} for context and
#'   \code{\link[GENOVA]{anchors}}. \code{\link[GENOVA]{anchors_filter_oob}} for
#'   general out of bounds filtering of anchors.
#'
#' @export
anchors_shift <- function(IDX, anchors, rel_pos, shift = 1) {

  # Translate indices to chromosomes
  chrom   <- IDX[match(anchors, IDX[, V4]), V1]
  shifted <- IDX[match(anchors + shift + max(rel_pos), IDX[, V4]), V1]

  # Check chromosome is same after shift
  inbounds <- matrix(chrom == shifted, ncol = 2)
  inbounds <- apply(inbounds, 1, all)

  # Shift upstream unless out of bounds
  shifted <- ifelse(inbounds, 1, -1) * shift

  anchors + c(shifted, shifted)
}

#' Filter anchors that become out of bounds.
#'
#' Discards anchors for which a lookup will result in any lookup bins that are
#' out of bounds.
#'
#' @inheritParams anchors_shift
#'
#' @return A \code{anchors} object with two columns in \code{matrix} format.
#'
#' @details The resulting matrix contains indices to the Hi-C matrix slot in the
#'   GENOVA experiment.
#'
#'   An index is considered out of bounds when that index plus the maximum, or
#'   minus the minimum, relative position would belong to a different
#'   chromosome.
#'
#' @seealso \code{\link[GENOVA]{anchors}}.
#'
#' @export
anchors_filter_oob <- function(IDX, anchors, rel_pos) {
  type <- attr(anchors, "type")
  class(IDX) <- "data.frame"

  # ATA has slightly different oob rules
  if (type == "TADs") {
    left  <- IDX[match(anchors[, 1], IDX[, 4]), 1]
    right <- IDX[match(anchors[, 2], IDX[, 4]), 1]
    keep <- left == right
    anchors <- anchors[keep, ]
    attr(anchors, "type") <- type
    return(anchors)
  }

  # Match idx +/- relative position to chrom
  plus  <- pmax(match(anchors + max(rel_pos), IDX[, 4]), 1, na.rm = TRUE)
  minus <- pmax(match(anchors + min(rel_pos), IDX[, 4]), 1, na.rm = TRUE)
  plus  <- IDX[plus, 1]
  minus <- IDX[minus, 1]

  # Check wether chromosomes have changed
  inbounds <- matrix(plus == minus, ncol = 2)
  inbounds <- apply(inbounds, 1, all)
  if (sum(inbounds) < 1) {
    stop("No suitable anchors left after out-of-bounds filtering.",
         call. = FALSE)
  }

  # Return anchors that are not out of bounds
  anchors <- anchors[inbounds, , drop = FALSE]
  anchors
}


#' Coerce to anchors
#'
#' Function to coerce to \code{anchors} if possible.
#'
#' @param x A two column \code{matrix} or \code{data.frame} with integers.
#'
#' @return An object of class \code{anchors}.
#'
#' @details The resulting \code{anchors} has the default '\code{type}' attribute
#'   \code{"custom"}.
#' @export
#'
#' @seealso \code{\link[GENOVA]{anchors}}
#'
#' @examples
#' \dontrun{
#' as_anchors(matrix(1:20, 2))
#' }
as_anchors <- function(x) {
  if (is.null(dim(x)) || length(dim(x)) > 2) {
    stop(paste0("An object of class ", class(x),
                " cannot be converted to anchors.\nUse a ",
                "'matrix' or 'data.frame' object with 2 integer ",
                "columns instead."),
         call. = FALSE)
  }
  if (dim(x)[2] != 2) {
    stop(paste0("Anchors require 2 columns, not ", dim(x)[2], "."),
         call. = FALSE)
  }
  x <- as.matrix(x)
  if (storage.mode(x) != "integer") {
    stop(paste0("Cannot convert ", storage.mode(x),
                " type to anchors. Use integers instead."),
         call. = FALSE)
  }
  class(x) <- c("anchors", "matrix")
  attr(x, "type") <- "custom"
  x
}

#' Test object is anchors
#'
#' Function to check if an object is \code{anchors}.
#'
#' @param x Any \R object.
#'
#' @return A \code{logical} vector of length 1.
#' @export
#'
#' @seealso \code{\link[GENOVA]{anchors}}
#'
#' @examples
#' \dontrun{
#' m <- matrix(1:20, 10)
#' is_anchors(m) # FALSE
#'
#' m <- as_anchors(m)
#' is_anchors(m) # TRUE
#' }
is_anchors <- function(x) {
  inherits(x, "anchors")
}

# Need subsetting function to allow attribute inheritance.
#' @export
#' @keywords internal
`[.anchors` <- function(x, i, j, ..., drop = FALSE) {
  # Treat as matrix
  y <- x
  class(y) <- "matrix"
  y <- y[i, j, ..., drop = drop]
  if (is.null(dim(y))) {
    return(y)
  }
  # Ensure resulting object inherits attributes from parent
  attrs <- setdiff(names(attributes(x)),
                   names(attributes(y)))
  x_attr <- attributes(x)[attrs]
  # Take particular care that direction is subsetted as rows
  if ("dir" %in% names(x_attr)) {
    x_attr[["dir"]] <- rle(inverse.rle(attr(x, "dir"))[i])
  }
  if ("group" %in% names(x_attr)) {
    x_attr[["group"]] <- rle(inverse.rle(attr(x, "group"))[i])
  }
  attributes(y) <- c(attributes(y), x_attr)
  y
}

names_from_bedpe <- function(bedpe) {
  paste0(names_from_bed(bedpe[, 1:3]), ";",
         names_from_bed(bedpe[, 4:6]))
}

names_from_bed <- function(bed) {
  paste0(bed[,1], ":",
         format(bed[, 2], scientific = FALSE, trim = TRUE), "-",
         format(bed[, 3], scientific = FALSE, trim = TRUE))
}
