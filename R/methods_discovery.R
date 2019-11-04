# Class Documentation -----------------------------------------------------------

#' @title Discovery class
#' @name discovery
#'
#' @description The discovery class contains the results of the analysis
#'   functions within GENOVA.
#'
#' @details Next to the results, the \code{discovery} objects may also contain
#'   useful metadata on the analysis that was run, such as which resolution was
#'   used.
#'
#'   Functions that generate \code{discovery} objects are the following:
#'   \describe{
#'    \item{\code{\link[GENOVA]{PESCAn}}}{\code{PESCAn_discovery} objects}
#'     \item{\code{\link[GENOVA]{APA}}}{\code{APA_discovery} objects}
#'     \item{\code{\link[GENOVA]{ATA}}}{\code{ATA_discovery} objects}
#'     \item{\code{\link[GENOVA]{ARA}}}{\code{ARA_discovery} objects}
#'     \item{\code{\link[GENOVA]{RCP}}}{\code{RCP_discovery} objects}
#'     \item{\code{\link[GENOVA]{compartment_score}}}{\code{CS_discovery} objects}
#'     \item{\code{\link[GENOVA]{saddle}}}{\code{saddle_discovery} objects}
#'     \item{\code{\link[GENOVA]{insulation_score}}}{\code{IS_discovery} objects}
#'     \item{\code{\link[GENOVA]{insulation_domainogram}}}{\code{domainogram_discovery} objects}
#'     \item{\code{\link[GENOVA]{virtual_4C}}}{\code{virtual4C_discovery} objects}
#'   }
#'
#' @section Operations: \subsection{Subsetting}{\code{discovery} objects can be
#'   subsetted by using \code{subset(discovery, i)} wherein \code{i} is an
#'   \code{integer} or \code{character} corresponding to the intended
#'   sample(s).} \subsection{Combining}{\code{discovery} objects of the same
#'   type can be combined by using \code{\link[GENOVA]{bundle}(discovery_A,
#'   discovery_B)}.} \subsection{Splitting}{\code{discovery} objects can be
#'   split to individual samples by using
#'   \code{\link[GENOVA]{unbundle}(discovery)}.}
#'   \subsection{Visualisation}{\code{discovery} objects can be visualised using
#'   \code{\link[GENOVA]{visualise}(discovery)}}
#'   \subsection{Quantification}{\code{discovery} objects can be quantified
#'   using \code{\link[GENOVA]{quantify}(discovery)}}
NULL

# Bundle documentation ----------------------------------------------------

#' @title Combine discovery objects
#' @name bundle
#'
#' @description \code{bundle} takes one or more GENOVA \code{discovery} objects
#'   and combines them along experiments.
#'
#' @param ... A number of \code{discovery} objects from the same analysis type.
#' @param collapse An optional character string to seperate names in case the
#'   experiment names are not unique. Not \code{NA_character_}.
#'
#' @details This function can be convenient when a the same analysis has been
#'   run on a sequence of experiments with parameters specific to the
#'   experiment. For example, consider the case when \code{ATAs} have been run
#'   on several single experiments with TAD boundary anchors that match that
#'   same experiment. The \code{bundle} function allows combining the resulting
#'   \code{discovery} objects for visualisation or quantification.
#'
#'   The arrays in slots of the \code{discovery} objects need to have compatible
#'   first and second dimensions, i.e. they must have had the same
#'   '\code{size_bp}' or '\code{size_bin}' arguments and have been run at the
#'   same resolutions.
#'
#' @seealso The \code{\link[GENOVA]{discovery}} class.
#'   \code{\link[GENOVA]{unbundle}} for splitting \code{discovery} objects.
#'
#' @return A \code{discovery}-class object of the same type.
#'
#' @examples
#' \dontrun{
#' # Running multiple analysis
#' ata1 <- ATA(WT_10kb, tads_wt)
#' ata2 <- ATA(KO_10kb, tads_ko)
#'
#' # Combining results
#' cata <- bundle(WT = ata1, KO = ata2)
#'
#' # Visualising the combined results
#' visualise(cata)
#' }
NULL

# Bundle functions --------------------------------------------------------

#' @rdname bundle
#' @export
bundle.ARMLA_discovery <- function(..., collapse = "_") {
  # Yeah groovy, baby!
  discos <- list(...)
  disco_names <- names(discos)
  if (is.null(disco_names)) {
    disco_names <- seq_along(discos)
  }

  # Check classes
  classes <- vapply(discos, function(x){class(x)[1]}, character(1))
  if (length(unique(classes)) > 1) {
    warning(paste0("It is not recommended to bundle discoveries from",
                   "different classes"), call. = FALSE)
  }

  # Check slots
  allslots <- table(sapply(discos, names))
  slots <- names(allslots)[allslots == max(allslots)]
  if (length(slots) == 0) {
    stop("No shared slots found between discovery objects", call. = FALSE)
  }
  if (!all(names(allslots) %in% slots)) {
    dropslots <- setdiff(names(allslots), slots)
    warning(paste0("'", paste0(dropslots, collapse = "', '"),
                   "' slots have been dropped"))
  }

  # Merge
  out <- lapply(slots, function(i) {
    slot <- lapply(discos, `[[`, i)
    if (class(slot[[1]]) == "array") {

      # Check dimensions
      dims <- lapply(slot, dim)
      dims <- do.call(rbind, dims)
      if (length(table(dims[,1:2])) > 1) {
        stop(paste0("Slot '", i, "' has incompatible dimensions,",
                    " have you ran the analysis at different resolutions or",
                    " with different size arguments?"),
             call. = FALSE)
      }
      new <- do.call(c, slot)
      dim(new) <- c(dims[1,1:2], sum(dims[, 3]))

      # Check dimnames
      dnames <- lapply(slot, dimnames)
      expnames <- do.call(c, lapply(dnames, `[[`, 3))
      if (length(unique(expnames)) < length(expnames)) {
        expnames <- unlist(lapply(seq_len(nrow(dims)), function(j) {
          paste0(dnames[[j]][[3]], collapse, disco_names[[j]])
        }))
      }

      # Set dimnames
      dimnames(new) <- list(dnames[[1]][[1]], dnames[[1]][[2]], expnames)
    } else if (class(slot[[1]]) == "list") {

      # Resolve potential naming conflicts
      dnames <- lapply(slot, names)
      expnames <- unlist(dnames)
      if (length(unique(expnames)) < length(expnames)) {
        expnames <- unlist(lapply(seq_along(slot), function(j) {
          paste0(dnames[[j]], collapse, disco_names[[j]])
        }))
      }
      new <- do.call(c, slot)
      names(new) <- expnames
    } else {
      stop("Slot '", i, "' is of an unrecognised class",
           call. = FALSE)
    }
    return(new)
  })

  # Set attributes
  names(out) <- slots
  extra_attr <- setdiff(names(attributes(discos[[1]])),
                        names(attributes(out)))
  attributes(out) <- c(attributes(discos[[1]])[extra_attr],
                       attributes(out))

  out
}

#' @rdname bundle
#' @export
bundle.domainogram_discovery <- function(..., collapse = "_"){
  discos <- list(...)
  
  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  chroms <- vapply(discos, attr, character(1), "chrom")
  if (lengths(unique(chroms)) > 1) {
    stop("Can only bundle domainograms at the same chromosome.", call. = FALSE)
  }
  res <- vapply(discos, attr, numeric(1), "resolution")
  if (length(unique(res)) > 1) {
    stop("Can only bundle domainograms with the same resolution.", 
         call. = FALSE)
  }
  
  # Combine and reorder
  out <- do.call(rbind, discos)
  out <- out[order(out$window, out$position, out$experiment), ]
  
  # Filter out duplicates
  dups <- duplicated(out[, c("window", "position", "experiment")])
  if (sum(dups) > 0) {
    message("Found duplicated insulation scores which are discarded.")
  }
  out <- out[!dups,]
  out
}

#' @rdname bundle
#' @export
bundle.IS_discovery <- function(..., collapse = "_"){
  discos <- list(...)

  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  res <- vapply(discos, attr, numeric(1), "resolution")
  if (length(unique(res)) > 1) {
    stop("Can only bundle insulation scores of the same resolution.",
         call. = FALSE)
  }
  # Grab colours
  cols <- lapply(discos, attr, "colours")
  cols <- unlist(cols)
  
  # Extract insulation scores
  dfs <- lapply(discos, `[[`, "insula_score")
  expnames <- lapply(lapply(dfs, colnames), tail, -4)
  
  # Merge insulation scores
  out <- dfs[[1]]
  if (length(dfs) > 1) {
    for (i in tail(seq_along(dfs), -1)) {
      out <- merge(out, dfs[[2]], by = c("chrom", "start", "end", "bin"))
    }
  }
  
  # Check output column names
  cnames <- tail(colnames(out), -4)
  if (!identical(cnames, unlist(expnames))) {
    newnames <- lapply(seq_along(expnames), function(i) {
      paste0(expnames[[i]], collapse, i)
    })
    colnames(out)[-c(1:4)] <- unlist(newnames)
    return(out)
  }
  setkey(out, "chrom", "start")
  
  structure(list(insula_score = out),
            PACKAGE = "GENOVA",
            colours = cols,
            class = "IS_discovery",
            resolution = attr(discos[[1]], "resolution"),
            window = attr(discos[[1]], "window"))
}

#' @rdname bundle
#' @export
bundle.virtual4C_discovery <- function(..., collapse = "_") {
  discos <- list(...)
  
  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  res <- vapply(discos, attr, numeric(1), "resolution")
  if (length(unique(res)) > 1) {
    stop("Can only bundle insulation scores of the same resolution.",
         call. = FALSE)
  }
  
  vps <- lapply(discos, attr, "viewpoint")
  vps <- do.call(rbind, vps)
  
  if (length(unique(vps[, 1])) > 1) {
    stop("Can only bundle virtual 4Cs with the",
         "viewpoint on the same chromosome", call. = FALSE)
  }

  vps <- vps[!duplicated(vps),]
  rownames(vps) <- NULL
  
  xlims <- unique(unlist(lapply(discos, attr, "xlim")))
  expnames <- unname(vapply(discos, 
                            function(disc){unique(disc$data$experiment)},
                     character(1)))
  
  datas <- lapply(discos, function(disc){disc$data})
  datas <- rbindlist(datas)
  datas <- datas[order(chromosome, mid)]
  
  structure(list(data = datas), 
            class = "virtual4C_discovery",
            'viewpoint' = vps, 
            'xlim' = xlims,
            'sample' = expnames,
            'resolution' = attr(explist[[1]], 'resolution'),
            package = "GENOVA")
}

#' @rdname bundle
#' @export
bundle.CS_discovery <- function(..., collapse = "_") {
  discos <- list(...)
  if (length(discos) < 2) {
    message("Attempting to bundle a single object. Input is returned.")
    return(discos[[1]])
  }
  
  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  res <- vapply(discos, attr, numeric(1), "resolution")
  if (length(unique(res)) > 1) {
    stop("Can only bundle compartment scores of the same resolution.",
         call. = FALSE)
  }
  signage <- vapply(discos, attr, logical(1), "signed")
  if (length(unique(signage)) > 1) {
    message("Attempting to bundle compartment scores with mixed signed status.", 
            " Setting signed status of output to `FALSE`.")
    signage <- FALSE
  } else {
    signage <- unique(signage)
  }
  
  dats <- lapply(discos, function(disc) {
    dat <- disc$compart_scores
    setkeyv(dat, c("chrom", "start"))
    dat[, party := inverse.rle(attr(disc, "partitioning"))]
    dat
  })
  dat <- dats[[1]]
  # I'm updating the object, so I can use a for-loop here, ok?
  for (i in tail(seq_along(dats), -1)) {
    dat <- merge(dat, dats[[i]], by = c("chrom", "start", "end", "bin"),
                 suffixes = c("_1", paste0("_", i)))
  }
  
  # Deal with partitioning
  parts <- dat[, startsWith(colnames(dat), "party"), with = FALSE]
  party <- parts[[1]]
  for (i in tail(seq_len(ncol(parts)), -1)) {
    if (all(parts[[i]] == party)) {
      next()
    }
    # Prioritise centromeres over p/q arm
    party <- ifelse(is.na(party) | endsWith(parts[[i]], "centro"), parts[[i]],
                    party)
  }
  party <- rle(party)
  
  dat <- dat[, !startsWith(colnames(dat), "party"), with = FALSE]
  setkey(dat, bin)
  
  # Grab colours
  cols <- lapply(discos, attr, "colours")
  cols <- unname(unlist(cols))
  
  structure(list(compart_scores = dat),
            package = "GENOVA",
            colours = cols,
            class = "CS_discovery",
            resolution = unique(res),
            partitioning = party,
            signed = signage)
}

#' @rdname bundle
#' @export
bundle.saddle_discovery <- function(..., collapse = "_") {
  discos <- list(...)
  if (length(discos) < 2) {
    message("Attempting to bundle a single object. Input is returned.")
    return(discos[[1]])
  }
  
  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  
  res <- unique(vapply(discos, attr, numeric(1), "resolution"))
  if (length(res) > 1) {
    warning("Attempting to bundle saddle discoveries called from different",
            " resolutions.", call. = FALSE)
    res <- max(res)
  }
  
  ranges <- t(vapply(discos, function(disco) {
    disco$saddle[, range(c(q1, q2), na.rm = TRUE)]
  }, integer(2)))
  if (nrow(ranges[!duplicated(ranges), , drop = FALSE]) > 1) {
    stop("Can only bundle saddle discoveries with the same number of bins.",
         call. = FALSE)
  }
  
  dats <- lapply(discos, `[[`, "saddle")
  dats <- rbindlist(dats)
  
  structure(list(saddle = dats),
            package = "GENOVA",
            resolution = res,
            class = "saddle_discovery")
}

#' @rdname bundle
#' @export
bundle.RCP_discovery <- function(..., collapse = "_") {
  discos <- list(...)
  if (length(discos) < 2) {
    message("Attempting to bundle a single object. Input is returned.")
    return(discos[[1]])
  }
  
  # Check for possible errors
  classes <- vapply(lapply(discos, class), `[`, character(1), 1)
  if (length(unique(classes)) > 1) {
    stop("Can only bundle discoveries of the same type.", call. = FALSE)
  }
  
  norms <- vapply(discos, attr, character(1), "norm")
  norms <- unique(unname(norms))
  
  if (lengths(norms) > 1) {
    warning("Attempting to bundle RCPs with different normalisations.")
  }
  
  raws <- lapply(discos, `[[`, "raw")
  raws <- rbindlist(raws)
  setkeyv(raws, "distance")
  
  smooths <- lapply(discos, `[[`, "smooth")
  smooths <- rbindlist(smooths)
  
  structure(list(raw = raws, smooth = smooths),
            class = "RCP_discovery",
            package = "GENOVA",
            norm = norms[1])
}

# Unbundle documentation --------------------------------------------------

#' @title Split discovery objects
#' @name unbundle
#'
#' @description \code{unbundle} takes a \code{discovery} object and splits these
#'   out to individual experiments.
#'
#' @param discovery A \code{discovery} object with more than 1 sample.
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return A \code{list} wherein each element is a \code{discovery} object for a
#'   single sample.
#'
#' @details In case the \code{discovery} contains incomplete samples with
#'   missing slots, \code{NULL} is returned.
#'
#' @seealso \code{\link[GENOVA]{bundle}} for the opposite of \code{unbundle}.
#'   The \code{\link[GENOVA]{discovery}} class.
#'
#' @examples
#' \dontrun{
#' # Getting a discovery object
#' apa <- APA(list(WT_20kb, KO_20kb), loops)
#'
#' # Splitting the results
#' split <- unbundle(apa)
#'
#' # Plotting the first result only
#' visualise(split[[1]])
#' }
NULL

# Unbundle functions ------------------------------------------------------

#' @rdname unbundle
#' @export
unbundle.ARMLA_discovery <- function(discovery, ...) {
  slotnames <- names(discovery)

  # Find experiment names
  expnames <- lapply(slotnames, function(i) {
    thisclass <- class(discovery[[i]])
    if (thisclass == "array") {
      thesenames <- dimnames(discovery[[i]])[[3]]
    } else if (thisclass == "list") {
      thesenames <- names(discovery[[i]])
    }
    thesenames
  })
  lens <- lengths(expnames)
  if (length(unique(lens)) > 1) {
    warning("Not all experiments were found in all slots.")
  }
  expnames <- expnames[[which.max(lens)]]
  expnames <- setNames(expnames, expnames)

  out <- lapply(expnames, function(ename) {
    subset(discovery, ename)
  })
}

#' @rdname unbundle
#' @export
unbundle.domainogram_discovery <- function(discovery, ...) {
  split(discovery, discovery$experiment)
}

#' @rdname unbundle
#' @export
unbundle.IS_discovery <- function(discovery, ...) {
  exps <- tail(colnames(discovery$insula_score), -4)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chrom", "start", "end", "bin", i)
  })
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    structure(list(insula_score = discovery$insula_score[, cols[[i]], 
                                                         with = FALSE]),
              PACKAGE = "GENOVA",
              colours = attr(discovery, "colours")[i],
              class = "IS_discovery",
              resolution = attr(discovery, "resolution"),
              window = attr(discovery, "window"))
  })
}

#' @rdname unbundle
#' @export
unbundle.virtual4C_discovery <- function(discovery, ...) {
  attris <- attributes(discovery)
  newdata <- split(discovery$data, discovery$data$experiment)
  lapply(setNames(seq_along(newdata), names(newdata)), function(i) {
    structure(list(data = newdata[[i]]), class = "virtual4C_discovery",
              xlim = attris$xlim,
              viewpoint = attris$viewpoint,
              sample = attris$sample[i],
              resolution = attris$resolution,
              package = attris$package)
  })
}

#' @rdname unbundle
#' @export
unbundle.CS_discovery <- function(discovery, ...) {
  exps <- tail(colnames(discovery$compart_scores), -4)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chrom", "start", "end", "bin", i)
  })
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    structure(list(compart_scores = discovery$compart_scores[, cols[[i]], 
                                                             with = FALSE]),
              PACKAGE = "GENOVA",
              colours = attr(discovery, "colours")[i],
              class = "CS_discovery",
              resolution = attr(discovery, "resolution"),
              signed = attr(discovery, "signed"),
              partitioning = attr(discovery, "partitioning"))
  })
}

#' @rdname unbundle
#' @export
unbundle.saddle_discovery <- function(discovery, ...) {
  dats <- split(discovery$saddle, discovery$saddle$exp)
  lapply(dats, function(dat) {
    structure(list(saddle = dat),
              package = "GENOVA",
              resolution = attr(discovery, "resolution"),
              class = "saddle_discovery")
  })
}

#' @rdname unbundle
#' @export
unbundle.RCP_discovery <- function(discovery, ...) {
  raw <- discovery$raw
  smooth <- discovery$smooth
  
  raw <- split(raw, raw$samplename)
  smooth <- split(smooth, smooth$samplename)
  
  nor <- attr(rcp, "norm")
  
  if (length(raw) != length(smooth)) {
    stop("Different number of samples found in the raw and smooth data.")
  }
  
  mapply(function(r, s) {
    structure(list(raw = r, smooth = s),
              class = "RCP_discovery",
              package = "GENOVA",
              norm = nor)
  }, r = raw, s = smooth, SIMPLIFY = FALSE)
}

# Utilities ---------------------------------------------------------------

#' @export
subset.ARMLA_discovery <- function(x, i, ...) {
  
  oldclass <- class(x)[[1]]
  if (!("signal" %in% names(x))) {
    warning(paste0("Invalid ", oldclass, "object: no 'signal' array found.\n",
                   "Returning 'NULL'"),
            call. = FALSE)
    return(NULL)
  }

  attris <- attributes(x)
  class(x) <- "list"

  out <- lapply(x, function(slot) {
    thisclass <- class(slot)
    if (thisclass == "array") {
      return(slot[,,i, drop = FALSE])
    } else if (thisclass == "list") {
      return(slot[i])
    }
  })
  extra_attr <- setdiff(names(attris), names(attributes(out)))
  attributes(out) <- c(attributes(out), attris[extra_attr])
  out
}

