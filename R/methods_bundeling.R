# As a rule of thumb, for any discovery the following is ideally true:
# identical(discovery, do.call(bundle, unbundle(discovery)))

# Class Documentation -----------------------------------------------------------

#' @title Discovery class
#' @name discovery
#'
#' @description The discovery classes contains the results of the analysis
#'   functions within GENOVA. Each discovery object is a list with the data as 
#'   list-elements and each discovery object has attributes that may be relevant 
#'   to downstream functions.
#'
#' @details Next to the results, the \code{discovery} objects may also contain
#'   useful metadata on the analysis that was run, such as which resolution was
#'   used. For details on a particular \code{discovery} class, please consult 
#'   the documentation of the functions that generate the class.
#'
#'   Functions that generate \code{discovery} objects are the following:
#'   \describe{ \item{\code{\link[GENOVA]{PESCAn}}}{\code{PESCAn_discovery}
#'   objects} \item{\code{\link[GENOVA]{APA}}}{\code{APA_discovery} objects}
#'   \item{\code{\link[GENOVA]{ATA}}}{\code{ATA_discovery} objects}
#'   \item{\code{\link[GENOVA]{ARA}}}{\code{ARA_discovery} objects}
#'   \item{\code{\link[GENOVA]{RCP}}}{\code{RCP_discovery} objects}
#'   \item{\code{\link[GENOVA]{compartment_score}}}{\code{CS_discovery} objects}
#'   \item{\code{\link[GENOVA]{saddle}}}{\code{saddle_discovery} objects}
#'   \item{\code{\link[GENOVA]{insulation_score}}}{\code{IS_discovery} objects}
#'   \item{\code{\link[GENOVA]{insulation_domainogram}}}{\code{domainogram_discovery}
#'    objects}
#'   \item{\code{\link[GENOVA]{virtual_4C}}}{\code{virtual4C_discovery} objects}
#'   \item{\code{\link[GENOVA]{direct_index}}}{\code{DI_discovery} objects}
#'   \item{\code{\link[GENOVA]{intra_inter_TAD}}}{\code{IIT_discovery} objects}}
#'
#' @section Operations: \subsection{Subsetting}{\code{discovery} objects can be
#'   subsetted by using \code{subset(discovery, i)} wherein \code{i} is an
#'   \code{integer} or \code{character} corresponding to the intended
#'   sample(s).} \subsection{Combining}{\code{discovery} objects of the same
#'   type can be combined by using \code{\link[GENOVA]{bundle}(discovery_A,
#'   discovery_B)}. Generally, discovery objects generated from different
#'   resolutions or specific to a particular genomic region can not be
#'   combined.} \subsection{Splitting}{\code{discovery} objects can be split to
#'   individual samples by using \code{\link[GENOVA]{unbundle}(discovery)}.}
#'   \subsection{Visualisation}{\code{discovery} objects can be visualised with
#'   ggplot2 using \code{\link[GENOVA]{visualise}(discovery)}. Alternatively,
#'   the different discovery types also have base R plotting methods by using
#'   \code{plot(discovery)}.} \subsection{Quantification}{\code{discovery}
#'   objects can be quantified using \code{\link[GENOVA]{quantify}(discovery)}}
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
#' @export
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
bundle <- function(..., collapse = "_") {
  UseMethod("bundle", list(...)[[1]])
}

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
      dim(new) <- c(head(dims[1, ], -1), sum(dims[, ncol(dims)]))

      # Check dimnames
      dnames <- lapply(slot, dimnames)
      dlen <- lengths(dnames)[1]
      expnames <- do.call(c, lapply(dnames, `[[`, dlen))
      if (length(unique(expnames)) < length(expnames)) {
        expnames <- unlist(lapply(seq_len(nrow(dims)), function(j) {
          paste0(dnames[[j]][[dlen]], collapse, disco_names[[j]])
        }))
      }
      dnames <- dnames[[1]]
      # Set dimnames
      dimnames(new) <- c(head(dnames, -1), list(expnames))
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
  out <- discos[[1]]$scores
  if (length(discos) > 1L) {
    for(i in 2:length(discos)) {
      out <- merge(out, discos[[i]]$scores, by = c("window", "position"))
    }
  }
  
  expnames <- lapply(discos, expnames)
  cnames <- tail(colnames(out), -2)
  if (!identical(cnames, unlist(expnames))) {
    newnames <- lapply(seq_along(expnames), function(i) {
      paste0(expnames[[i]], collapse, i)
    })
    colnames(out)[-c(1:2)] <- unlist(newnames)
  }

  structure(
    list(scores = out),
    class = c("domainogram_discovery", "genomescore_discovery", "discovery"),
    package = attr(discos[[1]], "package"),
    chrom = chroms[[1]],
    resolution = res[[1]]
  )
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
  dfs <- lapply(dfs, as.data.table)
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
  }
  setkey(out, "chrom", "start")
  
  structure(list(insula_score = as.data.frame(out)),
            PACKAGE = "GENOVA",
            colours = cols,
            class = c("IS_discovery", "genomescore_discovery", "discovery"),
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
  vps <- lapply(vps, `colnames<-`, c("chr", "start", "end", "exp"))
  vps <- do.call(rbind, vps)
  
  if (length(unique(vps[, 1])) > 1) {
    stop("Can only bundle virtual 4Cs with the",
         "viewpoint on the same chromosome", call. = FALSE)
  }

  # vps <- vps[!duplicated(vps),]
  rownames(vps) <- NULL
  
  xlims <- unique(unlist(lapply(discos, attr, "xlim")))
  expnames <- unlist(lapply(discos, function(disc) {
    tail(colnames(disc$data), -2)
  }))

  datas <- lapply(discos, function(disc){disc$data})
  if (any(duplicated(expnames))) {
    vps <- lapply(seq_along(discos), function(i) {
      vp <- attr(discos[[i]], "viewpoint")
      vp$exp <- paste0(vp$exp, collapse, i)
      vp
    })
    vps <- do.call(rbind, vps)
  
    datas <- lapply(seq_along(datas), function(i) {
      x <- datas[[i]]
      colnames(x) <- c(head(colnames(x), 2),
                       paste0(tail(colnames(x), -2), collapse, i))
      x
    })
  }
  
  out <- datas[[1]]
  if (length(datas) > 1) {
    for (i in tail(seq_along(datas), -1)) {
      out <- merge(out, datas[[i]], by = c("chromosome", "mid"))
    }
  }
  
  structure(
    list(data = out), 
    class = c("virtual4C_discovery", "genomescore_discovery", "discovery"),
    'viewpoint' = vps, 
    'xlim' = xlims,
    'sample' = expnames,
    'resolution' = attr(discos[[1]], 'resolution'),
    package = "GENOVA"
  )
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
    dat <- as.data.table(disc$compart_scores)
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
  
  structure(list(compart_scores = as.data.frame(dat)),
            package = "GENOVA",
            colours = cols,
            class = c("CS_discovery", "genomescore_discovery", "discovery"),
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
            class = c("saddle_discovery", "discovery"))
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

#' @rdname bundle
#' @export
bundle.DI_discovery <- function(..., collapse = "_"){
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
  dfs <- lapply(discos, `[[`, "DI")
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
  }
  
  structure(list(DI = out),
            PACKAGE = "GENOVA",
            colours = cols,
            class = c("DI_discovery", "genomescore_discovery", "discovery"),
            resolution = attr(discos[[1]], "resolution"))
}


#' @rdname bundle
#' @export
bundle.IIT_discovery <- function(..., collapse = "_") {
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
    stop("Can only bundle discoveries of the same resolution.",
         call. = FALSE)
  }
  
  dat <- lapply(discos, `[[`, "results")
  
  # Check if TADs are the same
  tads <- discos[[1]]$tads
  checktads <- vapply(discos, function(disc) {
    identical(tads, disc$tads)
  }, logical(1L))
  if (!all(checktads)) {
    warning("Cannot couple observations from different TAD calls.",
            "Returning a plain data.frame with observations.")
    newdata <- lapply(dat, melt.data.table, id.vars = c("x", "y"))
    newdata <- rbindlist(newdata)
    return(as.data.frame(newdata))
  }
  
  expnames <- unlist(lapply(dat, function(df) {
    tail(colnames(df), -2)
  }))
  if (any(duplicated(expnames))) {
    dat <- lapply(seq_along(dat), function(i) {
      df <- dat[[i]]
      colnames(df)[-c(1,2)] <- paste0(colnames(df)[-c(1,2)], collapse, i)
      df
    })
  }
  
  newdat <- dat[[1]]
  for (i in tail(seq_along(dat), -1)) {
    newdat <- merge.data.table(newdat, dat[[i]], by = c("x", "y"))
  }
  setkeyv(newdat, NULL)
  
  cols <- unlist(lapply(discos, attr, "colours"))
  
  structure(
    list(results = newdat,
         tads = tads),
    class = c("IIT_discovery", "discovery"),
    package = attr(discos[[1]], "package"),
    colours = unname(cols),
    resolution = res
  )
}

#' @rdname bundle
#' @export
bundle.chrommat_discovery <- function(..., collapse = "_") {
  discos <- list(...)
  disco_names <- names(discos)
  if (is.null(disco_names)) {
    disco_names <- seq_along(discos)
  }
  
  # Perform checks
  .check_disco_classes(discos)
  .check_disco_resolutions(discos)
  mode  <- .check_disco_attr(discos, "mode")
  slots <- .check_disco_slots(discos)
  
  # Merge discos
  out <- lapply(slots, function(i) {
    slot <- lapply(discos, `[[`, i)
    # Check is array and not matrix
    if (inherits(slot[[1]], "array") && length(dim(slot[[1]])) != 2) {
      new <- .merge_array_slot(arrays = slot, 
                               exp_dim = 3L, 
                               names = disco_names, 
                               slotname = i, 
                               collapse = collapse)
    } else {
      stop("Unrecognised slot type: '", i, "'.", call. = FALSE)
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
#' @export
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
unbundle <- function(discovery, ...) {
  UseMethod("unbundle", discovery)
}

# Unbundle functions ------------------------------------------------------

#' @rdname unbundle
#' @export
unbundle.ARMLA_discovery <- function(discovery, ...) {
  slotnames <- names(discovery)

  # Find experiment names
  expnames <- lapply(slotnames, function(i) {
    thisclass <- class(discovery[[i]])
    if (thisclass == "array") {
      thesenames <- tail(dimnames(discovery[[i]]), 1)[[1]]
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
  expnames <- expnames(discovery)
  lapply(setNames(expnames, expnames), function(i) {
    col <- c("window", "position", i)
    out <- discovery
    out$scores <- out$scores[, col]
    attr(out, "resolution") <- attr(discovery, "resolution")
    attr(out, "chrom") <- attr(discovery, "chrom")
    attr(out, "package") <- attr(discovery, "package")
    out
  })
}

#' @rdname unbundle
#' @export
unbundle.IS_discovery <- function(discovery, ...) {
  exps <- tail(colnames(discovery$insula_score), -4)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chrom", "start", "end", "bin", i)
  })
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    structure(list(insula_score = discovery$insula_score[, cols[[i]]]),
              PACKAGE = "GENOVA",
              colours = attr(discovery, "colours")[i],
              class = c("IS_discovery", "genomescore_discovery", "discovery"),
              resolution = attr(discovery, "resolution"),
              window = attr(discovery, "window"))
  })
}

#' @rdname unbundle
#' @export
unbundle.virtual4C_discovery <- function(discovery, ...) {
  exps <- tail(colnames(discovery$data), -2)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chromosome", "mid", i)
  })
  vp <- attr(discovery, "viewpoint")
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    thisvp <- vp[vp$exp == tail(cols[[i]], 1),]
    rownames(thisvp) <- NULL
    structure(
      list(data = discovery$data[, cols[[i]]]),
      package = "GENOVA",
      colours = attr(discovery, "colours")[i],
      class = c("virtual4C_discovery", "genomescore_discovery", "discovery"),
      resolution = attr(discovery, "resolution"),
      viewpoint = thisvp
    )
  })
}


#' @rdname unbundle
#' @export
unbundle.CS_discovery <- function(discovery, ...) {
  exps <- expnames(discovery)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chrom", "start", "end", "bin", i)
  })
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    structure(list(compart_scores = discovery$compart_scores[, cols[[i]]]),
              PACKAGE = "GENOVA",
              colours = attr(discovery, "colours")[i],
              class = c("CS_discovery", "genomescore_discovery", "discovery"),
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
              class = c("saddle_discovery", "discovery"))
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
              class = c("RCP_discovery", "discovery"),
              package = "GENOVA",
              norm = nor)
  }, r = raw, s = smooth, SIMPLIFY = FALSE)
}

#' @rdname unbundle
#' @export
unbundle.DI_discovery <- function(discovery, ...) {
  exps <- tail(colnames(discovery$DI), -4)
  cols <- lapply(setNames(exps, exps), function(i) {
    c("chrom", "start", "end", "bin", i)
  })
  
  out <- lapply(setNames(seq_along(exps), exps), function(i) {
    structure(list(DI = discovery$DI[, cols[[i]]]),
              PACKAGE = "GENOVA",
              colours = attr(discovery, "colours")[i],
              class = c("DI_discovery", "genomescore_discovery", "discovery"),
              resolution = attr(discovery, "resolution"),
              window = attr(discovery, "window"))
  })
}

#' @rdname unbundle
#' @export
unbundle.IIT_discovery <- function(discovery, ...) {
  
  dat <- discovery$results
  expnames <- tail(colnames(dat), -2)
  colours <- attr(discovery, "colours")
  
  lapply(setNames(seq_along(expnames), expnames), function(i) {
    grab <- c("x", "y", expnames[i])
    structure(
      list(
        results = dat[, ..grab],
        tads = discovery$tads
      ),
      class = c("IIT_discovery", "discovery"), 
      package = attr(discovery, "package"),
      resolution = attr(discovery, "resolution"),
      colours = colours[i]
    )
  })
}

#' @rdname unbundle
#' @export
unbundle.chrommat_discovery <- function(discovery, ...) {
  expnames <- dimnames(discovery$obs)[[3]]
  
  lapply(setNames(seq_along(expnames), expnames), function(i) {
    structure(list(
      obs = discovery$obs[, , i, drop = FALSE],
      exp = discovery$exp[, , i, drop = FALSE]
    ), 
    class = c("chrommat_discovery", "discovery"), 
    mode = attr(discovery, "mode"),
    package = attr(discovery, "package"),
    resolution = attr(discovery, "resolution"))
  })
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
      if (length(dim(slot)) == 3) {
        return(slot[,, i, drop = FALSE])
      } else if (length(dim(slot)) == 4) {
        return(slot[,,, i, drop = FALSE])
      }
    } else if (thisclass == "list") {
      return(slot[i])
    }
  })
  extra_attr <- setdiff(names(attris), names(attributes(out)))
  attributes(out) <- c(attributes(out), attris[extra_attr])
  out
}

# TODO: refactor the (un)bundle methods to use these general 
# checkers/constructors to reduce code duplication, increase consistency and
# improve maintainability.

# Checks if a list of objects have the same first classes
.check_disco_classes <- function(discos) {
  classes <- vapply(discos, function(x){class(x)[1]}, character(1))
  if (multiclass <- length(unique(classes)) > 1) {
    warning(paste0("It is not recommended to bundle discoveries from", 
                   " different classes"), call. = FALSE)
  }
  invisible(multiclass)
}

# Checks if a list of objects have the same resolution attributes
.check_disco_resolutions <- function(discos) {
  res <- unique(vapply(discos, attr, numeric(1), "resolution"))
  if (length(res) > 1) {
    stop("Can only bundle discoveries of the same resolution.",
         call. = FALSE)
  } else {
    return(invisible(TRUE))
  }
}

# Checks if a list of objects have the same list element names
.check_disco_slots <- function(discos) {
  allslots <- unlist(lapply(discos, names))
  allslots <- table(allslots)[unique(allslots)]
  # slots are 'valid' if they are in every disco
  validslots <- names(allslots)[allslots == max(allslots)]
  if (length(validslots) == 0) {
    stop("No shared slots found between discovery objects", call. = FALSE)
  }
  if (!all(names(allslots) %in% validslots)) {
    dropslots <- setdiff(names(allslots), validslots)
    warning(paste0("'", paste0(dropslots, collapse = "', '"),
                   "' slots have been dropped."))
  }
  invisible(validslots)
}

# Used to check attributes are the same
.check_disco_attr <- function(discos, attr_name) {
  attrs <- lapply(discos, attr, attr_name)
  ref <- attrs[[1]]
  ident <- vapply(tail(attrs, -1), identical, logical(1), y = ref)
  if (any(!ident)) {
    stop(paste0("Attribute '", attr_name, "' not identical among discoveries"),
         call. = FALSE)
  }
  return(invisible(ref))
}

# Merges an array with equal dimensions, except for 'exp_dim' which can be 1:n
.merge_array_slot <- 
  function(arrays, 
           exp_dim = 3, 
           names = character(),
           slotname = character(),
           collapse = "_") {
    
    # Check dimensions
    dims <- lapply(arrays, dim)
    dims <- do.call(rbind, dims)
    uni_dims <- apply(dims, 2, function(x){length(unique(x))})
    if (length(table(uni_dims[-exp_dim])) > 1) {
      stop(paste0("Slot '", slotname, "' has incompatible dimensions,",
                  " have you ran the analysis at different resolutions or",
                  " with different size arguments?"),
           call. = FALSE)
    }
    
    # Formulate new array
    new <- do.call(c, arrays)
    new_dim <- dims[1,]
    new_dim[exp_dim] <- sum(dims[, exp_dim])
    dim(new) <- new_dim
    
    # Check dimnames
    dnames <- lapply(arrays, dimnames)
    expnames <- do.call(c, lapply(dnames, `[[`, exp_dim))
    # Resolve non-unique expnames
    if (length(unique(expnames)) < length(expnames)) {
      expnames <- mapply(function(first, second) {
        browser()
        paste0(first[[exp_dim]], collapse, second)
      }, first = dnames, second = names)
    }
    
    # Set dimnames
    dnames <- dnames[[1]]
    dnames[[exp_dim]] <- unname(expnames)
    dimnames(new) <- dnames
    
    return(new)
  }