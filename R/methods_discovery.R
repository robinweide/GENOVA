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
#' # Running multiple analysis
#' ata1 <- ATA(WT_10kb, tads_wt)
#' ata2 <- ATA(KO_10kb, tads_ko)
#'
#' # Combining results
#' cata <- bundle(WT = ata1, KO = ata2)
#'
#' # Visualising the combined results
#' visualise(cata)
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

# Unbundle documentation --------------------------------------------------

#' @title Split discovery objects
#' @name unbundle
#'
#' @description \code{unbundle} takes a \code{discovery} object and splits these
#'   out to individual experiments.
#'
#' @param discovery A \code{discovery} object with more than 1 sample.
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
#' # Getting a discovery object
#' apa <- APA(list(WT_20kb, KO_20kb), loops)
#'
#' # Splitting the results
#' split <- unbundle(apa)
#'
#' # Plotting the first result only
#' visualise(split[[1]])
NULL

# Unbundle functions ------------------------------------------------------

#' @rdname unbundle
#' @export
unbundle.ARMLA_discovery <- function(discovery) {
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

# Utilities ---------------------------------------------------------------

#' @export
subset.ARMLA_discovery <- function(discovery, i) {
  oldclass <- class(discovery)[[1]]
  if (!("signal" %in% names(discovery))) {
    warning(paste0("Invalid ", oldclass, "object: no 'signal' array found.\n",
                   "Returning 'NULL'"),
            call. = FALSE)
    return(NULL)
  }

  attris <- attributes(discovery)
  class(discovery) <- "list"

  out <- lapply(discovery, function(slot) {
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

