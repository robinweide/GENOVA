
#' @rdname load_contacts
#' @export
#' @param location A \code{character} vector of locations. Either chromosome
#'   names, e.g. \code{"chr3"}, or chromosomal locations e.g. 
#'   \code{"chr3:30000000-90000000"}
load_contacts_subset <- function(
    signal_path,
    indices_path = NULL,
    location = NULL,
    resolution = 10e3,
    sample_name = NULL,
    centromeres = NULL,
    colour = NULL,
    z_norm = FALSE,
    scale_bp = NULL,
    balancing = TRUE,
    verbose = TRUE
) {
  # Check location
  if (is.null(location)) {
    stop("Cannot read location = NULL subset of data.")
  }
  location <- as.character(location)
  if (any(grepl(":", location))) {
    # Assumed to be chromosomal coordinates
    location <- standardise_location(location)
  } else {
    # Assumed to be whole chromosomes
    location <- standardise_location(location, start = -Inf, end = Inf)
  }
  
  # Control data.table threads
  dt.cores <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(dt.cores))
  data.table::setDTthreads(1)
  
  if (is.null(sample_name)) {
    sample_name = basename(signal_path)
  }
  
  do_juicer <- FALSE
  do_cooler <- FALSE
  do_hicpro <- FALSE
  balanced  <- FALSE
  
  juicer_path <- NULL
  cooler_path <- NULL
  
  input_type = switch(
    tools::file_ext(tolower(signal_path)),
    'matrix' = 'hicpro',
    'cooler' = 'cooler',
    'cool'   = 'cooler',
    'mcool'  = 'cooler',
    'hic'    = 'juicer'
  )
  
  if (input_type == "juicer") {
    juicer_path <- signal_path
  }
  if (input_type == "cooler") {
    cooler_path <- signal_path
  }
  
  if (!file.exists(signal_path)) {
    stop(paste0("Cannot find file `", signal_path, "`."))
  }
  
  # Ensure appropriate stuff can be loaded
  if (!is.null(juicer_path)) {
    try_require('strawr', "load_contacts_subset", source = 'github')
    do_juicer <- TRUE
  } else if (!is.null(cooler_path)) {
    try_require('rhdf5', "load_contacts_subset", source = "BIOC")
    do_cooler <- TRUE
  } else if (!is.null(signal_path)) {
    if (file.exists(indices_path)) {
      stop("Cannot find fine `", indices_path, "`.")
    }
    do_hicpro <- TRUE
  } else {
    stop(
      paste0("Please provide either a .matrix/indices_path, ",
             "a .cooler or a .hic file.")
    )
  }
  
  if (verbose) {
    message("Reading data...")
  }
  
  if (do_juicer) {
    stop("Juicer not implemented yet.")
  }
  
  if (do_hicpro) {
    stop("HiC Pro not implemented yet.")
  }
  
  if (do_cooler) {
    ans <- load_cooler_subset(
      cooler_path, location,
      balancing = balancing,
      scale_bp = scale_bp,
      resolution = resolution
    )
    balanced = balancing
  }

  signal <- ans$SIG
  if (!is.data.table(signal)) {
    stop("Signal is not a data.table.")
  }
  if (is.null(key(signal))) {
    setkeyv(signal, c("V1", "V2"))
  }
  
  index <- ans$ABS
  if (!is.data.table(index)) {
    stop("Index is not a data.table")
  }
  if (is.null(key(index))) {
    setkeyv(index, c("V1", "V2"))
  }
  res <- as.numeric(median(index$V3 - index$V2))
  
  # Check bins
  signal_range <- signal[, range(V1, V2)]
  index_range  <- index[, range(V4)]
  
  if (signal_range[2] > index_range[2]) {
    stop("Undefined HiC-bins")
  }
  if (abs(signal_range[2] - index_range[2]) > signal_range[2] * 0.1) {
    stop("Undefined HiC-bins")
  }
  
  if (z_norm) {
    if (verbose) {message("Running Z-score normalisation...")}
    signal <- zscore_hic(signal, index)
  }
  
  if (!all(signal[, V1 <= V2])) {
    signal[, c("V1", "V2") := list(pmin(V1, V2), pmax(V1, V2))]
    setkeyv(signal, c("V1", "V2"))
  }
  
  if (!is.null(centromeres)) {
    if (inherits(centromeres, "data.frame")) {
      centromeres <- clean_centromeres(centromeres, index)
    } else if (
      is.logical(centromeres) && length(centromeres) == 1 && !centromeres
    ) {
      centromeres <- data.table(
        chrom = unique(index[[1]]),
        start = -1,
        end   = -1
      )
    } else {
      stop(paste0(
        "Don't know how to interpret `", typeof(centromeres), "` as centromeres"
      ))
    }
  } else {
    centromeres <- find_centromeres(signal, index)
  }
  
  # Construct output
  if (is.null(colour)) {
    colour <- "black"
  }
  
  structure(
    list(
      MAT = signal,
      IDX = index,
      CHRS = unique(index$V1),
      CENTROMERES = centromeres
    ),
    class = "contacts",
    znorm = z_norm,
    samplename = sample_name,
    colour = colour,
    resolution = res,
    rmChrom = TRUE,
    balanced = balanced,
    scale_cis = FALSE,
    package = "GENOVA"
  )
}


load_cooler_subset <- function(
    file,
    location  = NULL,
    balancing = TRUE,
    scale_bp  = NULL,
    resolution = 10e3
) {
  
  if (is.null(location)) {
    stop("Cannot find NULL location in file.")
  }
  
  bins    <- "bins"
  pixels  <- "pixels"
  indexes <- "indexes" 
  
  ls <- rhdf5::h5ls(file)
  is_mcool <- any(grepl("resolutions", ls$group))
  
  if (is_mcool) {
    
    resos <- tstrsplit(ls$group[grepl("resolutions", ls$group)], "/")
    resos <- unique(resos[[3]])
    resos <- resos[!is.na(resos)]
    pick  <- which(as.numeric(resos) == resolution)
    if (length(pick) != 1) {
      stop("Couldn't find appropriate resolution.")
    }
    
    bins    <- paste0("/resolutions/", resos[pick], "/", bins)
    pixels  <- paste0("/resolutions/", resos[pick], "/", pixels)
    indexes <- paste0("/resolutions/", resos[pick], "/", indexes)
    
  }
  
  idx <- as.data.table(rhdf5::h5read(file, bins))
  idx[, bin := seq_len(.N)]
  
  if ("KR" %in% colnames(idx) && balancing) {
    idx[, weight := 1 / ifelse(is.finite(KR), KR, 1)]
  } else {
    if (balancing) {
      warning("No balancing vector found.")
    }
    idx[, weight := 1]
  }
  
  idx <- idx[, c("chrom", "start", "end", "bin", "weight"), with = FALSE]
  offset <- rhdf5::h5read(file, indexes)$bin1_offset + 1
  
  start <- bed2idx(idx, location, mode = "start")
  end   <- bed2idx(idx, location, mode = "end")
  # valid_bins <- unlist(Map(`:`, start, end))
  
  o_start <- offset[start]
  o_end   <- offset[end]
  run <- o_end - o_start
  
  # Read locations
  sig <- lapply(seq_along(start), function(i) {
    as.data.table(rhdf5::h5read(
      file, pixels, start = o_start[i], count = run[i]
    ))
  })
  rhdf5::h5closeAll()
  sig <- rbindlist(sig)
  
  # Cool file is 0-indexed
  sig[, bin1_id := bin1_id + 1L]
  sig[, bin2_id := bin2_id + 1L]
  sig <- sig[in_ranges(bin1_id, start, end) & in_ranges(bin2_id, start, end)]
  
  # Do balancing
  if (balancing) {
    weight <- idx$weight
    sig[, count := weight[bin1_id] * weight[bin2_id] * count]
  }
  
  # Do scaling
  if (!is.null(scale_bp)) {
    sig[, count := count * scale_bp / sum(count)]
  }
  
  setnames(sig, 1:3, paste0("V", 1:3))
  set(idx, j = "weight", value = NULL)
  keep <- in_ranges(idx$bin, start, end)
  idx <- idx[keep]
  setnames(idx, 1:4, paste0("V", 1:4))
  set(idx, j = "V1", value = as.character(idx$V1))
  setkeyv(sig, c("V1", "V2"))
  setkeyv(idx, c("V4", "V1", "V2"))
  list(SIG = sig, ABS = idx)
}

in_ranges <- function(value, start, end) {
  Reduce(
    `|`,
    lapply(seq_along(start), function(i) {
      between(value, start[i], end[i])
    })
  )
}
