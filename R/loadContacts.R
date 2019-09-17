loadContacts = function(signalPath = NULL, indicesPath = NULL, juicerPath = NULL, resolution = 10e3, name){

  doJuicer = F
  doHiCpro = F

  # what are we dealing with?
  if(!is.null(juicerPath)){

    # a juicerpath is given. Does it exist?
    if(!file.exists(juicerPath)){
      stop("juicerPath doesn't point to an existing .hic-file.")
    }

    # is strawr installed?
    try_require('strawr', "loadContacts", source = 'github')

    doJuicer = T

    # did you also load in hicpro?
    if(!is.null(signalPath)){
      warning('Both hicpro- and juicebox-files provided. Will use juicebox.')
    }
  } else if(!is.null(signalPath)){

    # a signalPath is given. Does it exist?
    if(!file.exists(signalPath)){
      stop("signalPath doesn't point to an existing file.")
    }

    if(!file.exists(indicesPath)){
      stop("indicesPath doesn't point to an existing file.")
    }

    doHiCpro = T
  } else {
    stop('please provide either signalPath and indicesPath or a juicerPath.')
  }

  if(doJuicer){
    sig_ind = loadJuicer(juicerPath, resolution)
  }

  if(doHiCpro){
    sig_ind = loadHiCpro(signalPath, indicesPath)
  }

}


try_require('rhdf5', 'loadCooler')


# taken from ggplot
try_require <- function(package, fun, source = NULL) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }

  if(source == 'BIOC'){
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from Bioconductor and try again.", call. = FALSE)
  } else   if(source == 'github'){
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install from github and try again.", call. = FALSE)
  } else {
    stop("Package `", package, "` required for `", fun , "`.\n",
         "Please install and try again.", call. = FALSE)
  }

}

