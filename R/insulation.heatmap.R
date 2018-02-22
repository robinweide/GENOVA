# taken from https://github.com/robinweide/RHWlib

#' Plot insulation-heatmaps.
#'
#' Takes a list of insulation-scores and plots a sorted heatmap and/or
#' average profile.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param insulationList A named list of results from genome.wide.insulation().
#' @param bed A BED-df with TAD-calls. Only the first three columns
#' (chrom, start & stop) are important.
#' @param focus Sort on which sample? (give the index-number of the list)
#' @param sortWidth Percentage of columns to sort on (aligned on center)
#' @param whatToPlot Do you want to plot a profile, heatmap or both?
#' @param profileFunct Function to make profile-plots (mean, median, sum)
#' @param flank Amount of flanking bp.
#' @param zlim The heatmap-zlims c(min,max).
#' @param profileCols Vector of line-colors for profiles
#' @param profileZlim The profile-zlims c(min,max).
#' @param heatmapCols A vector of colors for the heatmap-gradient.
#' @param rmNonVariateRows Ignore rows with an sd of 0.
#' @param verbose Should this function be chatty?
#' @return A plot plus an (invisible) dataframe of the underlying matrix.
#' @examples
#' # Get the insulation score with window-size 25 of two experiments.
#' WT_10kb_ins = genome.wide.insulation(hic = WT_10kb, window.size = 25)
#' SCC4_10kb_ins = genome.wide.insulation(hic = SCC4_10kb, window.size = 25)
#'
#' # Store insulation-scores in a named list
#' inList = list(WT = WT_10kb_ins, SCC4 = SCC4_10kb_ins )
#'
#' # Plot the heatmap and profile and store the underlying df in df.out.
#' df.out = insulation.heatmap(inList, bed = WT_TADs,
#'                             zlim = c(-.5,0.25), profileZlim = c(-.75,-.1))
#' @export
insulation.heatmap <- function(insulationList, bed,  focus = 1, sortWidth = 10,
                               whatToPlot = 'both',  profileFunct = mean,
                               flank=500e3, zlim = c(-1,1), profileCols = NULL,
                               profileZlim = NULL, heatmapCols = c("#f03b20",
                                                                   "#ffeda0",
                                                                   "white",
                                                                   "#31a354"),
                               rmNonVariateRows = T, verbose = F ){

  if(!whatToPlot %in% c('profile', 'heatmap', 'both')){
    stop('whatToPlot should be "profile", "heatmap" or "both".')
    }

  require(ggplot2)
  require(reshape2)
  # Get matrixes
  sampleList <- list()
  NVrows = c()
  for(i in 1:length(insulationList)){
    insdDat <- insulationList[[i]]
    insName <-  names(insulationList)[i]
    results <- align.insulation( insdDat, bed, flank.length=flank,
                                 verbose = verbose)
    #! place zlims and remove infs !
    # results[results < zlim[1]] = zlim[1]
    # results[results > zlim[2]] = zlim[2]
    NVrows = c(NVrows, which((apply(results, 1, sd) == 0) ))
    sampleList[[i]] <- results
    names(sampleList)[i] <- insName
  }
  sampleNames <- names(insulationList)

  perSample = NULL

  # sorting on focus
  foundOrder = NULL

  # remove nonVariateRows
  if(rmNonVariateRows){
    TNV = table(NVrows)
    nonVariateRows = names(TNV)[TNV == length(insulationList)]
    nonVariateRows = as.numeric(nonVariateRows)

    for(i in 1:length(sampleList)){
      sampleList[[i]] = sampleList[[i]][-nonVariateRows,]
    }
  }

  if(focus > length(insulationList)){
    stop("Focus is outside of the number of samples!")
  } else {
    focusYoungGrasshopper <- sampleList[[focus]]

    # which of columns to sort on?
    perSample = ncol(focusYoungGrasshopper)
    nSortCol <-  ((perSample)/100)*sortWidth
    halfnSortCol <- round(nSortCol/2)
    midpoint <- ((ncol(focusYoungGrasshopper)-1)/2)+1
    sortCols <- (midpoint-halfnSortCol):(midpoint+halfnSortCol)

    rowSums <- apply(as.matrix(focusYoungGrasshopper)[, sortCols],1, mean)
    foundOrder <- base::order(rowSums,decreasing = T) # in/de flipped for ggplot
  }

  # order individual samples
  for(i in 1:length(sampleList)){
    sampleList[[i]] <- sampleList[[i]][foundOrder,]
  }

  if(verbose){message("Melting...")}
  dfList <- NULL
  for(i in 1:length(sampleList)){
    melted <- reshape2::melt(t(sampleList[[i]]))
    melted$Var1 <- as.numeric(as.factor(melted$Var1))
    melted$Var2 <- as.numeric(as.factor(melted$Var2))
    melted$sample <- names(sampleList)[i]
    dfList[[i]] <- melted
  }
  df <- data.table::rbindlist(dfList)
  df$sample <- factor(df$sample, levels = sampleNames)

  # plot profiles
  groupedDF <- dplyr::group_by(df, sample, Var1)
  groupedDF[is.na(groupedDF$value), "value"] = 0
  minVal = min(groupedDF[groupedDF$value != -Inf, "value"])
  groupedDF[groupedDF$value == -Inf, "value"] = minVal
  profDF <- dplyr::summarize(groupedDF, val = profileFunct(value))

  if(is.null(profileCols)){
    profileCols <- rep('black', length(insulationList))
  }



  # finding ticks
  flankingSize = (perSample -1 )/2
  upTickPos = flankingSize/2
  downTickPos = perSample - upTickPos
  upTickLab = -1*((flank/2)/1e3)
  downTickLab = (flank/2)/1e3

  tickLabs = c(paste0(upTickLab, 'kb'),  "", paste0(downTickLab, 'kb'))
  tickPoss = c(upTickPos+1,  perSample - flankingSize ,  downTickPos )+1

  PRFLS <- NULL
  if(is.null(profileZlim)){

    #profileZlim = range(profDF$val) + c(-0.05, 0.05) * diff(range(profDF$val))
    profileZlim = zlim
    PRFLS <- ggplot2::ggplot(profDF, ggplot2::aes(Var1, val, col = sample)) +
      ggplot2::geom_line()  +
      ggplot2::scale_x_continuous(expand=c(0,0),
                                  breaks = tickPoss,
                                  labels = tickLabs) +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = profileZlim ) +
      ggplot2::labs( x = '', y = '') + ggplot2::guides(col = F) +
      ggplot2::scale_colour_manual(values = profileCols)+
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::facet_grid(. ~ sample ) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank())
  } else {
    PRFLS <- ggplot2::ggplot(profDF, ggplot2::aes(Var1, val, col = sample)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(expand=c(0,0), breaks = tickPoss,
                                  labels = tickLabs) +
      ggplot2::scale_y_continuous(expand=c(0,0), limits = profileZlim ) +
      ggplot2::labs( x = '', y = '') + guides(col = F) +
      ggplot2::scale_colour_manual(values = profileCols) +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::facet_grid(. ~ sample ) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank())
  }

  # plot heatmaps

  df$value[df$value < zlim[1]] = zlim[1]
  df$value[df$value > zlim[2]] = zlim[2]

  HTMPS <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = value)) +
    ggplot2::geom_raster(interpolate = T) +
    ggplot2::facet_grid(. ~ sample) +
    ggplot2::scale_x_continuous(expand=c(0,0), breaks =tickPoss,
                                labels = tickLabs) +
    ggplot2::scale_y_continuous(expand=c(0,0), breaks = NULL) +
    ggplot2::labs( x = '', y = '') +
    ggplot2::guides(fill = F) + ggplot2::theme_linedraw()+
    ggplot2::scale_fill_gradientn(colours = heatmapCols)+
    ggplot2::theme(strip.text = ggplot2::element_blank() ,
          strip.background = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin = ggplot2::unit( c(0,0,0,0) , units = "lines" ))



  if(whatToPlot == 'profile'){
    plot(PRFLS)
    invisible(df)
  } else if(whatToPlot == 'heatmap'){
    plot(HTMPS)
    invisible(df)
  } else {
    PRFLS <- NULL
    if(is.null(profileZlim)){
      profileZlim = zlim
      PRFLS <- ggplot2::ggplot(profDF, ggplot2::aes(Var1, val, col = sample)) +
        ggplot2::geom_line()  +
        ggplot2::scale_x_continuous(expand=c(0,0), breaks = tickPoss,
                                    labels = c("","","")) +
        ggplot2::scale_y_continuous(expand=c(0,0), limits = profileZlim ) +
        ggplot2::labs( x = '', y = '') + ggplot2::guides(col = F) +
        ggplot2::scale_colour_manual(values = profileCols)+
        ggplot2::theme(aspect.ratio = 1) +
        ggplot2::facet_grid(. ~ sample ) + ggplot2::theme_linedraw()+
        ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))+
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank())
    } else {
      PRFLS <- ggplot2::ggplot(profDF, ggplot2::aes(Var1, val, col = sample)) +
        ggplot2::geom_line() +
        ggplot2::scale_x_continuous(expand=c(0,0), breaks = tickPoss,
                                    labels = c("","","")) +
        ggplot2::scale_y_continuous(expand=c(0,0), limits = profileZlim ) +
        ggplot2::labs( x = '', y = '') + guides(col = F) +
        ggplot2::theme_linedraw()+
        ggplot2::scale_colour_manual(values = profileCols) +
        ggplot2::theme(aspect.ratio = 1) +
        ggplot2::facet_grid(. ~ sample ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank())
    }


    grid::grid.newpage()
    grid::grid.draw(rbind(ggplot2::ggplotGrob(PRFLS),
                          ggplot2::ggplotGrob(HTMPS),
                          size = "last"))
    invisible(df)
  }

}




align.insulation.chrom <- function( ins.data, bed, flank = 10 ){
  n <- findInterval(bed[,2]-1, ins.data[,2])
  bool <- n > flank & n + flank < nrow(ins.data)
  boolF <- which(bool == F)
  boolD <- diff(boolF)

  rows2addFront = which(boolD > 1)
  rows2addBack = length(boolF)-rows2addFront

  n <- n[n > flank & n + flank < nrow(ins.data)]

  #create a vector for the flanking sequence
  add <- -flank:flank
  #repeat it for every element in n
  add.long <- rep(add, length(n))
  #repeat every element in n for the length of add
  n.long <- rep(n, each=length(add))

  #add the selected value and the flanking sequences
  i.sel <- n.long + add.long

  #put the selected elements in matrix
  align.mat <- matrix(ins.data[i.sel,4], nrow=length(n), ncol=length(add),
                      byrow=T)

  # add zeroes for first boundaries too close
  if(length(rows2addFront) > 0){
    align.mat = rbind(matrix(rep(0, ncol(align.mat) * rows2addFront),
                             ncol =  ncol(align.mat)),
                      align.mat)
  }

  if(length(rows2addBack) > 0){
    align.mat = rbind(align.mat,
                      matrix(rep(0, ncol(align.mat) * rows2addBack),
                             ncol =  ncol(align.mat)))
  }

  return(align.mat)
}

align.insulation <- function(ins.data, bed, flank.length = 200e3, verbose = F){
  #elaborate way of selecting the resolution of the Hi-C matrix
  res <- as.numeric(names(tail(sort(table(ins.data[,3]-ins.data[,2])),1)))
  flank = flank.length/res

  chrom.vec <- unique(ins.data[,1])
  align.mat <- NULL
  for( chr in chrom.vec){
    if(verbose){message("Analyzing ", chr, "\\r")}

    sub.mat <- align.insulation.chrom( ins.data[ins.data[,1]==chr,],
                                       bed[bed[,1]==chr,], flank = flank )
    align.mat <- rbind(align.mat, sub.mat)
  }
  align.mat
}
