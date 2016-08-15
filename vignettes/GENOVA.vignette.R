## ------------------------------------------------------------------------
library("GENOVA")

## ----loading, eval = F, message=FALSE------------------------------------
#  WT <- construct.experiment(ICEDpath = '~/Desktop/myExperiment/WT_10000_iced.matrix',
#                            BEDpath = '~/Desktop/myExperiment/WT_10000_abs.bed',
#                            SAMPLENAME = "WT",
#                            COLOR = 'blue')
#  
#  DKO <- construct.experiment(ICEDpath = '~/Desktop/myExperiment/DKO_10000_iced.matrix',
#                             BEDpath = '~/Desktop/myExperiment/WT_10000_abs.bed',
#                             SAMPLENAME = "DKO",
#                             COLOR = 'red')
#  
#  # Parse BEDPE of HICCUPS-loops
#  hiccupsBEDPE <- read.bedpe('~/MEGA/Work/PhD/projects/HiCvis/DKO.hiccups')
#  
#  # Parse BEDPE of HiCseg-TADs
#  hicsegBEDPE <- read.bedpe('~/MEGA/Work/PhD/projects/HiCvis/DKO.HiCseg')

## ----APA, eval = F,message=FALSE, fig.align='center'---------------------
#  # Run APA
#  apaResult.WT <- APA(experiment = WT, loop.bed = hiccupsBEDPE)
#  apaResult.DKO <- APA(experiment = DKO, loop.bed = hiccupsBEDPE)
#  
#  # Create list of APA-results
#  apaResultList <- list(WT = apaResult.WT, DKO = apaResult.DKO)
#  
#  # Visualise APA with ggplot
#  visualise.APA.ggplot(apaResultList, title = 'Comparing WT and DKO: APA')

## ----4C, eval = F,message=FALSE, fig.align='center'----------------------
#  # Run virt4C
#  virt4Cdat <- virtual.4C(experiment = WT, loop.bed = hiccupsBEDPE)
#  
#  # Visualise with base R
#  visualise.virtual.4C(virt4Cdat)

## ----stackedTAD, eval = F, message=FALSE, fig.align='center'-------------
#  # Run stackedTAD
#  stackedTAD.dat.WT <- stackedTAD(experiment = WT,tad.bed = hicsegBEDPE, verbose = F)
#  stackedTAD.dat.DKO <- stackedTAD(experiment = DKO,tad.bed = hicsegBEDPE, verbose = F)
#  
#  # Create list of stackedTAD-results
#  stackedTADResultList <- list(WT = stackedTAD.dat.WT, DKO = stackedTAD.dat.DKO)
#  
#  # Visualise with ggplot
#  visualise.stacked.TAD.ggplot(stackedTADResultList, focus = 1, title = 'Comparing WT and DKO: stacked TAD')

## ----RCP, eval =F ,message=FALSE, fig.align='center'---------------------
#  # Run RCP
#  dat.RCP <- RCP(list(WT, DKO), chromsToUse = c('chr19','chr18'), verbose = F)
#  
#  # Visualise with ggplot
#  visualise.RCP.ggplot(dat.RCP)

