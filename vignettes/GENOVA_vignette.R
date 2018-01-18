## ----global_options, include=FALSE-----------------------------------------
knitr::opts_chunk$set(fig.pos = '!h')

## ----echo = F--------------------------------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(4, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(1.5,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

## ----echo=FALSE, out.width='100%', fig.align='center'----------------------
knitr::include_graphics('/DATA/users/r.vd.weide/github/GENOVA/t1logo')

## ---- echo=T, warning=FALSE, error=F, results='hide'-----------------------
library(GENOVA)

## ----peakEXP3, collapse=F, results='markup', echo = F----------------------
str(Hap1_WT_40kb, width = 60,   vec.len=1, strict.width = 'wrap')

## ----J2G, eval=FALSE, highlight=FALSE--------------------------------------
#  # Convert data from Sanborn et al. normalised at 10kb resoltion:
#  juicerToGenova.py -C ucsc.hg19_onlyRealChromosomes.noChr.chromSizes \
#  -JT ~/bin/juicer/AWS/scripts/juicebox_tools.7.0.jar \
#  -H ~/Downloads/Sanborn_Hap1_combined_30.hic \
#  -R 10000 \
#  -force TRUE \
#  -norm KR \
#  -O Sanborn_Hap1_combined_30.hic_10kb_KR

## ---- echo=F---------------------------------------------------------------
options(scipen = 1)

## ----RCPPLOT1, message=FALSE, fig.wide= T , fig.cap= "RCP. Every facet shows the RCP of one chromosome."----
# Plot RCP: combined
visualise.RCP.ggplot(RCPdata = RCP_out135, 
                     smooth = T, # use a LOESS smoothing
                     combine = F) # Don't merge data from all chromosomes

## ----RCPPLOT2, message=FALSE,   fig.small= T , fig.cap= "RCP. All data combined in one plot."----
# Plot RCP: per-chromosome
visualise.RCP.ggplot(RCPdata = RCP_out135, 
                     smooth = F, # do not use a LOESS smoothing
                     combine = T) # Merge data from all chromosomes

## ----RCPBED, message=FALSE,   fig.small= T , fig.cap= "RCP with BEDs. We can also add BEDs as sites to compute the RCP."----
CTCF = read.delim('data/CTCF_WT_motifs.bed', h = F)
SMC1 = read.delim('data/SMC1_WT_peaks.narrowPeak', h = F)

RCP_out = RCP(experimentList = list(Hap1_WT_40kb), 
               bedList =  list("CTCF" = CTCF, 
                               'Cohesin' =SMC1), 
               chromsToUse = c('chr1'))


visualise.RCP.ggplot(RCP_out)

## ----RCPBED2, message=FALSE,   fig.small= T , fig.cap= "RCP with BEDs. Multiple samples and BEDs can be combined."----
CTCF = read.delim('data/CTCF_WT_motifs.bed', h = F)
SMC1 = read.delim('data/SMC1_WT_peaks.narrowPeak', h = F)

RCP_out = RCP(experimentList = list(Hap1_WT_40kb, Hap1_WAPL_40kb ), 
               bedList =  list("CTCF" = CTCF, 
                               'Cohesin' =SMC1), 
               chromsToUse = c('chr1'))


visualise.RCP.ggplot(RCP_out)

## ----CTCF------------------------------------------------------------------
CTCF = read.delim('data/CTCF_WT_motifs.bed', h = F)
SMC1 = read.delim('data/SMC1_WT_peaks.narrowPeak', h = F)

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(CTCF, 3), caption = 'A data.frame holding a standard BED6 format.'
)

## ----bigwrig, eval=F, echo = F---------------------------------------------
#  library(devtools)
#  install_github(repo ='bigwrig', username =  'jayhesselberth')

## ----biomart---------------------------------------------------------------
## Gene stable ID & Transcript stable ID & Chromosome/scaffold name &
## Transcript start (bp) & Transcript end (bp) & Exon region start (bp) &
## Exon region end (bp) & Strand
martExport = read.delim('data/mart_export.txt.gz', stringsAsFactors = F)
colnames(martExport) = c('ENSG','ENST','chrom' , # change column names
                         'txStart' , 'txEnd' , 
                         'exonStart' , 'exonEnd' , 'strand')
martExport$chrom = gsub(martExport$chrom, # add chr-prefix
                        pattern = '^',
                        replacement = 'chr') 
martExport$strand = ifelse(martExport$strand == 1, '+',"-") # 1/-1 to +/-

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(martExport[,-c(1,2)], 5), caption = 'A data.frame holding the needed columns for plotting genes.'
)

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(WT_TADs[,1:3], 5), caption = 'A data.frame holding a standard BED3 format.'
)

## ----ATA, cache=T----------------------------------------------------------
ATA_Hap1_WT   <- ATA(experiment = Hap1_WT_10kb, verbose = F,
                    tad.bed = WT_TADs) 

ATA_Hap1_WAPL <-  ATA(experiment = Hap1_WAPL_10kb,verbose = F,
                    tad.bed = WT_TADs)

## ---- echo=F---------------------------------------------------------------
options(scipen = 1e9)

## ---- echo=F---------------------------------------------------------------
options(scipen = 1)

## ----APArun, cache=T-------------------------------------------------------
APA_Hap1_WT_extended   <- APA(experiment = Hap1_WT_10kb,
                              loop.bed = WT_Loops_extended)

APA_Hap1_WAPL_extended <- APA(experiment = Hap1_WAPL_10kb,
                              loop.bed = WT_Loops_extended)

## ----SE1-------------------------------------------------------------------
superEnhancers = read.delim('data/homerSuperEnhancers.txt',
                            h = F, 
                            comment.char = "#")

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(superEnhancers[,1:6], 5), caption = "A data.frame holding the output of homer's findPeaks -style super."
)

## ----sesh, echo = F--------------------------------------------------------
sessionInfo()

## ----REF, echo = F, results = 'hide'---------------------------------------
sessionInfo()

