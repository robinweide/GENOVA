---
title: "GENOVA: explore the Hi-C's"
author:
- name: Robin H. van der Weide
  affiliation:
  - Division of Gene Regulation, the Netherlands Cancer Institute
- name: Elzo de Wit
  affiliation: Division of Gene Regulation, the Netherlands Cancer Institute
  email: e.d.wit@nki.nl
date: "2018-01-03"
output:
  BiocStyle::pdf_document
vignette: >
  %\VignetteIndexEntry{GENOVA: explore the Hi-C's}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}

bibliography: GENOVA.bib
---
# Loading data


```r
library(GENOVA)
```


## Index-based
GENOVA expects two input files: the signal- and the index-file. Signal-files have three columns (bin1, bin2, contactCount) and index-files have four (chromosome, start, end, bin). These are the default output of the Hi-C mapping pipeline HiC-Pro [@Servant2015], where they are called \*.matrix and \*.bed. The files are expected to be genome-wide and may be corrected with ICE-normalisation.

## Juicebox
We added a convenience script **juicerToGENOVA.py**, to load files from Juicerbox (.hic files). This allows for a fast conversion to signal- and index-files from, for example, data from Sanborn et al.[-@Sanborn2015]:


```text
# Convert data from Sanborn et al. normalised at 10kb resoltion:
juicerToGenova.py -C ucsc.hg19_onlyRealChromosomes.noChr.chromSizes \
-JT ~/bin/juicer/AWS/scripts/juicebox_tools.7.0.jar \
-H ~/Downloads/Sanborn_Hap1_combined_30.hic \
-R 10000 \
-force TRUE \
-norm KR \
-O Sanborn_Hap1_combined_30.hic_10kb_KR 
```

## Recommended resolutions
To ensure computational strain and time is kept to a minimum, we recommend different resolutions for different functions. More experienced users are free to deviate, while keeping in mind that these datasets are quite memory-heavy. 

Function   | Resolution
---------- | ----------
APA | 10kb-20kb
ATA  | 10kb-40kb
cisTotal.perChrom | 500kb-1Mb
chromosomeMatrix | 500kb-1Mb
RCP | 40kb-500kb
intra.inter.TAD.contacts | 20kb - 40kb 
PE-SCAn | 20kb-40kb
hic.matrixplox | $\frac{width\ in\ bp\ of\ window}{500}$ 
: (\#tab:table) Recommended resolutions. These will provide optimal resource/result tradeoffs.

## construct.experiment
Every Hi-C experiment will be stored in an experiment-object. This is done by invoking the `construct.experiment` function. Inside, several sanity checks will be performed and the data is normalised to the total sum of reads. You can also add centromere-information in the form of a three-column data.frame:

```r
# Please make sure that the chromosome-names match.
centromeres = read.delim('data/hg19_cytobandAcen.bed', 
                         sep = '\t', 
                         h = F, 
                         stringsAsFactors = F)

head(centromeres)
##      V1        V2        V3
## 1  chr1 121500000 128900000
## 2 chr10  38000000  42300000
## 3 chr11  51600000  55700000
## 4 chr12  33300000  38200000
## 5 chr13  16300000  19500000
## 6 chr14  16100000  19100000
```

For this example, we are going to use the Hi-C maps of WT and $\Delta$WAPL Hap1 cells from Haarhuis et al. [-@Haarhuis2017]. Since the genome-wide analyses do not need very high-resolution data, we will construct both 10kb and 100kb resolution experiment-objects. 

```r
Hap1_WT_10kb <- construct.experiment(
                                signalPath = 'data/WT_10000_iced.matrix', 
                                indicesPath = 'data/WT_10000_abs.bed', 
                                name = "WT",
                                centromeres = centromeres,
                                color = "black", 
                                comments = "This is an optional memo-field.")

Hap1_WAPL_10kb <- construct.experiment(
                                signalPath = 'data/WAPL_10000_iced.matrix', 
                                indicesPath = 'data/WAPL_10000_abs.bed', 
                                name = "WAPL",
                                centromeres = centromeres,
                                color = "red")

Hap1_WT_100kb <- construct.experiment(
                                signalPath = 'data/WT_100000_iced.matrix', 
                                indicesPath = 'data/WT_100000_abs.bed', 
                                name = "WT",
                                centromeres = centromeres,
                                color = "black", 
                                comments = "This is an optional memo-field.")

Hap1_WAPL_100kb <- construct.experiment(
                                signalPath = 'data/WAPL_100000_iced.matrix', 
                                indicesPath = 'data/WAPL_100000_abs.bed', 
                                name = "WAPL",
                                centromeres = centromeres,
                                color = "red")

Hap1_WT_1MB <- construct.experiment(
                                signalPath = 'data/WT_1000000_iced.matrix', 
                                indicesPath = 'data/WT_1000000_abs.bed', 
                                name = "WT",
                                centromeres = centromeres,
                                color = "black", 
                                comments = "This is an optional memo-field.")

Hap1_WAPL_1MB <- construct.experiment(
                                signalPath = 'data/WAPL_1000000_iced.matrix', 
                                indicesPath = 'data/WAPL_1000000_abs.bed', 
                                name = "WAPL",
                                centromeres = centromeres,
                                color = "red")
```

In this object are several slots:

```
## List of 9
## $ ICE :Classes 'data.table' and 'data.frame': 67934954 obs.
##    of 3 variables:
## ..$ V1: int [1:67934954] 1 1 ...
## ..$ V2: int [1:67934954] 1 7 ...
## ..$ V3: num [1:67934954] 1556 ...
## ..- attr(*, ".internal.selfref")=<externalptr>
## ..- attr(*, "sorted")= chr [1:2] "V1" ...
## $ ABS :'data.frame': 30971 obs. of 4 variables:
## ..$ V1: chr [1:30971] "chrM" ...
## ..$ V2: int [1:30971] 0 0 ...
## ..$ V3: int [1:30971] 16571 100000 ...
## ..$ V4: int [1:30971] 1 2 ...
## $ NAME : chr "WT"
## $ RES : num 1e+05
## $ CHRS : chr [1:25] "chrM" ...
## $ COL : chr "black"
## $ COMM : chr "This is an optional memo-field."
## $ MASK : logi(0)
## $ CENTROMERES:'data.frame': 24 obs. of 3 variables:
## ..$ V1: chr [1:24] "chr1" ...
## ..$ V2: int [1:24] 121500000 38000000 ...
## ..$ V3: int [1:24] 128900000 42300000 ...
```

# Genome-wide analyses
## Cis-quantification
Another important quality-metric is the fraction of *cis*-interactions. Work by the group of Amons Tanay showed that the expected amount of intra-chromosomal contacts is 90-93% [@Olivares-Chauvet2016]. Assuming that any extra inter-chromosomal contacts are due to dubris/noise, the user might aspire to get the *cis*-percentages as close to 90% as possible.

```r
cisChrom_out <- cisTotal.perChrom(Hap1_WT_1MB)
```

\begin{smallfigure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/cis-1} \caption{Fraction of cis-contacts per chromosome.}(\#fig:cis)
\end{smallfigure}

Of course, one can also inspect the results per chromosome more closely:

```r
plot(cisChrom_out$perChrom, las=2)
abline(h = cisChrom_out$genomeWide)
```

\begin{figure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/cis2-1} \caption{Fraction of cis-contacts per chromosome. Chromosomes 9, 15, 19 \& 22 have translocations, which leads to the appearance of more trans-interactions than similar-sized chromosomes.}(\#fig:cis2)
\end{figure}

## chromosome plots
To find possible translocations and/or flawed mapping, we can plot the genome-wide enrichment of interactions between all combinations of chromosomes. The values in the matrix are $log10(observed/expected)$.

```r
par(pty ='s')
# Lets remove mitochondrial and Y-chromosomal contacts
chromosomeMatrix(Hap1_WT_1MB, remove = c("chrM","chrY"))
```

\begin{figure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/chromMat-1} \caption{Chromosome matrix. The two known translocations of Hap1 cells are easily seen (15-19 \& 9-22).}(\#fig:chromMat)
\end{figure}

## RCP
The Relative Contact Probability computes the contact probability as a function of genomic distance, as described in [@Lieberman-Aiden2009]. This can be computed for a specific set of chromosomes or genome-wide. To be able to ignore centromeric contacts (which have a abberant RCP), centromeric information is need. This is taken from the experiment-object or found emperically by comparing trans-interactions.


```r
RCP_out <- RCP(experimentList = list(Hap1_WT_100kb, Hap1_WAPL_100kb), 
               chromsToUse = c('chr1','chr2','chr3'))
```

The user can decide to plot the RCP per chromosome. If the data is sparse, a LOESS-smooting could be convenient. It takes the color and name from the experiment-objects. 

```r
# Plot RCP: combined
visualise.RCP.ggplot(RCPdata = RCP_out, 
                     smooth = T, # use a LOESS smoothing
                     combine = F) # Don't merge data from all chromosomes
```

\begin{figure*}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/RCPPLOT1-1} \caption{RCP. Every facet shows the RCP of one chromosome.}(\#fig:RCPPLOT1)
\end{figure*}

It is also possible to combine all available data into a genome-wide RCP-plot.

```r
# Plot RCP: per-chromosome
visualise.RCP.ggplot(RCPdata = RCP_out, 
                     smooth = T, # use a LOESS smoothing
                     combine = T) # Merge data from all chromosomes
```

\begin{smallfigure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/RCPPLOT2-1} \caption{RCP. All data combined in one plot.}(\#fig:RCPPLOT2)
\end{smallfigure}


## matrix plots


# TADs
GENOVA has a large repetoire of functions to analyse TADs. We use the TAD-calls of WT Hap1 20-kb matrices from Haarhuis et al. [-@Haarhuis2017], generated with HiCseg [@Levy-Leduc2014].

```r
# Lets load in some TAD-annotations from HiC-seg
WT_TADs = read.delim('data/WT_hicseg_TADs.bed', h = F)
```


```
##     V1     V2      V3
## 1 chr1  50000  120000
## 2 chr1 120000  210000
## 3 chr1 210000  240000
## 4 chr1 240000  610000
## 5 chr1 610000  900000
## 6 chr1 900000 1270000
```


## ATA

```r
ATA_Hap1_WT   <- ATA(experiment = Hap1_WT_10kb,
                tad.bed = WT_TADs) 
ATA_Hap1_WAPL <-  ATA(experiment = Hap1_WAPL_10kb,
                 tad.bed = WT_TADs)
```

We can use `visualise.ATA.ggplot` to combine the ATA-results.


```r
visualise.ATA.ggplot(stackedlist = list('WT' = ATA_Hap1_WT, 
                                        'WAPL' = ATA_Hap1_WAPL), # a named list
                     title = "Hap1 Hi-C vs WT TADs from HiCseg", 
                     zlim1 = c(0,75),
                     zlim2 = c(-5,5), 
                     focus = 1) # which entry to use as comparison
```

\begin{figure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/ATAplot-1} \caption{ATA. In the WAPL-knockout, we see a decrease of contacts within the TAD, but an increase at the corner.}(\#fig:ATAplot)
\end{figure}

## TAD+N

```r
TAD_N_WT   <- intra.inter.TAD.contacts(TAD = WT_TADs, 
                                       max.neighbor = 10, 
                                       exp = Hap1_WT_10kb)
TAD_N_WAPL <- intra.inter.TAD.contacts(TAD = WT_TADs, 
                                       max.neighbor = 10, 
                                       exp = Hap1_WAPL_10kb)
```

We can compute the enrichment of contacts between TADs with the `differential.TAD.dotplot`-function. 

```r
differential.TAD.dotplot(exp1 = TAD_N_WT, # denominator
                         exp2 = TAD_N_WAPL) # numerator
```

\begin{figure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/iiDIFF-1} \caption{Differential TAD-analysis. Experiment 2 (WAPL) has more interactions between neighbouring TADs}(\#fig:iiDIFF)
\end{figure}

# Loops
For this section, we are using the extended loops from Haarhuis et al. [-@Haarhuis2017]. These are the anchor-combinations of the merged loop-calls of WT Hap1 5-, 10- and 20-kb matrices, generated with HICCUPS [@Rao2014].

```r
WT_Loops = read.delim('data/WT_3Mb_extended_loops.bed', h = F)
```


```
##      V1     V2     V3    V4      V5      V6
## 1 chr11 875000 900000 chr11 2020000 2025000
## 2 chr11 875000 900000 chr11 2162500 2187500
## 3 chr11 875000 900000 chr11 2020000 2025000
## 4 chr11 875000 900000 chr11 2020000 2030000
## 5 chr11 875000 900000 chr11 1940000 1945000
## 6 chr11 875000 900000 chr11 1965000 1970000
```


## APA
Explain smalltreshold

```r
APA_Hap1_WT   <- APA(experiment = Hap1_WT_10kb,
                     smallTreshold = 300e3,
                     loop.bed = WT_Loops)
APA_Hap1_WAPL <- APA(experiment = Hap1_WAPL_10kb,
                     smallTreshold = 300e3,
                     loop.bed = WT_Loops)
```


We can use `visualise.APA.ggplot` to combine the APA-results.

```r
visualise.APA.ggplot(APAlist = list('WT' = APA_Hap1_WT, 
                                    'WAPL' = APA_Hap1_WAPL), # a named list
                     title = "Hap1 Hi-C vs extended loops", 
                     zTop = c(1,12), # set the zlims of the upper row
                     zBottom = c(-8.33,8.33),
                     focus = 1) # which item in APAlist to use as comparison
```

\begin{figure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/APAplot-1} \caption{APA. In the WAPL-knockout, we see an increase of contacts at the loop.}(\#fig:APAplot)
\end{figure}


## PE-SCAn

```r
superEnhancers = read.delim('data/homerSuperEnhancers.txt', h = F, comment.char = "#")

WT_PE_OUT = PESCAn(exp = Hap1_WT_10kb, bed = superEnhancers[superEnhancers$V2=='chr2',2:4])
WAPL_PE_OUT = PESCAn(exp = Hap1_WAPL_10kb, bed = superEnhancers[superEnhancers$V2=='chr2',2:4])
visualise.PESCAn.ggplot(list(WT = WT_PE_OUT, WaplKO = WAPL_PE_OUT), resolution = 10e3, smooth = F)
```

\begin{smallfigure}
\includegraphics{/tmp/Rtmpvs8F7z/preview-88f7488130c3.dir/GENOVA_vignette_files/figure-latex/PESCAn-1} \caption{PE-SCAn. There is a pairwise enrichment of contacts between Superenhancers, compared to shifted regions in both the WT and Wapl samples.}(\#fig:PESCAn)
\end{smallfigure}

# Session info

```
## R version 3.4.3 (2017-11-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.3 LTS
## 
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_2.2.1     viridis_0.4.0     viridisLite_0.2.0 bindrcpp_0.2     
## [5] GENOVA_0.9.5      BiocStyle_2.6.1  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.14           pillar_1.0.1           compiler_3.4.3        
##  [4] GenomeInfoDb_1.14.0    plyr_1.8.4             XVector_0.18.0        
##  [7] bindr_0.1              bitops_1.0-6           tools_3.4.3           
## [10] zlibbioc_1.24.0        digest_0.6.12          evaluate_0.10.1       
## [13] tibble_1.4.1           gtable_0.2.0           pkgconfig_2.0.1       
## [16] rlang_0.1.6            yaml_2.1.16            parallel_3.4.3        
## [19] gridExtra_2.3          GenomeInfoDbData_1.0.0 dplyr_0.7.4           
## [22] stringr_1.2.0          knitr_1.18             S4Vectors_0.16.0      
## [25] IRanges_2.12.0         stats4_3.4.3           rprojroot_1.3-1       
## [28] grid_3.4.3             glue_1.2.0             data.table_1.10.4-3   
## [31] R6_2.2.2               rmarkdown_1.8.5        bookdown_0.5          
## [34] reshape2_1.4.2         magrittr_1.5           codetools_0.2-15      
## [37] backports_1.1.2        scales_0.4.1           htmltools_0.3.6       
## [40] BiocGenerics_0.24.0    GenomicRanges_1.30.1   assertthat_0.2.0      
## [43] colorspace_1.3-2       labeling_0.3           stringi_1.1.5         
## [46] RCurl_1.95-4.9         lazyeval_0.2.0         munsell_0.4.3
```

# References

