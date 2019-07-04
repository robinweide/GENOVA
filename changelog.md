# Changelog
All notable changes to this project will be documented in this file.
To generate this file until 0.9.98, I have used [Github](https://github.com/robinweide/GENOVA/compare/8edef82..9951098). 

## [0.9.995] - 20-3-2019

### Added
- HiCseg.callTAD: a function to call TADs with HiCseq (Levy-Leduc) per chromosome-arm or windowed chromosomes.
- construct.experiment: added an option to get Z-score normalised experiment-objects.
- hic.matrixplot: added an option to use Z-score normalised experiment-objects.
- insulation.heatmap: added option to set leftlost bin to 0
- insulation.heatmap: added option to use a borders-file (also for peaks!)

### Changed
- APA: sped up by using a vapply loop (@MarijneMia)
- ATA: sped up by using a vapply loop
- select.subset: sped up the pos-lookup
- visualise.APA.ggplot: better divergent colors (same as Z-score colors in hic.matrixplot)
- visualise.ATA.ggplot: better divergent colors (same as Z-score colors in hic.matrixplot)
- genome.wide.insulation: added speudocounts
- compartment.score: uses a more robust handing of comparable tracks.
- genome.wide.insulation::matrix.insulation: shifted output half a bin upstream.
- genome.wide.insulation: forcing even window-sizes to get score between bins.
- insulation.heatmap: nicer plotting of profile and heatmap without grid
- insulation.heatmap: refactored the scripts.

## [0.9.98] - 27-02-2019
Github-hash: 9951098

### Added
- ATA: outputs also the used TADs.bed 
- HiC_matrixplot: smoothing of the regions with no data (e.g. the white stripes) has been implemented. This is done by filling these bins with the result of a Nadaraya/Watson normalization of the kernel. Set `smoothNA` to true to use and `smoothBandwidth` to tweak 
- HiC_matrixplot: If `chip.yMax` is NULL, a warning will be given for the used yMaxes.
- PE-SCAn: A maximal distance can be given with `maxDist`.
- PE-SCAn: verbosity can be limited off with `verbose = F`
- PE-SCAn: Catches errors when no matrices are found by cov2d.
- PE-SCAn: A treshold for the minimum amount of BED-entries per chromosome is added: `minComparables`.
- PE-SCAn: Shifted bed-entries bigger than their chromosome are now fixed (they "wrap around" to the start)
- PE-SCAn: Returns not only the O/E matrix, but a list with a O/E score-matrix (if shift != 0), otherwise an observed score-matrix, the underlying signal and background-matrices and the shift used.
- visualise.PESCAn.ggplot: several checks to see if the input-data is comparable (including shift).
- visualise.PESCAn.ggplot: colorscales for O/E and O are better (i.e. no log2-divergent scale for shift==0 matrices).
- cis.compartment.plot: a second experiment-object can be used. This is plotted in the lower-left corner. Note: only the first matrix is invisibly returned.
- cis.compartment.plot: smoothing of the regions with no data (e.g. the white stripes) has been implemented. This is done by filling these bins with the result of a Nadaraya/Watson normalization of the kernel. Set `smoothNA` to true to use and `smoothBandwidth` to tweak 
- getTADstats: a new function to get some stats on ATA-results. Still quite buggy.
- saddle: Forces user for either CS or chip.
- README now contains actual text.

### Changed
- HiC_matrixplot: `bed.col` and `bw.col` are merged to chip.col. This argument can take a vector of four colours, parallel to the `chip`-argument.
- HiC_matrixplot: `yMax` is renamed to `chip.ymax`.
- visualise.PESCAn.ggplot: works woth new output of PE-SCAn
- PE-SCAn: cov2d loops over a vectorised list without duplicates, which speeds up the code.
- PE-SCAn: cleaned up cov2d's return-code
- quantifyAPA: the `enrichmentType` argument lets you choose between pixel/mean(backgroundRegions) and the fraction of loops with more than 50% higher signal than the background.
- visualise.compartmentStrength: at least one bin is used (bugfix). 
- saddle: extra effort to use the CS is taken by matching the CS-vector and the eigen-vector. This will fix problems due to big centromeres or shaky/different centromere-calls.
- vignette: fixed typo's

## [0.9.971] - 13-06-2018
Github-hash: 8edef82

## [0.9.97] - 11-06-2018
Github-hash: 35772ec

## [0.9.96] - 22-05-2018
Github-hash: 7deebca



