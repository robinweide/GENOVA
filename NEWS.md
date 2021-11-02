# GENOVA (development version)

- The `tornado_insulation()` function now accepts a list of bed-like data.frames 
  (#266).
- Exposed the `strand` argument in `ARA()`, making it easier to do strand aware
  ARAs (#268)
- Renamed `cut.off` to `colour_lim` and allow `c(minimum, maximum)` 
  specification, rather than `cut.off = maximum` in `hic_matrixplot()` (#260).
- Added a `colour_bar` legend option in `*_matrixplot()` family (#260).
- The `par()` options now reset to original state after calling `*_matrixplot()`
  family of functions.
- Fixed a bug with `NA`s in `saddle()` (#263).
- Removed {bigwrig} as a dependency in favour of Bioconductor's {rtracklayer}.
- Fixed a bug wherein the correct resolution couldn't be found in `.mcool` files
  due to scientific formatting of large numbers (#248, #249).
- Fixed a bug in `insulation_score()` wherein formatting of sample names could
  throw errors (#246).
- In `load_contacts()`, you can now set `centromeres = FALSE` to signal that
  the data doesn't contain centromeres (#241).
- Fixed a balancing bug in .hic file loading (#225).

# GENOVA v1.0 27-01-2021 (The Cartographer)

This version marks a large refactoring of the code base in an effort to increase
consistency. However, it also breaks a lot of old code. It is the version that
was used during submission of the 
[GENOVA publication](https://doi.org/10.1093/nargab/lqab040).

## Data representation

- A loaded Hi-C experiment is now dubbed a <contacts> class object.
- Added support for loading data from the cooler and juicer pipelines. Juicer
  import requires the {strawr} dependency and cooler import requires the {Rhdf5} 
  dependency.
- Centromere information is automatically estimated while loading the data by
  looking for the largest stretch of empty bins.
- Added the `sync_indices()` function to harmonise data from various 
  pipelines.

## Analysis functions

- Results from analysis functions have been dubbed <discovery> class objects.
- For quantification of analysis, the `quantify()` generic was added with 
  methods for several <discovery> objects.
- Likewise, for visualising results, the `visualise()` generic was added with
  methods for several <discovery> objects. Additionaly, some base R equivalents
  are available as <discovery> methods to the `plot()` function.
- New `bundle()` and `unbundle()` functions to more easily combine <discovery>
  objects.
- New `ARA()` functions for aggregating on-diagonal regions, convenient for
  stripe analysis.
- New `CSCAn()` function for crosswise intersections of different sets of 
  regions.
- New `anchors_*()` family of functions for related analysis functions for 
  repeated lookup analyses, to make these analyses more flexible.
- New `anchors_extendedloops()` function to recapitulate the Haarhuis *et al*. 
  (2017) strategy for loop anchors.
- Standardised the `APA()`, `ATA()` and `PESCAn()` functions to use a common
  pattern for repeatedly looking up regions in the Hi-C map.
- Speed-up of the repeated lookup functions by making better use of {data.table}
  joins.
- The `insulation_score()` algorithm was sped up for larger data.
- TAD calling via `call_TAD_insulation()` was sped up and closer to the 
  description of the algorithm in Crane *et al.* (2015).
- Compartment strength is no longer calculated in the deprecated 
  `visualise.compartmentStrength()`, but in the `quantify.saddle_discovery()`
  method.
- Refactored the `RCP()` function.
  
## Plotting

- Added colour palettes for sequential and divergent values 
  (see `?GENOVA_colours`). The default sequential palette can be changed
  using `options("GENOVA.colour.palette" = {a character vector of colours})`.
- Wrappers for these colour palettes in {ggplot2} scales are in the 
  `scale_{colour/fill}_GENOVA()` and `scale_{colour/fill}_GENOVA_div()` 
  functions for sequential and divergent palettes respectively.
- Added the `pyramid()` and `pyramid_difference()` function for plotting a
  Hi-C map at 45 degree angles.
- The `pyramid()` plot can be annotated by the `add_bed_graph()`, 
  `add_bed_track()`, `add_ctcf_sites()` functions.
- The `pyramid()` plot can be annotated by compartment score, insulation score,
  directionality index, domainograms, virtual 4C or arbitrary ggplot2 layer 
  using the `as_track()` function.
- Most `image()`-based plotting now have the `rasterise` argument.
- Replaced the `visualise.{analysis}.ggplot()` pattern functions by a 
  `visualise()` generic with methods for <discovery> objects. 

## Renamed

Generally, functions have been renamed to prevent `dot.case` functions implying
S3 methods or `lowerCamelCase` implying S4 functions. Instead they now all use
the `snake_case`.

| Old name                     | New name                   |
|------------------------------|----------------------------|
| `hic.matrixplot()`           | `hic_matrixplot()`         |
| `construct.experiment()`     | `load_contacts()`          |
| `insulation.score()`         | `insulation_score()`       |
| `insulation.domainogram()`   | `insulation_domainogram()` |
| `fastDI()`                   | `direct_index()`           |
| `chromosomeMatrix()`         | `chromsome_matrix()`       |
| `compartment.score()`        | `compartment_score()`      |
| `insulation.callTAD()`       | `call_TAD_insulation()`    |
| `intra.inter.TAD.contacts()` | `intra_inter_TAD()`        |
| `select.subset()`            | `select_subset()`          |
| `trans.compartment.plot()`   | `trans_matrixplot()`       |


## Deprecated

- The `badBin.find()` function.
- The `badBin.plot()` function.
- The `rescale()` function in favour of `scales::rescale()`.
- The `resize_mat()` function.
- The `differential.TAD.dotplot/scatterplot()` function in favour of the 
  `visualise.IIT_discovery()` method.
- The `frequency.from.matrix()` function.
- The `chrom.comparison.plot()` function in favour of the
  `visualise.chrommat_discovery()` method.
- The `quantifyAPA()` function in favour of the `quantify.APA_discovery()` 
  method.
- The `cisTotal.perChrom()` and `cisTotal.scores()` functions in favour of the
  `cis_trans()` function.
- The `insulation.plot.dual()` and `insulation.plot.single()` functions in 
  favour of the `insulation_matrixplot()`.
  
## Miscelaneous

- Added the `get_test_data()` function and 40k/150k datasets of two small
  chromosomes to quickly try some functions with some data. It uses the 
  {pkgfilecache}

- Limit the number of threads used by {data.table} to 1 by default. We noticed
  that this is faster if functions are executed in a loop.
- Removed {reshape2} and {dplyr} dependencies.
- Replaced {bigwrig} dependency to {rtracklayer} for importing bigwig data.
- Added the `expnames()` function to glance which samples were involved in
  a <discovery> or <contacts> object.
- Added the `resolution()` function to glance at which resolution a <discovery>
  object was calculated or <contacts> object has.

# GENOVA 0.9.995 - 20-3-2019 (The Old Lighthouse)

## Added

- HiCseg.callTAD: a function to call TADs with HiCseq (Levy-Leduc) per chromosome-arm or windowed chromosomes.
- construct.experiment: added an option to get Z-score normalised experiment-objects.
- hic.matrixplot: added an option to use Z-score normalised experiment-objects.
- insulation.heatmap: added option to set leftlost bin to 0
- insulation.heatmap: added option to use a borders-file (also for peaks!)

## Changed

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

# GENOVA 0.9.98 - 27-02-2019
Github-hash: 9951098

## Added

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

## Changed

- HiC_matrixplot: `bed.col` and `bw.col` are merged to chip.col. This argument can take a vector of four colours, parallel to the `chip`-argument.
- HiC_matrixplot: `yMax` is renamed to `chip.ymax`.
- visualise.PESCAn.ggplot: works woth new output of PE-SCAn
- PE-SCAn: cov2d loops over a vectorised list without duplicates, which speeds up the code.
- PE-SCAn: cleaned up cov2d's return-code
- quantifyAPA: the `enrichmentType` argument lets you choose between pixel/mean(backgroundRegions) and the fraction of loops with more than 50% higher signal than the background.
- visualise.compartmentStrength: at least one bin is used (bugfix). 
- saddle: extra effort to use the CS is taken by matching the CS-vector and the eigen-vector. This will fix problems due to big centromeres or shaky/different centromere-calls.
- vignette: fixed typo's

# GENOVA 0.9.971 - 13-06-2018
Github-hash: 8edef82

# GENOVA 0.9.97 - 11-06-2018
Github-hash: 35772ec

# GENOVA 0.9.96 - 22-05-2018
Github-hash: 7deebca



