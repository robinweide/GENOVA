#' @keywords internal
"_PACKAGE"

#' @import utils graphics stats grDevices
#' @import methods
NULL

# data.table is generally careful to minimize the scope for namespace
# conflicts (i.e., functions with the same name as in other packages);
# a more conservative approach using @importFrom should be careful to
# import any needed data.table special symbols as well, e.g., if you
# run DT[ , .N, by='grp'] in your package, you'll need to add
# @importFrom data.table .N to prevent the NOTE from R CMD check.
# See ?data.table::`special-symbols` for the list of such symbols
# data.table defines; see the 'Importing data.table' vignette for more
# advice (vignette('datatable-importing', 'data.table')).
#
#' @import data.table
NULL


utils::globalVariables(c('.',
                         LETTERS[1:5],
                         'colour', 
                         'samplename',
                         'SUM', 
                         'DI',
                         'distance',
                         'value', 
                         'variable',
                         'contacts',
                         'facet',
                         'Var3',
                         'pixel',
                         'N',
                         'background',
                         'experiment',
                         'rowIDX',
                         'score',
                         'mid', 
                         'P',
                         'pos',
                         '..contrast',
                         '..grab',
                         'left',
                         'right',
                         'delta',
                         'ends',
                         'starts',
                         'values',
                         'ins',
                         'type',
                         'strength',
                         'feat',
                         'antoni',
                         'mu',
                         'ord',
                         'q1',
                         'q2',
                         'id',
                         'idx1',
                         'idx2',
                         'uni5p',
                         'i',
                         'j',
                         'term1', 'term2',
                         paste0('V',1:6),
                         "x", "xmin", "xmax",
                         'y',
                         'J', 
                         'unLog',
                         'experiment',
                         'ymin',
                         '..expnames', 
                         'position',
                         'insulation',
                         'chromosome',
                         "part",
                         'rcp',
                         'panel',
                         'dir_index',
                         "start",
                         'index',
                         'Var1',
                         'Var2',
                         'obs.exp.matrix',
                         'obsexp', 
                         'region', 
                         'signal',
                         'newbin',
                         'i.V3',
                         'V3.x', 'V3.y', 'n', '.x',
                         'chrom', 
                         'chrom_x',
                         'chrom_y',
                         'chr', 
                         'grp',
                         'error',
                         'bin',
                         'binned'))

