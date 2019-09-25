library(GENOVA); library(tictoc)
IDX = '/DATA/projects/Hap1/Med12/hicpro/5670/output/hic_results/matrix/1_wildtype_haploid/raw/10000/1_wildtype_haploid_10000_abs.bed'

tic('3 samples 50kb load')

WT = loadContacts('/DATA/projects/Hap1/Med12/hicpro/5670/output/hic_results/matrix/1_wildtype_haploid/iced/10000/1_wildtype_haploid_10000_iced.matrix', 
                  indicesPath = IDX, 
                  sampleName = 'WT')

MED12 = loadContacts('/DATA/projects/Hap1/Med12/hicpro/5670/output/hic_results/matrix/2_MED12_3_10/iced/10000/2_MED12_3_10_10000_iced.matrix', 
                     indicesPath = IDX, 
                     sampleName = 'MED12')

DKO = loadContacts('/DATA/projects/Hap1/Med12/hicpro/5670/output/hic_results/matrix/3_MED12_Wapl_3_1/iced/10000/3_MED12_Wapl_3_1_10000_iced.matrix', 
                   indicesPath = IDX, 
                   sampleName = 'DKO')

toc()

library(GenomicRanges)
loops = read.delim('/shared/dewit/oidBackUp/WAPL_Project/Hi-C/analysis/HICCUPS/WT_rep123.hiccups.bedpe_5k10k25k_merged.bedpe')
anchors = loops[,1:3]; anchors$CHR1 = paste0('chr',anchors$CHR1)
anchors = makeGRangesFromDataFrame(anchors,start.field = 'X1',end.field = 'X2', seqnames.field = 'CHR1')
anchors_shuffle = regioneR::circularRandomizeRegions(anchors)
loopList = lapply(list(anchors = anchors, shuffled = anchors_shuffle), as.data.frame)

RCP_out = RCP(explist = list('WT' = WT, 'MED12' = MED12, 'DKO' = DKO), bedlist = loopList)



RCP_out


visualise(RCP_out, metric = 'smooth')
visualise(RCP_out, metric = 'smooth', flipFacet = T)
visualise(RCP_out, metric = 'smooth', raw = T)

visualise(RCP_out, metric = 'both')
visualise(RCP_out, metric = 'both', flipFacet = T)
visualise(RCP_out, metric = 'both', raw = T)

visualise(RCP_out, metric = 'lfc', contrast = 'WT')
visualise(RCP_out, metric = 'lfc', contrast = 'WT', flipFacet = T)
visualise(RCP_out, metric = 'lfc', contrast = 'WT', raw = T)
