library(dplyr)
library(bsseq)
library(BiocParallel)

covariates <- read.csv('~/p-sy58-0/covariates.csv')
fnames <- list.files('~/p-sy58-0/remapping/meth_extract/human_XY/CpG_coverage_report_files', pattern='.txt', full.names=T)

covariates_WGBS <- covariates %>% filter(grepl('WGBS', Analyses))

hg38_XY_mask <- Biostrings::readDNAStringSet('~/p-sy58-0/genomes/human_XY/hg38_XYmask.fa')
chrX_loci <- findLoci('CG', hg38_XY_mask, include=c('X'))

bs.chrX <- read.bismark(fnames,
			loci=chrX_loci,
			rmZeroCov=T, 
			strandCollapse=T)
sampleNames(bs.chrX) <- covariates_WGBS[['Sample']]
sampleNames(bs.chrX)

saveRDS(bs.chrX, 'bs_chrX.Rds')
