suppressPackageStartupMessages(library(bsseq))
library(dplyr)
library(matrixStats)

setwd("~/p-sy58-0/plot_scripts")

covariates <- read.csv('../covariates.csv')
covariates_WGBS <- covariates %>% filter(WGBS_Conversion_rates >= 0.95)

samples_WGBS <- covariates_WGBS[['Sample']]

snps_df <- read.table(file='../1000G_phase1.chr5_minAF_0.01_BS_snps.hg38.tsv', sep='\t', header=T)
snps_gr <- makeGRangesFromDataFrame(snps_df, 
                                    seqnames.field='CHROM', 
                                    start.field='POS', 
                                    end.field='POS',
                                    ignore.strand=T)

bs.chr <- readRDS('bs_chrX.Rds')
seqlevels(bs.chr) <- 'chrX'

M <- getCoverage(bs.chr, type='M')
Cov <- getCoverage(bs.chr, type='Cov')
loci_gr <- getBSseq(bs.chr, type='gr')

sample_mask <- colnames(M) %in% samples_WGBS

snp_mask <- is.na(findOverlaps(loci_gr, snps_gr, select='arbitrary', ignore.strand=T))

loci_gr_nosnps <- loci_gr[snp_mask]
M_nosnps <- M[snp_mask, sample_mask]
Cov_nosnps <- Cov[snp_mask, sample_mask]

bs.chr.nosnps <- BSseq(gr=loci_gr_nosnps, 
                        M=M_nosnps, 
                        Cov=Cov_nosnps, 
                        sampleNames=samples_WGBS)
saveRDS(bs.chr.nosnps, file='bs_chrX_nosnps.Rds')
