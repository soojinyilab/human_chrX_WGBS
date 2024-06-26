library(bsseq)
library(dplyr)
library(matrixStats)

library(ggplot2)

setwd('~/p-sy58-0/DMR_scripts/data')

covariates <- read.csv('covariates.csv')
covariates_WGBS <- covariates %>% filter(WGBS_Conversion_rates >= 0.95)
samples_WGBS <- covariates_WGBS[['Sample']]

idx_F <- which(covariates_WGBS$Sex == 'F')
idx_NeuN <- which(covariates_WGBS$Cell == 'NeuN')

bs.chrX.nosnps <- readRDS('bs_chrX_nosnps.Rds')

Cov_nosnps <- getCoverage(bs.chrX.nosnps, type='Cov')
M_nosnps <- getCoverage(bs.chrX.nosnps, type='M')
meth_nosnps <- getMeth(bs.chrX.nosnps, type='raw', what='perBase')
loci_gr_nosnps <- getBSseq(bs.chrX.nosnps, type='gr')

mean_meth <- rowMeans(meth_nosnps)

min_cov <- 5
min_rep <- 0.8
num_samples <- length(samples_WGBS)
num_F <- length(idx_F)
num_NeuN <- length(idx_NeuN)
cov_mask_FNeuN <- rowSums(Cov_nosnps[,(idx_F & idx_NeuN)] >= min_cov) >= round(num_F * min_rep, 0)
cov_mask_MNeuN <- rowSums(Cov_nosnps[,(-idx_F & idx_NeuN)] >= min_cov) >= round((num_samples - num_F) * min_rep, 0)
cov_mask_FOlig2 <- rowSums(Cov_nosnps[,(idx_F & -idx_NeuN)] >= min_cov) >= round(num_F * min_rep, 0)
cov_mask_MOlig2 <- rowSums(Cov_nosnps[,(-idx_F & -idx_NeuN)] >= min_cov) >= round((num_samples - num_F) * min_rep, 0)
cov_mask <- which(cov_mask_FNeuN & cov_mask_MNeuN & cov_mask_FOlig2 & cov_mask_MOlig2)

filter_mask <- cov_mask
loci_gr_nosnps_filtered <- loci_gr_nosnps[filter_mask]
M_nosnps_filtered <- M_nosnps[filter_mask,]
Cov_nosnps_filtered <- Cov_nosnps[filter_mask,]

bs.chrX.nosnps.filtered <- BSseq(gr=loci_gr_nosnps_filtered, 
                                 M=M_nosnps_filtered, 
                                 Cov=Cov_nosnps_filtered, 
                                 sampleNames=samples_WGBS)
saveRDS(bs.chrX.nosnps.filtered, file='bs_chrX_nosnps_min5_percent80.Rds')
frac_meth_nosnps_filtered <- getMeth(bs.chrX.nosnps.filtered, type='raw', what='perBase')
write.table(data.frame('loci'=start(loci_gr_nosnps_filtered), frac_meth_nosnps_filtered), 
            file='chrX_CG_fracMeth_perBase_min5_percent80_loci.txt', 
            quote=F, sep='\t', row.names=F)
