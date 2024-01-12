library(bsseq)
library(DSS)

setwd("~/p-sy58-0/DMR_scripts/data")

covariates <- read.csv('covariates.csv')
genetic_PCs <- read.csv('top10_PC_scores.csv')
#bs.chrX <- readRDS('bs_chrX_nosnps_min5_percent80.Rds')
bs.chrX <- readRDS('bs_chrX_nosnps_min5_percent80_rm05.Rds')
#bs.chrX <- readRDS('bs_chrX_nosnps_min5_percent80_rm10.Rds')

sample_match <- match(colnames(bs.chrX), covariates[,1])
matched_covariates <- covariates[sample_match,]

genetic_PC_samples <- sapply(strsplit(colnames(bs.chrX), '_'), function(sample_split) paste(sample_split[1]))
genetic_PC_match <- pmatch(genetic_PC_samples, genetic_PCs[,1], duplicates.ok=T)
matched_genetic_PCs <- genetic_PCs[genetic_PC_match, 2:11]

# substitute mean conversion rate if missing for sample
if(sum(is.na(matched_covariates[,12]))>=1) {
  matched_covariates[is.na(matched_covariates[,12]),12] <- mean(na.omit(matched_covariates[,12]))
}

design <- data.frame(as.factor(matched_covariates[,2]), as.factor(matched_covariates[,3]), as.factor(matched_covariates[,4]), as.factor(matched_covariates[,10]), as.factor(matched_covariates[,6]), as.factor(matched_covariates[,11]), as.factor(matched_covariates[,9]), as.numeric(matched_covariates[,12]), matched_genetic_PCs)

colnames(design) <- c('diagnosis', 'cell_type', 'sex', 'age', 'bank', 'pmi', 'hemi', 'conv', paste0('pc', 1:10))

#X <- model.matrix(~diagnosis+cell_type+sex+age+bank+pmi+hemi+conv+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, design)
#dim(X)

DMLfit = DMLfit.multiFactor(bs.chrX, design, ~diagnosis+cell_type+sex+age+bank+pmi+hemi+conv+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+cell_type:sex)
colnames(DMLfit$X)

cell_DMLtest = DMLtest.multiFactor(DMLfit, coef='cell_typeOLIG2')
sex_DMLtest = DMLtest.multiFactor(DMLfit, coef='sexM')
interaction_DMLtest = DMLtest.multiFactor(DMLfit, coef='cell_typeOLIG2:sexM')

write.csv(cell_DMLtest, 'DSS_cell_raw_gPC_rm05.csv', quote=F, row.names=F)
write.csv(sex_DMLtest, 'DSS_sex_raw_gPC_rm05.csv', quote=F, row.names=F)
write.csv(interaction_DMLtest, 'DSS_sex_cell_interaction_raw_gPC_rm05.csv', quote=F, row.names=F)

head(cell_DMLtest)

# calculate n for Bonferroni correction
num_samples <- nrow(cell_DMLtest) - length(which(is.na(cell_DMLtest[,4])))

cell_dmp <- callDML(cell_DMLtest, delta=0, p.threshold=(0.05 / num_samples))
sex_dmp <- callDML(sex_DMLtest, delta=0, p.threshold=(0.05 / num_samples))
interaction_dmp <- callDML(interaction_DMLtest, delta=0, p.threshold=(0.05 / num_samples))

write.table(cell_dmp, 'DSS_cell_raw_gPC_rm05_DMPs.txt', quote=F, row.names=F)
write.table(sex_dmp, 'DSS_sex_raw_gPC__rm05_DMPs.txt', quote=F, row.names=F)
write.table(interaction_dmp, 'DSS_sex_cell_interaction_raw_gPC_rm05_DMPs.txt', quote=F, row.names=F)

cell_dmr <- callDMR(cell_DMLtest, p.threshold=(0.05 / num_samples), minCG=5, dis.merge=100)
sex_dmr <- callDMR(sex_DMLtest, p.threshold=(0.05 / num_samples), minCG=5, dis.merge=100)
interaction_dmr <- callDMR(interaction_DMLtest, p.threshold=(0.05 / num_samples), minCG=5, dis.merge=100)

write.table(cell_dmr, 'DSS_cell_raw_gPC_rm05_DMRs.txt', quote=F, row.names=F)
write.table(sex_dmr, 'DSS_sex_raw_gPC_rm05_DMRs.txt', quote=F, row.names=F)
write.table(interaction_dmr, 'DSS_sex_cell_interaction_raw_gPC_rm05_DMRs.txt', quote=F, row.names=F)

