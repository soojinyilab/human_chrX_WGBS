# Load libraries
#suppressPackageStartupMessages(library(sva))
#suppressPackageStartupMessages(library(limma))
#suppressPackageStartupMessages(library(DESeq2))
#suppressPackageStartupMessages(library(scales))
#suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(ggjoy))
#suppressPackageStartupMessages(library(knitr))
#suppressPackageStartupMessages(library(preprocessCore))
#suppressPackageStartupMessages(library(variancePartition))
#suppressPackageStartupMessages(library(doParallel))
#suppressPackageStartupMessages(library(Biobase))
#suppressPackageStartupMessages(library(DMwR))
#suppressPackageStartupMessages(library(DT))
#suppressPackageStartupMessages(library(randomForest))
#suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggpubr))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(plyr))
#suppressPackageStartupMessages(library(xlsx))
#suppressPackageStartupMessages(library(pheatmap))
#suppressPackageStartupMessages(library(edgeR))
#library(splitstackshape)
#library(tidyr)
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(GGally))
#library(circlize)
#library(colorspace)
#library(GetoptLong)
#library("RColorBrewer")
#library(png)
#library(grDevices)
#library("gridExtra")
suppressPackageStartupMessages(library(ggsci))
#library(tidyverse)
#library("gridExtra")
#library(ggjoy)
#library("scatterplot3d")


setwd("~/p-sy58-0/var_explained")
source("Utility_Functions.R")

#meth_mat <- read.table(file="plot_scripts/full_DMRs_smooth_avgMeth_perRegion.txt",
#                        check.names=F, row.names=1, header=T, sep='\t')
#meth_mat <- read.table(file="plot_scripts/chrX_smooth_CG_perBase_min3_mean5_rmExtreme_combined.txt",
#                       check.names=FALSE, row.names=1,header=T, sep='\t')
meth_mat <- read.table(file="~/p-sy58-0/DMR_scripts/data/chrX_smooth_CG_fracMeth_perBase_min5_percent80_loci.txt",
			check.names=F, row.names=1, header=T, sep='\t')
#meth_mat <- subset(meth_mat, select=-c(loci))
colnames(meth_mat) <- sapply(strsplit(colnames(meth_mat), split='X', fixed=T), function(sample_splitX) (sample_splitX[length(sample_splitX)]))
sample_names <- colnames(meth_mat)
meth_mat <- data.matrix(meth_mat)

colnames(meth_mat)
dim(meth_mat)
is.na(meth_mat)
pd <- read.table(file="~/p-sy58-0/covariates.csv", row.names=1, header=T,sep=',')
covar_names <- c('Diagnosis', 'Sex', 'Cell', 'BrainBank', 'Hemisphere', 'AgeClass', 'PMIClass', 'WGBS_Conversion_rates')
pd <- pd[sample_names, covar_names]

###Var explained###
var <- VarExp(meth_mat, pd, 80, FALSE)

var_plot_input <- data.frame(Type="WGBS", eff=names(var), prop=var)
var_plot_input <- var_plot_input[var_plot_input$eff != "resid",]
var_plot_input$eff <- factor(var_plot_input$eff, levels=c("Diagnosis", "Sex", "Cell", "BrainBank", "Hemisphere", "AgeClass", "PMIClass", "WGBS_Conversion_rates"), labels=c("Diagnosis", "Sex", "Cell", "Brain bank", "Hemisphere", "Age", "PMI", "Conversion rate"))
var_plot_input$percent <- var_plot_input$prop * 100
var_plot_input

###Plotting###
var_barplot <- ggplot(var_plot_input, aes(x=reorder(eff, -prop), y=prop)) + 
  geom_bar(stat="identity",colour="black", aes(fill=Type)) +
  #ggtitle("Variance of methylation explained") +
  #coord_cartesian(ylim=c(0,9))+
  ylab("Variance explained (%)") +
  theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
  #geom_text(aes(label=round(percent,2), y=prop+0.04),position = position_dodge(0.9),vjust=1.6, size=4)+
  scale_fill_aaas() +
  theme(
    legend.position = "none",#element_blank(),#c(0.77,0.89),#"inside",#element_blank(),#"bottom",
    plot.background = element_blank(),
    axis.text.x = element_text(color="black", size=11, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 11),
    strip.text = element_text(size = 11),
    text = element_text(size = 11),
    legend.title = element_text(size=11),
    legend.text = element_text(size=11),
    axis.title.x = element_blank(),#element_text(size = 11),
    axis.title.y = element_text(size = 11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

var_barplot

########################
########################
pca.Sample <- prcomp(t(meth_mat))
dd <- data.frame(pca.Sample$x)
head(dd)
dd["Diagnosis"] <- pd$Diagnosis
dd["Sex"] <- pd$Sex
dd["Cell"] <- pd$Cell
dd["BrainBank"] <- pd$BrainBank
dd["Hemisphere"] <- pd$Hemisphere
dd["Age"] <- pd$AgeClass
dd["PMI"] <- pd$PMIClass
dd["BS_conversion"] <- pd$WGBS_Conversion_rates

PCi <- dd
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)

PCi$Sex <- factor(PCi$Sex, levels=c("M", "F"),
                        labels=c("Male", "Female"))

pca_a <- ggscatter(PCi, x = "PC1", y = "PC2",
                  color = "Cell",
                  shape="Sex",
                  size = 2.4) +
  xlab(paste("PC1 (",round(variance[1],1),"% )")) + 
  ylab(paste("PC2 (",round(variance[2],1),"% )")) +
  #ggtitle("Principal component analysis") +
  theme_classic() +
  theme_bw() + theme(
    legend.position='right', 
    plot.background = element_blank(),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 11),
    strip.text = element_text(size = 11),
    text = element_text(size = 11),
    axis.title.x = element_text(size = 11, margin = margin(t=-20)),
    axis.title.y = element_text(size = 11),
    legend.title = element_text(size=11),
    legend.text = element_text(size=11),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pca_a
var_barplot

combined_plot <- plot_grid(pca_a, var_barplot, nrow = 1, rel_widths = c(1.0, .8), labels=c("B", "C"), vjust=1, align="h") +
	theme(plot.background=element_rect(fill="white"), panel.border=F)
combined_plot
save_plot("PCA_and_Var_explained.png", combined_plot, base_width=7, base_height=3.5)
#save_plot("PCA_and_Var_explained.pdf", combined_plot,
#          base_height = 4.2,
#          base_width = 9.1)
