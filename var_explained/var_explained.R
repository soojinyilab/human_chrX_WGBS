# Load libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggsci))

setwd("~/p-sy58-0/var_explained")
source("utility_functions.R")

meth_mat <- read.table(file="~/p-sy58-0/DMR_scripts/data/chrX_smooth_CG_fracMeth_perBase_min5_percent80_loci.txt",
			check.names=F, row.names=1, header=T, sep='\t')
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
  ylab("Variance explained (%)") +
  theme_classic() + 
  scale_y_continuous(labels = scales::percent) +
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

combined_plot <- plot_grid(pca_a, var_barplot, nrow = 1, rel_widths = c(1.0, .8), labels=c("B", "C"), vjust=1, align="h") +
	theme(plot.background=element_rect(fill="white"), panel.border=F)
save_plot("PCA_and_Var_explained.png", combined_plot, base_width=7, base_height=3.5)
