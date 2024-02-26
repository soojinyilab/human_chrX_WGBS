
##### RNAseq - DESeq #####
l=ls()
l=ls()
l=l[which(!ls() %in% c("","",""))]
rm(list=l)

genemaptab<-read.table("gene_maptab.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F)
head(genemaptab)
genecounttab<-read.table("gene_counttab_filtered.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F)
table(row.names(genecounttab) %in% (genemaptab$gene_id))
genecounttab=genecounttab[(row.names(genecounttab) %in% genemaptab[genemaptab$chr=="X","gene_id"]),]
head(genecounttab)
names(genecounttab)

metatab<-read.table("de_metatab.txt",header=T,sep="\t",stringsAsFactors=F)
row.names(metatab)=metatab$id
metatab$celltype_sex=as.factor(paste(metatab$Cell_type,metatab$sex,sep="_"))
metatab
metatab$id %in% names(genecounttab)

runlist=list(gene=genecounttab)
library(DESeq2)

for (i in 1:length(runlist)) {
  #i=1   i=2
  runname=names(runlist)[i]
  counttab=runlist[[i]]
  
  mcounttab=as.matrix(counttab[,metatab$id])
  mcounttab=mcounttab[apply(mcounttab,1,sum)!=0,]
  head(mcounttab)
  dim(mcounttab)
  table(colnames(mcounttab)==row.names(metatab))
  
  # DESeq2
  metatab$Cell_type
  metatab$sex
  
  dds <- DESeq2::DESeqDataSetFromMatrix(mcounttab, metatab, formula(~ Cell_type + sex + Cell_type:sex))
  design(dds)
  
  # get the model matrix
  mod_mat <- model.matrix(design(dds), colData(dds))
  mod_mat
  
  # Define coefficient vectors for each condition
  NeuN_F <- colMeans(mod_mat[dds$Cell_type == "NeuN" & dds$sex == "F", ])
  NeuN_M <- colMeans(mod_mat[dds$Cell_type == "NeuN" & dds$sex == "M", ])
  OLIG2_F <- colMeans(mod_mat[dds$Cell_type == "OLIG2" & dds$sex == "F", ])
  OLIG2_M <- colMeans(mod_mat[dds$Cell_type == "OLIG2" & dds$sex == "M", ])
  
  # run
  dds <- DESeq(dds)
  DESeq2::plotDispEsts(dds)
  resultsNames(dds)
  
  reslist=list( OLIG2vNeuN.All=data.frame(results(dds, contrast = (OLIG2_F + OLIG2_M) - (NeuN_F + NeuN_M))),
                OLIG2vNeuN.Female=data.frame(results(dds, contrast = (OLIG2_F - NeuN_F))),
                OLIG2vNeuN.Male=data.frame(results(dds, contrast = (OLIG2_M - NeuN_M))),
                FvM.All=data.frame(results(dds, contrast = (OLIG2_F + NeuN_F) - (OLIG2_M + NeuN_M))),
                FvM.OLIG2=data.frame(results(dds, contrast = (OLIG2_F - OLIG2_M))),
                FvM.NeuN=data.frame(results(dds, contrast = (NeuN_F - NeuN_M))),
                OLIG2vNeuN.x.FvM=data.frame(results(dds, contrast = (OLIG2_F - NeuN_F) - (OLIG2_M - NeuN_M)))
  )
  
  desummary=data.frame()
  for (j in 1:length(reslist)) {
    #j=1
    compname=names(reslist)[j]
    tempresult=reslist[[j]]
    head(tempresult)
    head(genemaptab)
    if (runname=="gene") {
      tempresult2=merge(tempresult,genemaptab[,c("gene_id","chr","start","end","strand","gene_name","gene_biotype")],by.x=0,,by.y="gene_id",all.x=T)
      #row.names(tempresult2)=tempresult2$Row.names
      names(tempresult2)[1]="gene_id"
    }else if (runname=="transcript") {
      tempresult2=merge(tempresult,txmaptab[,c("transcript_id","chr","start","end","strand","gene_id","gene_name","gene_biotype","transcript_name")],by.x=0,,by.y="transcript_id",all.x=T)
      #row.names(tempresult2)=tempresult2$Row.names
      names(tempresult2)[1]="transcript_id"
    }
    head(tempresult2)
    tempresult2=tempresult2[order(tempresult2$pvalue,decreasing=F),]
    head(tempresult2)
    #write.table(tempresult2,paste("deseq2.",runname,".",compname,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=T)
    desummary[compname,"num_features_tested"]=length(row.names(tempresult2))
    desummary[compname,"num_features_pval<0.05"]=sum(tempresult2$pvalue<0.05,na.rm=T)
    desummary[compname,"num_features_padj<0.05"]=sum(tempresult2$padj<0.05,na.rm=T)

    cat(runname,"_",compname,"\n",sep="")
  }
  desummary
  #write.table(desummary,paste("deseq2.",runname,".summary.txt",sep=""),quote=F,sep="\t",row.names=T,col.names=T)
}

##### RNAseq - DESeq #####
