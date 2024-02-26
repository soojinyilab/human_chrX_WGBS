
##### Methylation Calculations for Promoters and Genebodies #####
l=ls()
l=ls()
l=l[which(!ls() %in% c("","",""))]
rm(list=l)

# read in samples and their covariates
metatab<-read.table("samples_covariates.txt",header=T,sep="\t",stringsAsFactors=F)
rownames(metatab)=metatab$Sample
metatab=metatab[grep("WGBS",metatab$Analyses),]
metatab$group=paste(metatab$Cell,metatab$Sex,sep="_")
metatab=metatab[order(metatab$group),]
table(metatab$group)
head(metatab)

# read in methylation data (ie. BSseq object)
mettab<-readRDS("bs_chrX.Rds")
cpgpositions=paste(mettab@rowRanges@seqnames,mettab@rowRanges@ranges@start,sep="_")
samplenames=mettab@colData@rownames
table(samplenames %in% rownames(metatab))

metcounttab=mettab@assays@data$M   #methylation counts table
colnames(metcounttab)=samplenames
rownames(metcounttab)=cpgpositions
head(metcounttab)

metcovtab=mettab@assays@data$Cov   #methylation coverage table
colnames(metcovtab)=samplenames
rownames(metcovtab)=cpgpositions
head(metcovtab)
colnames(metcovtab)==colnames(metcounttab)

metfractab=metcounttab/metcovtab   #methylation fraction table
head(metfractab)
dim(metfractab)
colnames(metfractab)==colnames(metcounttab)

#for all groups, minimum coverage of 5 in at least 80% of the group's samples 
covcut=5
grpperccut=0.8
metcovpasstab=(metcovtab>=covcut)
grp1samples=metatab[metatab$group=="NeuN_F","Sample"]
grp2samples=metatab[metatab$group=="NeuN_M","Sample"]
grp3samples=metatab[metatab$group=="OLIG2_F","Sample"]
grp4samples=metatab[metatab$group=="OLIG2_M","Sample"]
grp1length=length(grp1samples)
grp2length=length(grp2samples)
grp3length=length(grp3samples)
grp4length=length(grp4samples)

head(metcovpasstab)
grp1pass=(apply(metcovpasstab[,grp1samples],1,sum)>(grpperccut*grp1length))
grp2pass=(apply(metcovpasstab[,grp2samples],1,sum)>(grpperccut*grp2length))
grp3pass=(apply(metcovpasstab[,grp3samples],1,sum)>(grpperccut*grp3length))
grp4pass=(apply(metcovpasstab[,grp4samples],1,sum)>(grpperccut*grp4length))
pass=((grp1pass)&(grp2pass)&(grp3pass)&(grp4pass))
pass
table(pass)
table(names(pass)==row.names(metcounttab))

metfractab=metfractab[pass,]
head(metfractab)
dim(metfractab)

#remove cpg positions of known SNPs
snptab<-read.table("1000G_phase1.chrX_minAF_0.01_BS_snps.hg38.tsv",header=T,sep="\t",stringsAsFactors=F)
snptab$chrpos=paste(sub("chr","",snptab$CHROM),snptab$POS,sep="_")
head(snptab)
summary(snptab$AF)

table(snptab$chrpos %in% row.names(metfractab))
snptab[(snptab$chrpos %in% row.names(metfractab)),]
table(row.names(metfractab) %in% snptab$chrpos)
metfractab=metfractab[!(row.names(metfractab) %in% snptab$chrpos),]
metfracdf=data.frame(metfractab,check.names=F)
metfracdf$chr=sub("_\\S+","",row.names(metfracdf))
metfracdf$pos=as.numeric(sub("^\\S+\\_","",row.names(metfracdf)))
head(metfracdf)
dim(metfracdf)

### Calculate Gene Level methylation data
l=ls()
l=ls()
l=l[which(!ls() %in% c("metfractab","metatab","metfracdf"))]
rm(list=l)
head(metfracdf)

# read in gene info file, originally generated from gtf
genetab<-read.table("gene_maptab.txt",header=T,sep="\t",stringsAsFactors=F,check.names=F)
genetab=genetab[genetab$chr=="X",]
genetab=genetab[(genetab$gene_biotype=="protein_coding")|(genetab$gene_name=="XIST"),]
head(genetab)

NFsamples=metatab[metatab$group=="NeuN_F","Sample"]
NMsamples=metatab[metatab$group=="NeuN_M","Sample"]
OFsamples=metatab[metatab$group=="OLIG2_F","Sample"]
OMsamples=metatab[metatab$group=="OLIG2_M","Sample"]
Fsamples=c(NFsamples,OFsamples)
Msamples=c(NMsamples,OMsamples)

# genebody calculations
gb_names=c("DMgenebody")
genegenebodytab=genetab
head(genegenebodytab)
for (k in 1:length(gb_names)) {
  #k=1
  gb_name=gb_names[k]
  genegenebodytab[genegenebodytab$strand=="+","start"]
  genegenebodytab[,paste(gb_name,"_start",sep="")]=genetab[,"start"]
  genegenebodytab[,paste(gb_name,"_end",sep="")]=genetab[,"end"]
  head(genegenebodytab)
  
  for (i in 1:length(row.names(genegenebodytab))) {
    #i=1
    head(metfracdf)
    tempstart=genegenebodytab[i,paste(gb_name,"_start",sep="")]
    tempend=genegenebodytab[i,paste(gb_name,"_end",sep="")]
    tempmetfrac=metfracdf[(metfracdf$pos>=tempstart)&(metfracdf$pos<=tempend),]
    head(tempmetfrac)
    
    if (length(row.names(tempmetfrac))>0) {
      tempAF=apply(tempmetfrac[,Fsamples],2,mean,na.rm=T)
      tempAM=apply(tempmetfrac[,Msamples],2,mean,na.rm=T)
      tempNF=apply(tempmetfrac[,NFsamples],2,mean,na.rm=T)
      tempNM=apply(tempmetfrac[,NMsamples],2,mean,na.rm=T)
      tempOF=apply(tempmetfrac[,OFsamples],2,mean,na.rm=T)
      tempOM=apply(tempmetfrac[,OMsamples],2,mean,na.rm=T)
      
      listA=list('FvM'=c(tempAF),
                 'NFvNM'=c(tempNF),
                 'OFvOM'=c(tempOF))
      listB=list('FvM'=c(tempAM),
                 'NFvNM'=c(tempNM),
                 'OFvOM'=c(tempOM))
      
      for (j in 1:length(listA)) {
        #j=1   j=2
        temptest=wilcox.test(listA[[j]],listB[[j]],digits.rank = 7)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_pval",sep="")]=signif(temptest$p.value,3)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_mean",sub("v[^v]+$","",names(listA)[j]),sep="")]=signif(mean(listA[[j]],na.rm=T),3)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_mean",sub("^[^v]+v","",names(listA)[j]),sep="")]=signif(mean(listB[[j]],na.rm=T),3)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_log2FC",sep="")]=signif(log2(mean(listA[[j]],na.rm=T)/mean(listB[[j]],na.rm=T)),3)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_diff",sep="")]=signif(mean(listA[[j]],na.rm=T)-mean(listB[[j]],na.rm=T),3)
        genegenebodytab[i,paste(gb_name,"_",names(listA)[j],"_normdiff",sep="")]=signif((mean(listA[[j]],na.rm=T)-mean(listB[[j]],na.rm=T))/((mean(listA[[j]],na.rm=T)+mean(listB[[j]],na.rm=T))),3)
      }
    }
    if (i%%100==0) {
      cat(i," ",sep="")
    }
  }
}
head(genegenebodytab)
write.table(genegenebodytab,"genebody_methylation.txt",quote=F,sep="\t",row.names=T,col.names=T)

# genepromoter calculations
promot_names=c("DMp0.5k")
promot_uplens=c(500)
promot_downlens=c(500)

genepromotertab=genetab
head(genepromotertab)
for (k in 1:length(promot_names)) {
  #k=1
  p_name=promot_names[k]
  p_up=promot_uplens[k]
  p_down=promot_downlens[k]
  genepromotertab[genepromotertab$strand=="+","start"]
  genepromotertab[genetab$strand=="+",paste(p_name,"_start",sep="")]=genetab[genetab$strand=="+","start"]-p_up
  genepromotertab[genetab$strand=="+",paste(p_name,"_end",sep="")]=genetab[genetab$strand=="+","start"]+p_down
  genepromotertab[genetab$strand=="-",paste(p_name,"_start",sep="")]=genetab[genetab$strand=="-","end"]-p_down
  genepromotertab[genetab$strand=="-",paste(p_name,"_end",sep="")]=genetab[genetab$strand=="-","end"]+p_up
  head(genepromotertab)
  
  for (i in 1:length(row.names(genepromotertab))) {
    #i=1
    head(metfracdf)
    tempstart=genepromotertab[i,paste(p_name,"_start",sep="")]
    tempend=genepromotertab[i,paste(p_name,"_end",sep="")]
    tempmetfrac=metfracdf[(metfracdf$pos>=tempstart)&(metfracdf$pos<=tempend),]
    head(tempmetfrac)
    
    if (length(row.names(tempmetfrac))>0) {
      tempAF=apply(tempmetfrac[,Fsamples],2,mean,na.rm=T)
      tempAM=apply(tempmetfrac[,Msamples],2,mean,na.rm=T)
      tempNF=apply(tempmetfrac[,NFsamples],2,mean,na.rm=T)
      tempNM=apply(tempmetfrac[,NMsamples],2,mean,na.rm=T)
      tempOF=apply(tempmetfrac[,OFsamples],2,mean,na.rm=T)
      tempOM=apply(tempmetfrac[,OMsamples],2,mean,na.rm=T)
      
      listA=list('FvM'=c(tempAF),
                 'NFvNM'=c(tempNF),
                 'OFvOM'=c(tempOF))
      listB=list('FvM'=c(tempAM),
                 'NFvNM'=c(tempNM),
                 'OFvOM'=c(tempOM))
      
      for (j in 1:length(listA)) {
        #j=1   j=2
        temptest=wilcox.test(listA[[j]],listB[[j]],digits.rank = 7)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_pval",sep="")]=signif(temptest$p.value,3)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_mean",sub("v[^v]+$","",names(listA)[j]),sep="")]=signif(mean(listA[[j]],na.rm=T),3)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_mean",sub("^[^v]+v","",names(listA)[j]),sep="")]=signif(mean(listB[[j]],na.rm=T),3)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_log2FC",sep="")]=signif(log2(mean(listA[[j]],na.rm=T)/mean(listB[[j]],na.rm=T)),3)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_diff",sep="")]=signif(mean(listA[[j]],na.rm=T)-mean(listB[[j]],na.rm=T),3)
        genepromotertab[i,paste(p_name,"_",names(listA)[j],"_normdiff",sep="")]=signif((mean(listA[[j]],na.rm=T)-mean(listB[[j]],na.rm=T))/((mean(listA[[j]],na.rm=T)+mean(listB[[j]],na.rm=T))),3)
      }
    }
    if (i%%100==0) {
      cat(i," ",sep="")
    }
  }
}
head(genepromotertab)
write.table(genepromotertab,"genepromoter_methylation.txt",quote=F,sep="\t",row.names=T,col.names=T)

##### Methylation Calculations for Promoters and Genebodies #####
