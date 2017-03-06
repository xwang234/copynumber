#!/usr/bin/env Rscript
#original positions: GRCh38.p7
#FHIT:chr3:59000000-62300000
#CDKN2A:chr9:21967752-21995043
#WWOX:chr16:78099413-79212667
#chr9p
#hg19,by liftover
#chr3:58985726-62285675
#chr9:21967751-21995042
#chr16:78133310-79246564
#chr9:0-49000000
regions=data.frame(chr=c("3","9","16","9"),start=c(58985726,21967751,78133310,0),end=c(62285675,21995042,79246564,49000000))
rownames(regions)=c("FHIT","CDKN2A","WWOX","Chr9p")








wgstumors1=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
             "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868","SRR1001635")
wgstumors2=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
wgstumors3=paste0("T",c(1:6,8:18))
freecdir1="/fh/scratch/delete30/dai_j/freec"
freecdir2="/fh/scratch/delete30/dai_j/henan/freec"
freecdir3="/fh/scratch/delete30/dai_j/escc/freec"
library(GenomicRanges)

extractregions=function(wgstumors,freecdir)
{
  gr_regions=GRanges(seqnames=regions$chr,ranges=IRanges(start=regions$start,end=regions$end))
  cnres=bafres=ratiores=cbind(regions,data.frame(matrix(NA,nrow=nrow(regions),ncol=length(wgstumors))))
  colnames(ratiores)[(ncol(regions)+1):ncol(ratiores)]=wgstumors
  colnames(cnres)=colnames(bafres)=colnames(ratiores)
  for (i in 1:length(wgstumors))
  {
    ratiofile=paste0(freecdir,"/",wgstumors[i],"/ploid2degree3force0/",wgstumors[i],".pileup.gz_ratio.txt")
    ratiodata=read.table(file=ratiofile,header=T,sep="\t")
    #remove outliars caused by extra high coverage 
    ratiocutoff1=0 #as.numeric(quantile(ratiodata$Ratio[ratiodata$Ratio!=-1],0.001))
    ratiocutoff2=as.numeric(quantile(ratiodata$Ratio[ratiodata$Ratio!=-1],0.999))
    ratiodata=ratiodata[ratiodata$Ratio!=-1 & ratiodata$Ratio >= ratiocutoff1 & ratiodata$Ratio <= ratiocutoff2 & ratiodata$BAF!=-1
                        & ratiodata$UncertaintyOfGT<20 , ]
    gr_ratiodata=GRanges(seqnames=ratiodata$Chromosome,ranges=IRanges(start=ratiodata$Start,width=1000),ratio=ratiodata$Ratio,
                         baf=ratiodata$BAF,uncertain=ratiodata$UncertaintyOfGT,copynumber=ratiodata$CopyNumber)
    for (j in 1:nrow(regions))
    {
      olap=subsetByOverlaps(gr_ratiodata,gr_regions[j,])
      ratios=mcols(olap)["ratio"]
      ratios=ratios[,1]
      ratiores[j,ncol(regions)+i]=median(ratios,na.rm=T)
      bafs=mcols(olap)["baf"]
      bafs=bafs[,1]
      uncertainties=mcols(olap)["uncertain"]
      uncertainties=uncertainties[,1]
      bafres[j,ncol(regions)+i]=median(bafs,na.rm=T)
      copynumbers=mcols(olap)["copynumber"]
      copynumbers=copynumbers[,1]
      cnres[j,ncol(regions)+i]=median(copynumbers,na.rm=T)
    }
  }
  return(result=list(ratiores=ratiores,bafres=bafres,cnres=cnres))
}

res1=extractregions(wgstumors=wgstumors1,freecdir=freecdir1)
res2=extractregions(wgstumors=wgstumors2,freecdir=freecdir2)
res3=extractregions(wgstumors=wgstumors3,freecdir=freecdir3)

for (i in 1:4)
{
  print(sum(res1$ratiores[i,4:ncol(res1$cnres)]<0.85))
}

for (i in 1:4)
{
  print(sum(res2$ratiores[i,4:ncol(res2$cnres)]<0.85))
}
#visualization
library("ComplexHeatmap")
library(circlize)
cols=c("blue","green","red")
cols1=c("black","blue","gray","green", "red")
call_delloh=function(cn,baf,opt=0,cutoff1=0.8,cutoff2=0.75)
{
  res=0
  if (cn<=cutoff1)
  {
    res=1
  }else
  {
    if (opt==1) #check baf
    {
      if (baf>=cutoff2)
      {
        res=1
      }
    }
  }
  return(res)
}

#use copynumber
call_delloh1=function(cn,baf,opt=0,cutoff2=0.75)
{
  res=0
  if (cn<2)
  {
    res=1
  }else
  {
    if (cn==2 & opt==1) #check baf
    {
      if (baf>=cutoff2)
      {
        res=1
      }
    }
  }
  return(res)
}

draw_4regions=function(res)
{
  
  #cnmat=as.matrix(res$ratiores[,4:ncol(res$ratiores)])
  cnmat=as.matrix(res$cnres[,4:ncol(res$ratiores)])
  bafmat=as.matrix(res$bafres[,4:ncol(res$bafres)])
  mat=matrix(NA,nrow=nrow(bafmat),ncol=ncol(bafmat))
  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
      if (i==nrow(cnmat))
      {
        mat[i,j]=call_delloh1(cn=cnmat[i,j],baf=bafmat[i,j],opt=1)
      }else
      {
        mat[i,j]=call_delloh1(cn=cnmat[i,j],baf=bafmat[i,j]) #not consider loh
      }
    }
  }

  rownames(mat)=rownames(res$ratiores)
  colnames(mat)=colnames(res$ratiores)[4:ncol(res$ratiores)]
  print(Heatmap(mat,cluster_rows = FALSE,name="del event",show_column_names=F,show_column_dend=F,
          rect_gp = gpar(col = "white", lty = 1, lwd = 1),
          heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes")))
  )
  return(mat)
}
delres1=draw_4regions(res=res1)
delres2=draw_4regions(res=res2)
delres3=draw_4regions(res=res3)


#include all the data
mat=as.matrix(delres1)
mat=cbind.data.frame(mat,as.matrix(delres2))

Heatmap(mat[c(2,4),],column_dend_height = unit(30, "mm"),
        cluster_rows = T,name="del events",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),
        heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes")))
mat1=matrix(c(7,8,3,8),nrow=2,byrow=T)
fisher.test(mat1)
meancnbaf=function(wgstumors,freecdir)
{
  chrs=c(1:22,"X","Y")
  chrs1=paste0("chr",chrs)
  ratiores=bafres=data.frame(matrix(NA,ncol=length(wgstumors),nrow=length(chrs)+1))
  colnames(ratiores)=colnames(bafres)=wgstumors
  rownames(ratiores)=rownames(bafres)=c(chrs1,"all")
  for (i in 1:length(wgstumors))
  {
    ratiofile=paste0(freecdir,"/",wgstumors[i],"/ploid2degree3force0/",wgstumors[i],".pileup.gz_ratio.txt")
    ratiodata=read.table(file=ratiofile,header=T,sep="\t")
    ratiocutoff1=0 #as.numeric(quantile(ratiodata$Ratio[ratiodata$Ratio!=-1],0.001))
    ratiocutoff2=as.numeric(quantile(ratiodata$Ratio[ratiodata$Ratio!=-1],0.999))
    ratiodata=ratiodata[ratiodata$Ratio!=-1 & ratiodata$Ratio >= ratiocutoff1 & ratiodata$Ratio <= ratiocutoff2 & ratiodata$BAF!=-1
                        & ratiodata$UncertaintyOfGT<20, ]
    for (j in 1:length(chrs))
    {
      ratiodata1=ratiodata[ratiodata$Chromosome==chrs[j],]
      ratiores[j,i]=mean(ratiodata1$Ratio,na.rm=T)
      bafres[j,i]=mean(ratiodata1$BAF,na.rm=T)
    }
    ratiores[nrow(ratiores),i]=mean(ratiodata$Ratio,na.rm=T)
    bafres[nrow(bafres),i]=mean(ratiodata$BAF,na.rm=T)
  }
  return(result=list(ratiores=ratiores,bafres=bafres))
}

meanres1=meancnbaf(wgstumors=wgstumors1,freecdir=freecdir1)
meanres2=meancnbaf(wgstumors=wgstumors2,freecdir=freecdir2)
meanres3=meancnbaf(wgstumors=wgstumors3,freecdir=freecdir3)

cutoff=4.589971e-60 #from compute_cnvlength.R
extractregions1=function(wgstumors,freecdir,pcutoff)
{
  gr_regions=GRanges(seqnames=regions$chr,ranges=IRanges(start=regions$start,end=regions$end))
  cnres=widthres=cbind(regions,data.frame(matrix(NA,nrow=nrow(regions),ncol=length(wgstumors))))
  colnames(cnres)[(ncol(regions)+1):ncol(cnres)]=wgstumors
  colnames(widthres)=colnames(cnres)
  for (i in 1:length(wgstumors))
  {
    ratiofile=paste0(freecdir,"/",wgstumors[i],"/ploid2degree3force0/",wgstumors[i],".pvalue.txt")
    ratiodata=read.table(file=ratiofile,header=T,sep="\t")
    
    ratiodata=ratiodata[ratiodata[,10]<=pcutoff & ratiodata$uncertainty<20, ]
    gr_ratiodata=GRanges(seqnames=ratiodata$chr,ranges=IRanges(start=ratiodata$start,end=ratiodata$end),copynumber=ratiodata$copy.number)
    for (j in 1:nrow(regions))
    {
      olap=subsetByOverlaps(gr_ratiodata,gr_regions[j,])
      cns=mcols(olap)["copynumber"]
      cns=cns[,1]
      cnres[j,ncol(regions)+i]=median(cns,na.rm=T)
      widthres[j,ncol(regions)+i]=sum(width(olap))
    }
  }
  return(result=list(cnres=cnres,widthres=widthres))
}

res11=extractregions1(wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
res22=extractregions1(wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)

draw_4regions1=function(res)
{
  
  #cnmat=as.matrix(res$ratiores[,4:ncol(res$ratiores)])
  cnmat=as.matrix(res$cnres[,4:ncol(res$cnres)])
  mat=matrix(0,nrow=4,ncol=ncol(cnmat))
  rownames(mat)=rownames(cnmat)
  colnames(mat)=colnames(cnmat)
  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
      if (i<4)
      {
        if (!is.na(cnmat[i,j]) & cnmat[i,j]<2) mat[i,j]=1
      }else #chr9p
      {
        if (!is.na(cnmat[i,j]) & cnmat[i,j]<=2) mat[i,j]=1
      }
    }
  }

  print(Heatmap(mat,cluster_rows = FALSE,name="del event",show_column_names=T,show_column_dend=F,
                rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes")))
  )
  return(mat)
}

delres11=draw_4regions1(res=res11)
delres22=draw_4regions1(res=res22)
#include all the data
mat=as.matrix(delres11)
mat=cbind.data.frame(mat,as.matrix(delres22))
colnames(mat)=c(paste0("US-EA",1:16),paste0("CH-EA",1:10))

Heatmap(mat,column_dend_height = unit(30, "mm"),
        cluster_rows = T,name="del events",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),
        heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes")))
