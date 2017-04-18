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
  regions$chr=as.character(regions$chr)
  regions$chr=gsub("chr","",regions$chr)
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


#consider the overlap of gistic regions in two nature papers
library(gdata)

readgisticcna=function(gisticfile,sheet=1)
{
  if(grepl(".xls",gisticfile))
  {
    gistictable=read.xls(xls=gisticfile,sheet=sheet,header=F)
  }else
  {
    gistictable=read.table(gisticfile,header=F,sep="\t")
  }
  idxNA=is.na(gistictable[1,])
  gistictable=gistictable[,!idxNA]
  res=data.frame(matrix(NA,nrow=ncol(gistictable)-1,ncol=8))
  colnames(res)=c('chr','start','end','cytoband','qvalue','residualq','numgenes','genes')
  for (i in 2:ncol(gistictable))
  {
    tmp=as.character(gistictable[4,i])
    tmp1=unlist(strsplit(tmp,":"))
    tmp2=unlist(strsplit(tmp1[2],"-"))
    res[i-1,1]=tmp1[1]
    res[i-1,2]=tmp2[1]
    res[i-1,3]=tmp2[2]
    res[i-1,4]=as.character(gistictable[1,i])
    res[i-1,5]=as.character(gistictable[2,i])
    res[i-1,6]=as.character(gistictable[3,i])
    genes=as.character(gistictable[5,i])
    n=1
    for (j in 6:nrow(gistictable))
    {
      gene=as.character(gistictable[j,i])
      if (gene!="")
      {
        n=n+1
        genes=paste0(genes,',',gene)
      }
    }
    res[i-1,7]=n
    res[i-1,8]=genes
  }
  #output=paste0(gisticdir,'/',prefix,'.gistic.cna.txt')
  return(res)
}

#only process the region column
readgisticcna1=function(gisticfile)
{
  
  gistictable=read.table(gisticfile,header=T,sep="\t")
  
  idxNA=is.na(gistictable[1,])
  gistictable=gistictable[,!idxNA]
  res=data.frame(matrix(NA,nrow=nrow(gistictable),ncol=8))
  colnames(res)=c('chr','start','end','cytoband','qvalue','residualq','numgenes','genes')
  for (i in 1:nrow(gistictable))
  {
    tmp=as.character(gistictable[i,4])
    tmp1=unlist(strsplit(tmp,":"))
    tmp2=unlist(strsplit(tmp1[2],"-"))
    res[i,1]=tmp1[1]
    res[i,2]=tmp2[1]
    res[i,3]=tmp2[2]
    res[i,4]=as.character(gistictable[i,1])
    res[i,5]=as.character(gistictable[i,2])
    res[i,6]=as.character(gistictable[i,3])
    genes=as.character(gistictable[i,5])
    res[i,8]=genes
    genes=unlist(strsplit(genes,",",fixed = T))
    res[i,7]=length(genes)
  }
  #output=paste0(gisticdir,'/',prefix,'.gistic.cna.txt')
  return(res)
}

gisticfile="/fh/fast/dai_j/CancerGenomics/Chinese_EA/nature20805-s2Table 2.xlsx"
tcga_amp=readgisticcna(gisticfile,1)
tcga_del=readgisticcna(gisticfile,2)
tcga_scc_amp=readgisticcna(gisticfile,3)
tcga_scc_del=readgisticcna(gisticfile,4)
                           
gisticfile="/fh/fast/dai_j/CancerGenomics/Chinese_EA/ng.3659_ampgistic.txt"
uk_amp=readgisticcna1(gisticfile)
gisticfile="/fh/fast/dai_j/CancerGenomics/Chinese_EA/ng.3659_delgistic.txt"
uk_del=readgisticcna1(gisticfile)

library(GenomicRanges)
overlapgistic=function(gistic1,gistic2)
{
  GR_gistic1=GRanges(seqnames = as.character(gistic1$chr),ranges = IRanges(start=as.numeric(gistic1$start),end=as.numeric(gistic1$end)))
  GR_gistic2=GRanges(seqnames = as.character(gistic2$chr),ranges = IRanges(start=as.numeric(gistic2$start),end=as.numeric(gistic2$end)))
  res=data.frame(matrix(NA,nrow=0,ncol=5))
  colnames(res)=c("chr","start","end","cytoband","genes")
  for (i in 1:length(GR_gistic1))
  {
    olap=intersect(GR_gistic2,GR_gistic1[i])
    olap=reduce(olap)
    if (length(olap)>0)
    {
      for (j in 1:length(olap))
      {
        res=rbind.data.frame(res,data.frame(chr=seqnames(olap)[j],start=start(olap)[j],end=end(olap)[j],
                                            cytoband=as.character(gistic1$cytoband[i]),genes=as.character(gistic1$genes[i])))
        
      }
    }
  }
  GR_res=GRanges(seqnames = res$chr,ranges = IRanges(start=res$start,end=res$end),cytoband=res$cytoband)
  GR_gistic=reduce(GR_res) #merge segments
  res1=data.frame(matrix(NA,nrow=0,ncol=ncol(res)))
  colnames(res1)=colnames(res)
  for (i in 1:length(GR_gistic))
  {
    for (j in 1:length(GR_res))
    {
      olap=subsetByOverlaps(GR_res[j],GR_gistic[i])
      if (length(olap)>0)
      {
        tmp=data.frame(chr=seqnames(GR_gistic)[i],start=start(GR_gistic)[i],end=end(GR_gistic)[i],
                       cytoband=res$cytoband[j],genes=res$genes[j])
        res1=rbind.data.frame(res1,tmp)
        break
      }
    }
  }
  
  #sort
  chrs=paste0("chr",c(1:22,"X","Y"))
  res2=NULL
  for (chr in chrs)
  {
    tmp=res1[res1$chr==chr,]
    if (nrow(tmp)>0)
    {
      if (nrow(tmp)==1)
      {
        res2=rbind.data.frame(res2,tmp)
      }else
      {
        idx=order(tmp$start)
        res2=rbind.data.frame(res2,tmp[idx,])
      }
    }
  }
  return(res2)
}

#uniq regions in gistic1
uniqgistic1=function(gistic1,gistic2)
{
  comgistic=overlapgistic(gistic1,gistic2)
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(gistic1)))
  colnames(res)=colnames(gistic1)
  GR_gistic1=GRanges(seqnames = gistic1$chr,ranges = IRanges(start=as.numeric(gistic1$start),end=as.numeric(gistic1$end)))
  GR_comgistic=GRanges(seqnames = comgistic$chr,ranges = IRanges(start=comgistic$start,end=comgistic$end))
  for (i in 1:nrow(gistic1))
  {
    olap=subsetByOverlaps(GR_comgistic,GR_gistic1[i])
    if (length(olap)==0)
    {
      res=rbind.data.frame(res,gistic1[i,])
    }
  }
  #sort
  chrs=paste0("chr",c(1:22,"X","Y"))
  res1=NULL
  for (chr in chrs)
  {
    tmp=res[res$chr==chr,]
    if (nrow(tmp)>0)
    {
      if (nrow(tmp)==1)
      {
        res1=rbind.data.frame(res1,tmp)
      }else
      {
        idx=order(tmp$start)
        res1=rbind.data.frame(res1,tmp[idx,])
      }
    }
  }
  return(res1)
}

amp_tcga_ea_regions=uniqgistic1(tcga_amp,tcga_scc_amp)
#rows 27->17
amp_tcga_escc_regions=uniqgistic1(tcga_scc_amp,tcga_amp)
#rows 22->13
del_tcga_ea_regions=uniqgistic1(tcga_del,tcga_scc_del)
#rows 37->12
del_tcga_escc_regions=uniqgistic1(tcga_scc_del,tcga_del)
#rows 43->16

# amp_regions=overlapgistic(tcga_amp,uk_amp)
# del_regions=overlapgistic(tcga_del,uk_del)

extractregions2=function(regions=regions,wgstumors,freecdir,pcutoff)
{
  #use specific regions
  regions$chr=as.character(regions$chr)
  regions$chr=gsub("chr","",regions$chr)
  regions$start=as.numeric(regions$start)
  regions$end=as.numeric(regions$end)
  gr_regions=GRanges(seqnames=regions$chr,ranges=IRanges(start=regions$start,end=regions$end))
  cnres=widthres=cbind(regions,data.frame(matrix(NA,nrow=nrow(regions),ncol=length(wgstumors))))
  colnames(cnres)[(ncol(regions)+1):ncol(cnres)]=wgstumors
  colnames(widthres)=colnames(cnres)
  
  #sort
  chrs=c(1:22,"X","Y")
  res1=NULL
  for (chr in chrs)
  {
    tmp1=cnres[cnres$chr==chr,]
    if (nrow(tmp1)>0)
    {
      if (nrow(tmp1)==1)
      {
        res1=rbind.data.frame(res1,tmp1)
      }else
      {
        idx=order(tmp1$start)
        res1=rbind.data.frame(res1,tmp1[idx,])
      }
    }
  }
  cnres=widthres=res1
  #add rownames
  cnres$cytoband=as.character(cnres$cytoband)
  idxdup=duplicated(cnres$cytoband)
  rownames(cnres)[!idxdup]=cnres$cytoband[!idxdup]
  dupsign=0
  if (sum(idxdup)==1)
  {
    rownames(cnres)[idxdup]=paste0(cnres$cytoband[idxdup],"(1)")
  }
  if (sum(idxdup)>1)
  {
    idxdup=which(idxdup==T)
    rownames(cnres)[idxdup[1]]=paste0(cnres$cytoband[idxdup[1]],"(1)")
    dupsign=1
    for (i in 2:length(idxdup))
    {
      if(cnres$cytoband[idxdup[i]]==cnres$cytoband[idxdup[i-1]])
      {
        dupsign=dupsign+1
        rownames(cnres)[idxdup[i]]=paste0(cnres$cytoband[idxdup[i]],"(",dupsign,")")
      }else
      {
        dupsign=1
        rownames(cnres)[idxdup[i]]=paste0(cnres$cytoband[idxdup[i]],"(1)")
      }
    }
  }
  rownames(widthres)=rownames(cnres)
  for (i in 1:length(wgstumors))
  {
    ratiofile=paste0(freecdir,"/",wgstumors[i],"/ploid2degree3force0/",wgstumors[i],".pvalue.txt")
    ratiodata=read.table(file=ratiofile,header=T,sep="\t")
    
    ratiodata=ratiodata[ratiodata[,10]<=pcutoff & ratiodata$uncertainty<20, ]
    gr_ratiodata=GRanges(seqnames=ratiodata$chr,ranges=IRanges(start=ratiodata$start,end=ratiodata$end),copynumber=ratiodata$copy.number)
    for (j in 1:nrow(regions))
    {
      olap=subsetByOverlaps(gr_ratiodata,gr_regions[j,])
      if (length(olap)>0)
      {
        cns=mcols(olap)["copynumber"]
        cns=cns[,1]
        cnres[j,ncol(regions)+i]=median(cns,na.rm=T)
        widthres[j,ncol(regions)+i]=sum(width(olap)) 
      }
    }
  }
  return(result=list(cnres=cnres,widthres=widthres))
}

#more aggressive, call a region del if it has a del seg. 2:have both, 1:have amp,-1:have del
extractregions3=function(regions=regions,wgstumors,freecdir,pcutoff)
{
  #use specific regions
  regions$chr=as.character(regions$chr)
  regions$chr=gsub("chr","",regions$chr)
  regions$start=as.numeric(regions$start)
  regions$end=as.numeric(regions$end)
  gr_regions=GRanges(seqnames=regions$chr,ranges=IRanges(start=regions$start,end=regions$end))
  cnres=widthres=cbind(regions,data.frame(matrix(NA,nrow=nrow(regions),ncol=length(wgstumors))))
  colnames(cnres)[(ncol(regions)+1):ncol(cnres)]=wgstumors
  colnames(widthres)=colnames(cnres)
  
  #sort
  chrs=c(1:22,"X","Y")
  res1=NULL
  for (chr in chrs)
  {
    tmp1=cnres[cnres$chr==chr,]
    if (nrow(tmp1)>0)
    {
      if (nrow(tmp1)==1)
      {
        res1=rbind.data.frame(res1,tmp1)
      }else
      {
        idx=order(tmp1$start)
        res1=rbind.data.frame(res1,tmp1[idx,])
      }
    }
  }
  cnres=widthres=res1
  #add rownames
  cnres$cytoband=as.character(cnres$cytoband)
  idxdup=duplicated(cnres$cytoband)
  rownames(cnres)[!idxdup]=cnres$cytoband[!idxdup]
  dupsign=0
  if (sum(idxdup)==1)
  {
    rownames(cnres)[idxdup]=paste0(cnres$cytoband[idxdup],"(1)")
  }
  if (sum(idxdup)>1)
  {
    idxdup=which(idxdup==T)
    rownames(cnres)[idxdup[1]]=paste0(cnres$cytoband[idxdup[1]],"(1)")
    dupsign=1
    for (i in 2:length(idxdup))
    {
      if(cnres$cytoband[idxdup[i]]==cnres$cytoband[idxdup[i-1]])
      {
        dupsign=dupsign+1
        rownames(cnres)[idxdup[i]]=paste0(cnres$cytoband[idxdup[i]],"(",dupsign,")")
      }else
      {
        dupsign=1
        rownames(cnres)[idxdup[i]]=paste0(cnres$cytoband[idxdup[i]],"(1)")
      }
    }
  }
  rownames(widthres)=rownames(cnres)
  for (i in 1:length(wgstumors))
  {
    ratiofile=paste0(freecdir,"/",wgstumors[i],"/ploid2degree3force0/",wgstumors[i],".pvalue.txt")
    ratiodata=read.table(file=ratiofile,header=T,sep="\t")
    
    ratiodata=ratiodata[ratiodata[,10]<=pcutoff & ratiodata$uncertainty<20, ]
    gr_ratiodata=GRanges(seqnames=ratiodata$chr,ranges=IRanges(start=ratiodata$start,end=ratiodata$end),copynumber=ratiodata$copy.number)
    for (j in 1:nrow(regions))
    {
      olap=subsetByOverlaps(gr_ratiodata,gr_regions[j,])
      if (length(olap)>0)
      {
        cns=mcols(olap)["copynumber"]
        cns=cns[,1]
        if (any(cns>2) & any(cns<2))
        {
          cnres[j,ncol(regions)+i]=2 #have both amp and del
        }else
        {
          if (any(cns>2))
          {
            cnres[j,ncol(regions)+i]=1 #have only amp
          }
          if (any(cns<2))
          {
            cnres[j,ncol(regions)+i]=-1 #have only del
          }
        }
        widthres[j,ncol(regions)+i]=sum(width(olap)) 
      }
    }
  }
  return(result=list(cnres=cnres,widthres=widthres))
}

ampres1=extractregions2(regions=amp_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres2=extractregions2(regions=amp_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
library(ComplexHeatmap)
draw_4regions2=function(cnmat,opt="amp",opt1="ordercolumn")
{
  
  #cnmat=as.matrix(res$ratiores[,4:ncol(res$ratiores)])
  cnmat=as.matrix(cnmat)
  mat=matrix(0,nrow=nrow(cnmat),ncol=ncol(cnmat))
  rownames(mat)=rownames(cnmat)
  colnames(mat)=colnames(cnmat)
  for (i in 1:nrow(mat))
  {
    for (j in 1:ncol(mat))
    {
       if (!is.na(cnmat[i,j]) & cnmat[i,j]>2 & opt=="amp") mat[i,j]=1
       if (!is.na(cnmat[i,j]) & cnmat[i,j]<2 & opt!="amp") mat[i,j]=1
      
    }
  }
  if (opt=="amp")
  {
    if (opt1=="ordercolumn")
    {
      print(Heatmap(mat,column_dend_height = unit(30, "mm"),cluster_rows = T,name="gain event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="gain event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
    
          
  }else
  {
    if (opt1=="ordercolumn")
    {
      mat1=mat[,sample(1:ncol(mat))]
      print(Heatmap(mat1,column_dend_height = unit(30, "mm"),cluster_rows = T,name="loss event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="loss event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
    
  }
  
  return(mat)
}
ampcnmat=cbind.data.frame(ampres1$cnres[,6:ncol(ampres1$cnres)],ampres2$cnres[,6:ncol(ampres2$cnres)])
colnames(ampcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(ampcnmat)[(length(wgstumors1)+1):ncol(ampcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawamp=draw_4regions2(ampcnmat,opt="amp",opt1="nordercolumn")

delres1=extractregions2(regions=del_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions2(regions=del_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres1$cnres[,6:ncol(delres1$cnres)],delres2$cnres[,6:ncol(delres2$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
numdel=apply(delcnmat,1,function(x){
  sum(x<2,na.rm=T)
})

idx3regions=which(del_regions$genes %in% c("FHIT","WWOX","CDKN2A,CDKN2B,C9orf53,CDKN2B-AS1"))
drawdel=draw_4regions2(delcnmat[idx3regions,],opt="del",opt1="ordercolumn")
drawdel=draw_4regions2(delcnmat[numdel>4,],opt="del",opt1="nordercolumn")
drawdel=draw_4regions2(delcnmat[numdel>4,],opt="del",opt1="ordercolumn")
idxrows=which(rownames(delcnmat) %in% c("4q22.1","7q31.1","20p12.1","16q23.1","16p13.3","8p12","22q11.1"))
drawdel=draw_4regions2(delcnmat[idxrows,],opt="del",opt1="ordercolumn")

drawdel=draw_4regions2(delcnmat,opt="del",opt1="ordercolumn")
drawdel=draw_4regions2(delcnmat,opt="del",opt1="nordercolumn")

##

library(ComplexHeatmap)
#work with extractregions3
draw_4regions3=function(cnmat,opt="amp",opt1="ordercolumn",numamp=0)
{
  
  #cnmat=as.matrix(res$ratiores[,4:ncol(res$ratiores)])
  cnmat=as.matrix(cnmat)
  mat=matrix(0,nrow=nrow(cnmat),ncol=ncol(cnmat))
  rownames(mat)=rownames(cnmat)
  colnames(mat)=colnames(cnmat)
  if (opt != "display") #need to process
  {
    for (i in 1:nrow(mat))
    {
      for (j in 1:ncol(mat))
      {
        if (!is.na(cnmat[i,j]) & cnmat[i,j] != -1 & opt=="amp") mat[i,j]=1
        if (!is.na(cnmat[i,j]) & cnmat[i,j] !=1 & opt=="del") mat[i,j]=1
        if (opt=="combined") #first numamp rows of amp then del
        {
          if (i<=numamp)
          {
            rownames(mat)[i]=paste0("A",rownames(cnmat)[i])
            if (!is.na(cnmat[i,j]) & cnmat[i,j] != -1 ) mat[i,j]=1
          }else
          {
            rownames(mat)[i]=paste0("D",rownames(cnmat)[i])
            if (!is.na(cnmat[i,j]) & cnmat[i,j] !=1 ) mat[i,j]=1
          }
        }
      }
    }
  }else
  {
    mat=cnmat
  }
  if (opt=="amp")
  {
    if (opt1=="ordercolumn")
    {
      #mat1=mat[,sample(1:ncol(mat))]
      mat1=cbind.data.frame(mat[,sample(1:length(wgstumors1))],mat[,sample((length(wgstumors1)+1):ncol(mat))])
      print(Heatmap(mat1,column_dend_height = unit(30, "mm"),cluster_rows = T,name="gain event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="gain event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
  }
  if (opt=="del")
  {
    if (opt1=="ordercolumn")
    {
      #mat1=mat[,sample(1:ncol(mat))]
      mat1=cbind.data.frame(mat[,sample(1:length(wgstumors1))],mat[,sample((length(wgstumors1)+1):ncol(mat))])
      print(Heatmap(mat1,column_dend_height = unit(30, "mm"),cluster_rows = T,name="loss event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="loss event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
  }
  if (opt=="combined")
  {
    if (opt1=="ordercolumn")
    {
      #mat1=mat[,sample(1:ncol(mat))]
      mat1=cbind.data.frame(mat[,sample(1:length(wgstumors1))],mat[,sample((length(wgstumors1)+1):ncol(mat))])
      print(Heatmap(mat1,column_dend_height = unit(30, "mm"),cluster_rows = T,name="event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
  }
  if (opt=="display")
  {
    if (opt1=="ordercolumn")
    {
      #mat1=mat[,sample(1:ncol(mat))]
      mat1=cbind.data.frame(mat[,sample(1:length(wgstumors1))],mat[,sample((length(wgstumors1)+1):ncol(mat))])
      print(Heatmap(mat1,column_dend_height = unit(30, "mm"),cluster_rows = T,name="event",show_column_names=T,show_column_dend=T,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }else
    {
      print(Heatmap(mat,cluster_rows = T,name="event",show_column_names=T,cluster_columns = F,
                    rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                    heatmap_legend_param = list(at = c(0,1), labels = c("no", "yes"))))
    }
  }
  
  return(mat)
}

#only use tcga_ea
ampres1=extractregions3(regions=amp_tcga_ea_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres2=extractregions3(regions=amp_tcga_ea_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
ampcnmat=cbind.data.frame(ampres1$cnres[,9:ncol(ampres1$cnres)],ampres2$cnres[,9:ncol(ampres2$cnres)])
colnames(ampcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(ampcnmat)[(length(wgstumors1)+1):ncol(ampcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawamp=draw_4regions3(ampcnmat,opt="amp",opt1="ordercolumn")

#overlap ea with escc
amp_regions=overlapgistic(tcga_amp,tcga_scc_amp) #10 regions
ampres11=extractregions3(regions=amp_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres22=extractregions3(regions=amp_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
ampcnmat=cbind.data.frame(ampres11$cnres[,6:ncol(ampres11$cnres)],ampres22$cnres[,6:ncol(ampres22$cnres)])
colnames(ampcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(ampcnmat)[(length(wgstumors1)+1):ncol(ampcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawamp=draw_4regions3(ampcnmat,opt="amp",opt1="ordercolumn")
num=apply(drawamp,1,function(x){
  sum(x==1,na.rm=T)
})
draw_4regions3(drawamp[num>0,],opt="display",opt1="ordercolumn")

#overlap with uk_amp
amp_regions1=overlapgistic(amp_regions,uk_amp) # 6 regions
ampres11=extractregions3(regions=amp_regions1,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres22=extractregions3(regions=amp_regions1,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
ampcnmat=cbind.data.frame(ampres11$cnres[,6:ncol(ampres11$cnres)],ampres22$cnres[,6:ncol(ampres22$cnres)])
colnames(ampcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(ampcnmat)[(length(wgstumors1)+1):ncol(ampcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawamp=draw_4regions3(ampcnmat,opt="amp",opt1="ordercolumn")


#only use tcga_ea
delres1=extractregions3(regions=del_tcga_ea_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=del_tcga_ea_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres1$cnres[,9:ncol(delres1$cnres)],delres2$cnres[,9:ncol(delres2$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawdel=draw_4regions3(delcnmat,opt="del",opt1="ordercolumn")
idxrow=which(! rownames(delcnmat) %in% c("9p23","4q34.3"))
drawdel=draw_4regions3(delcnmat[idxrow,],opt="del",opt1="ordercolumn")

#overlap ea with escc
del_regions=overlapgistic(tcga_del,tcga_scc_del)
delres1=extractregions3(regions=del_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=del_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres1$cnres[,6:ncol(delres1$cnres)],delres2$cnres[,6:ncol(delres2$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawdel=draw_4regions3(delcnmat,opt="del",opt1="ordercolumn")
num=apply(drawdel,1,function(x){
  sum(x==1,na.rm=T)
})
draw_4regions3(drawdel,opt="display",opt1="ordercolumn")
draw_4regions3(drawdel[num>10, ],opt="display",opt1="ordercolumn")

#overlap with uk_del
del_regions1=overlapgistic(del_regions,uk_del)
delres11=extractregions3(regions=del_regions1,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres22=extractregions3(regions=del_regions1,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres11$cnres[,6:ncol(delres11$cnres)],delres22$cnres[,6:ncol(delres22$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawdel=draw_4regions3(delcnmat,opt="del",opt1="ordercolumn")
num=apply(drawdel,1,function(x){
  sum(x==1,na.rm=T)
})
draw_4regions3(drawdel,opt="display",opt1="ordercolumn")
draw_4regions3(drawdel[num>10, ],opt="display",opt1="ordercolumn")


##use overlap ea and scc,combine amp and del
ampres1=extractregions3(regions=amp_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres2=extractregions3(regions=amp_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat1=cbind.data.frame(ampres1$cnres[,6:ncol(ampres1$cnres)],ampres2$cnres[,6:ncol(ampres2$cnres)])
delres1=extractregions3(regions=del_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=del_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat2=cbind.data.frame(delres1$cnres[,6:ncol(delres1$cnres)],delres2$cnres[,6:ncol(delres2$cnres)])
cnmat=rbind(cnmat1,cnmat2)
colnames(cnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(cnmat)[(length(wgstumors1)+1):ncol(cnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
num=apply(drawres,1,function(x){
  sum(x!=0,na.rm=T)
})
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
drawres1=draw_4regions3(drawres[num>9,],opt="display",opt1="ordercolumn")


##only use tcga_ea,combine amp and del
ampres1=extractregions3(regions=amp_tcga_ea_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres2=extractregions3(regions=amp_tcga_ea_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat1=cbind.data.frame(ampres1$cnres[,9:ncol(ampres1$cnres)],ampres2$cnres[,9:ncol(ampres2$cnres)])
delres1=extractregions3(regions=del_tcga_ea_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=del_tcga_ea_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat2=cbind.data.frame(delres1$cnres[,9:ncol(delres1$cnres)],delres2$cnres[,9:ncol(delres2$cnres)])
cnmat=rbind(cnmat1,cnmat2)
colnames(cnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(cnmat)[(length(wgstumors1)+1):ncol(cnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
num=apply(drawres,1,function(x){
  sum(x!=0,na.rm=T)
})
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
drawres1=draw_4regions3(drawres[num>1,],opt="del",opt1="ordercolumn")


##use overlap ea and scc and uk,combine amp and del
ampres1=extractregions3(regions=amp_regions1,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres2=extractregions3(regions=amp_regions1,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat1=cbind.data.frame(ampres1$cnres[,6:ncol(ampres1$cnres)],ampres2$cnres[,6:ncol(ampres2$cnres)])
delres1=extractregions3(regions=del_regions1,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=del_regions1,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
cnmat2=cbind.data.frame(delres1$cnres[,6:ncol(delres1$cnres)],delres2$cnres[,6:ncol(delres2$cnres)])
cnmat=rbind(cnmat1,cnmat2)
colnames(cnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(cnmat)[(length(wgstumors1)+1):ncol(cnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
num=apply(drawres,1,function(x){
  sum(x!=0,na.rm=T)
})
drawres=draw_4regions3(cnmat,opt="combined",opt1="ordercolumn",numamp = nrow(cnmat1))
drawres1=draw_4regions3(drawres[num>9,],opt="display",opt1="ordercolumn")
draw_4regions3(drawres1[c(1,4,5,6,7),],opt="display",opt1="ordercolumn")


#combine ea and escc
delres1=extractregions3(regions=rbind(del_tcga_ea_regions,del_tcga_escc_regions),wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions3(regions=rbind(del_tcga_ea_regions,del_tcga_escc_regions),wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres1$cnres[,9:ncol(delres1$cnres)],delres2$cnres[,9:ncol(delres2$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
drawdel=draw_4regions3(delcnmat,opt="del",opt1="ordercolumn")





delres1=extractregions2(regions=del_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
delres2=extractregions2(regions=del_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)
delcnmat=cbind.data.frame(delres1$cnres[,6:ncol(delres1$cnres)],delres2$cnres[,6:ncol(delres2$cnres)])
colnames(delcnmat)[1:length(wgstumors1)]=paste0("US-EA",1:length(wgstumors1))
colnames(delcnmat)[(length(wgstumors1)+1):ncol(delcnmat)]=paste0("CH-EA",1:length(wgstumors2))
numdel=apply(delcnmat,1,function(x){
  sum(x<2,na.rm=T)
})

idx3regions=which(del_regions$genes %in% c("FHIT","WWOX","CDKN2A,CDKN2B,C9orf53,CDKN2B-AS1"))
drawdel=draw_4regions2(delcnmat[idx3regions,],opt="del",opt1="ordercolumn")


#use regions generated by china EA
gisticfile="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1/del_genes.conf_95.txt"
henan_del=readgisticcna(gisticfile)
gisticfile="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1/amp_genes.conf_95.txt"
henan_amp=readgisticcna(gisticfile)
henan_del_regions=henan_del[,c(1:4,8)]
henan_amp_regions=henan_amp[,c(1:4,8)]
ampres11=extractregions2(regions=henan_amp_regions,wgstumors=wgstumors1,freecdir=freecdir1,pcutoff = cutoff)
ampres22=extractregions2(regions=henan_amp_regions,wgstumors=wgstumors2,freecdir=freecdir2,pcutoff = cutoff)