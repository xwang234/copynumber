#!/usr/bin/env Rscript

library(GenomicRanges)
computep=function(x,y)
{
  print(mean(x))
  print(mean(y))
  res=2*(1-pnorm(abs((mean(x)-mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))))
}
chrs=c(paste0("chr",c(1:22)),"chrX","chrY")
chrds=c(1:22,"X","Y")
cnvlength=function(tumors,freecdir,subfolder,pcutoff=NULL)
{
  res=data.frame(matrix(0,nrow=length(tumors),ncol=3*length(chrds)+1))
  colnames(res)=c("sample",paste0("gain",chrds),paste0("loss",chrds),paste0("loh",chrds))
  for (i in 1:length(tumors))
  {
    tumor=tumors[i]
    res[i,1]=tumor
    cnvfile=paste0(freecdir,'/',tumor,'/',subfolder,'/',tumor,'.pvalue.txt')
    cnvtable=read.table(cnvfile,header=T,sep="\t")
    cnvtable[,1]=as.character(cnvtable[,1])
    pval1col=10 #pvalue column
    pval2col=11 #pvalue column
    statcol=5 #status column
    somaticcol=8 #somatic or germline
    cnvtable[,somaticcol]=as.character(cnvtable[,somaticcol])
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="gain" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="gain" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      
      res[i,k+1]=0
      if (sum(idx1)>0)
      {
        tmptable1=tmptable[idx1,]
        gr_tmptable1=GRanges(seqnames=tmptable1$chr,ranges=IRanges(start=tmptable1$start,end=tmptable1$end))
        gr_tmptable1=reduce(gr_tmptable1) #combine the possible overlaps
        res[i,k+1]=sum(as.numeric(width(gr_tmptable1)))
      }
    }
    
    #idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      res[i,k+length(chrds)+1]=0
      if (sum(idx1)>0)
      {
        tmptable1=tmptable[idx1,]
        gr_tmptable1=GRanges(seqnames=tmptable1$chr,ranges=IRanges(start=tmptable1$start,end=tmptable1$end))
        gr_tmptable1=reduce(gr_tmptable1) #combine the possible overlaps
        res[i,k+length(chrds)+1]=sum(as.numeric(width(gr_tmptable1)))
      }
    }
    
    #idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      res[i,k+2*length(chrds)+1]=0
      if (sum(idx1)>0)
      {
        tmptable1=tmptable[idx1,]
        gr_tmptable1=GRanges(seqnames=tmptable1$chr,ranges=IRanges(start=tmptable1$start,end=tmptable1$end))
        gr_tmptable1=reduce(gr_tmptable1) #combine the possible overlaps
        res[i,k+2*length(chrds)+1]=sum(as.numeric(width(gr_tmptable1)))
      }
    }
  }
  return(res)
}

#count cnv numbers
cnvnum=function(tumors,freecdir,subfolder,pcutoff=NULL)
{
  res=data.frame(matrix(0,nrow=length(tumors),ncol=3*length(chrds)+1))
  colnames(res)=c("sample",paste0("gain",chrds),paste0("loss",chrds),paste0("loh",chrds))
  for (i in 1:length(tumors))
  {
    tumor=tumors[i]
    res[i,1]=tumor
    cnvfile=paste0(freecdir,'/',tumor,'/',subfolder,'/',tumor,'.pvalue.txt')
    cnvtable=read.table(cnvfile,header=T,sep="\t")
    cnvtable[,1]=as.character(cnvtable[,1])
    pval1col=10 #pvalue column
    pval2col=11 #pvalue column
    statcol=5 #status column
    somaticcol=8 #somatic or germline
    cnvtable[,somaticcol]=as.character(cnvtable[,somaticcol])
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="gain" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="gain" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      res[i,k+1]=0
      if (sum(idx1)>0)
      {
        res[i,k+1]=sum(idx1)
      }
    }
    
    #idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="loss" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      res[i,k+length(chrds)+1]=0
      if (sum(idx1)>0)
      {
        res[i,k+length(chrds)+1]=sum(idx1)
      }
    }
    
    #idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    if (is.null(pcutoff))
    {
      idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic"
    }else
    {
      idx=cnvtable[,statcol]=="normal" & cnvtable[,pval1col]<=pcutoff & cnvtable[,somaticcol]=="somatic"
    }
    tmptable=cnvtable[idx,]
    for (k in 1:length(chrds))
    {
      chr=chrds[k]
      idx1=tmptable[,1]==chr
      res[i,k+2*length(chrds)+1]=0
      if (sum(idx1)>0)
      {
        res[i,k+2*length(chrds)+1]=sum(idx1)
      }
    }
  }
  return(res)
}

chrlen=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,
         107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566
)
plotcnv=function(normals,tumors,freecdir,subfolder,selchrids,cnvtype,outfig) #plot loss this time
{
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  selchrlen=chrlen[selchrids]
  selchrstart=rep(0,length(selchrlen))
  for (i in 1:(length(selchrlen)-1))
  {
    
    selchrstart[i+1]=selchrstart[i]+selchrlen[i]
  }
  if (cnvtype=="gain")
  {
    mycolor="red"
  }else
  {
    mycolor="black"
  }
  xmin=0
  xmax=selchrstart[length(selchrids)]+selchrlen[length(selchrids)]+10000
  ymin=0
  ymax=length(normals)+1
  
  png(outfig, width = 12, height = 8, units = 'in', res=300)
  tmp=par("mar")
  par(mar=tmp+c(0,1.5,0,0))
  plot(c(xmin,xmax),c(ymin,ymax),type="n",xlab='chromosome',ylab='sample',xaxt = 'n',cex.lab=1.5,cex.axis=2)
  xlabpos=selchrstart[1:length(selchrids)]+selchrlen[1:length(selchrids)]/2
  xlabels=chrs[selchrids]
  axis(1,at=xlabpos,labels = xlabels,cex.axis=1.5)
  #plot(c(xmin,xmax),c(ymin,ymax),type="n",xlab='chromosome',ylab='sample',cex.lab=1.5,cex.axis=2)
  for (i in 1:length(normals))
  {
    normal=normals[i]
    tumor=tumors[i]
    res[i,1]=tumor
    cnvfile=paste0(freecdir,'/',normal,'/',subfolder,'/',tumor,'.pvalue.txt')
    cnvtable=read.table(cnvfile,header=T,sep="\t")
    cnvtable[,1]=as.character(cnvtable[,1])
    pvalcol=10 #pvalue column
    statcol=5 #status column
    idx=cnvtable[,statcol]==cnvtype & cnvtable[,pvalcol]<=0.05 #for loss
    tmptable=cnvtable[idx,]
    for (j in 1:length(selchrids))
    {
      selchr=selchrids[j]
      idx1=tmptable[,1]==selchr
      if (sum(idx1)>0)
      {
        tmptable1=tmptable[idx1,]
        for (k in 1:nrow(tmptable1))
        {
          xs=c(tmptable1[k,2],tmptable1[k,3])+selchrstart[j]
          ys=c(i,i)
          lines(xs,ys,col=mycolor,lwd=8)
        }
      }
    }
  }
  dev.off()
}

readboundaries=function(rawtable,opt=0)
{
  #split boundaries into chr start end
  rawtable$boundaries=as.character(rawtable$boundaries)
  chr=data.frame(chr=matrix(NA,nrow=nrow(rawtable),ncol=1))
  start=data.frame(start=matrix(NA,nrow=nrow(rawtable),ncol=1))
  end=data.frame(end=matrix(NA,nrow=nrow(rawtable),ncol=1))
  for (i in 1:nrow(rawtable))
  {
    tmp=as.character(rawtable$boundaries[i])
    tmp1=unlist(strsplit(tmp,":"))
    tmp2=unlist(strsplit(tmp1[2],"-"))
    chr[i,1]=tmp1[1]
    if (opt==0)
    {
      start[i,1]=as.integer(tmp2[1])
      end[i,1]=as.integer(tmp2[2])
    }else
    {
      start[i,1]=as.numeric(tmp2[1])*10^6
      end[i,1]=as.numeric(tmp2[2])*10^6
    }
    
  }
  res=cbind(rawtable,chr,start,end)
  return(res)
}

#count Y in the all results table to generate consensus table
countY=function(cnvtable,ycol)
{
  for (i in ycol)
  {
    cnvtable[,i]=as.character(cnvtable[,i])
  }
  
  res=rep(0,nrow(cnvtable))
  for (i in 1:nrow(cnvtable))
  {
    idx=is.na(cnvtable[i,])
    if(sum(idx)>0)
    cnvtable[i,idx]=rep("N",sum(idx))
    for (j in ycol)
    {
      res[i]=res[i]+as.integer(cnvtable[i,j]=="Y")
    }
  }
  return(res)
}

addstatus=function(reftable,cnvtype)
{
  status=data.frame(max(ncol=1,nrow=nrow(reftable)))
  colnames(status)="status"
  status=rep(cnvtype,nrow(reftable))
  res=cbind(reftable,status)
}

checkoverlap=function(cnvtable,reftable1)
{
  cnvtable[,"status"]=as.character(cnvtable[,"status"])
  reftable1[,"status"]=as.character(reftable1[,"status"])
  gr_cnvtable=GRanges(seqnames=cnvtable$chr,ranges=IRanges(start=cnvtable$start,end=cnvtable$end),status=cnvtable$status)
  gr_reftable1=GRanges(seqnames=reftable1$chr,ranges=IRanges(start=reftable1$start,end=reftable1$end),status=reftable1$status)
  options(warn=-1) #in case no overlap
  olap=subsetByOverlaps(gr_cnvtable,gr_reftable1)
  options(warn=0)
  #res=ifelse(length(olap)>0,1,0) #1 if overlap
  #assume only one overlap
  if (length(olap)>0)
  {
    res=max(as.numeric(width(olap)))/as.numeric(width(gr_reftable1))
    res=format(res,digits = 2)
    if (res>=1) res=1
  }else
  {
    res=0
  }
  
  return(res)
}

#count overlaps of each cnv in cnvtable with cnv results in freecdir/subfolder,
overlapcnvcount=function(reftable,normals,tumors,freecdir,subfolder)
{
  reftable[,"chr"]=as.character(reftable[,"chr"])
  res=data.frame(matrix(NA,nrow=nrow(reftable),ncol=4+length(normals)))
  res[,1:4]=reftable[,c("chr","start","end","status")]
  colnames(res)=c("chr","start","end","status",tumors)
  res[,1]=as.character(res[,1])
  res[,4]=as.character(res[,4])
  
  #column first
  for (i in 1:length(normals))
  {
    normal=normals[i]
    tumor=tumors[i]
    cnvfile=paste0(freecdir,'/',normal,'/',subfolder,'/',tumor,'.pvalue.txt')
    cnvtable=read.table(cnvfile,header=T,sep="\t")
    cnvtable[,1]=as.character(cnvtable[,1])
    if (!grepl(cnvtable[1,1],"chr"))
    {
      cnvtable[,1]=paste0("chr",cnvtable[,1])
    }
    pval1col=10 #pvalue column
    pval2col=11
    statcol=5 #status column
    somaticcol=8 #somatic or germline
    cnvtable[,somaticcol]=as.character(cnvtable[,somaticcol])
    idx=cnvtable[,pval1col]<=0.05 & cnvtable[,pval2col]<=0.05 & cnvtable[,somaticcol]=="somatic" #filter
    cnvtable=cnvtable[idx,]
    
    #row
    for (j in 1:nrow(reftable))
    {
      res[j,i+4]=checkoverlap(cnvtable,reftable[j,])
      
    }
  }
  return(res)
}

#read xiaohong's data  
readtable=function(file)
{
  tmptable=read.table(file,header=F,sep=" ")
  res=data.frame(matrix(NA,ncol=6,nrow=nrow(tmptable)))
  colnames(res)=c("chr","start","end","status","freq1","freq2")
  tmptable[,1]=as.character(tmptable[,1])
  tmptable[,2]=as.character(tmptable[,2])
  tmptable[,3]=as.character(tmptable[,3])
  res[,"chr"]=paste0("chr",gsub(":","",tmptable[,1]))
  res[,c("freq1","freq2")]=tmptable[,c(5,6)]
  for (i in 1:nrow(tmptable))
  {
    tmp=unlist(strsplit(tmptable[i,2],"–"))
    res[i,"start"]=as.integer(tmp[[1]])*10^6
    res[i,"end"]=as.integer(tmp[[2]])*10^6
    if (tmptable[i,3]=="Loss")
    {
      res[i,"status"]="loss"
    }else if (tmptable[i,3]=="Gain")
    {
      res[i,"status"]="gain"
    }else
    {
      res[i,"status"]="normal"
    }
  }
  return(res)
}

#henandata:
tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
freecdir="/fh/scratch/delete30/dai_j/henan/freec"
subfolder="ploid2degree3force0"
nullfile="/fh/scratch/delete30/dai_j/henan/freec/42A/ploid2degree3force0/42A.pvalue.txt"
nulldata=read.table(nullfile,header=T,sep="\t")
nulldata=nulldata[order(nulldata[,10]),]
#cutoff=nulldata[as.integer(0.01*nrow(nulldata)),10]
cutoff=nulldata[as.integer(0.05*nrow(nulldata)),10]

henanres=cnvlength(tumors,freecdir,subfolder)
write.table(henanres,file="/fh/scratch/delete30/dai_j/henan/freec/henancnvlengthcount.txt",row.names = F,sep="\t",quote=F)
henannumres=cnvnum(tumors,freecdir,subfolder)
write.table(henannumres,file="/fh/scratch/delete30/dai_j/henan/freec/henancnvnumcount.txt",row.names = F,sep="\t",quote=F)
#using pvalue cutoff from normal data
henanres1=cnvlength(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(henanres1,file="/fh/scratch/delete30/dai_j/henan/freec/henancnvlengthcount_stringent.txt",row.names = F,sep="\t",quote=F)
henannumres1=cnvnum(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(henannumres1,file="/fh/scratch/delete30/dai_j/henan/freec/henancnvnumcount_stringent.txt",row.names = F,sep="\t",quote=F)


#henannormaldata:
tumors=paste0(c(3,11,13,15,17,25,29,33,37,41)+1,"A")
freecdir="/fh/scratch/delete30/dai_j/henan/freec"
subfolder="ploid2degree3force0"
henannormalsres=cnvlength(tumors,freecdir,subfolder)
write.table(henannormalsres,file="/fh/scratch/delete30/dai_j/henan/freec/henannormalscnvlengthcount.txt",row.names = F,sep="\t",quote=F)
henannormalsnumres=cnvnum(tumors,freecdir,subfolder)
write.table(henannormalsnumres,file="/fh/scratch/delete30/dai_j/henan/freec/henannormalscnvnumcount.txt",row.names = F,sep="\t",quote=F)


#dulakdata
tumors=c('SRR1001842','SRR1002713','SRR999423','SRR1001466','SRR1002670','SRR1001823','SRR999489','SRR1002343','SRR1002722','SRR1002656','SRR1002929','SRR999438','SRR1001915','SRR999594','SRR1001868','SRR1001635')
freecdir="/fh/scratch/delete30/dai_j/freec"
subfolder="ploid2degree3force0"
dulakres=cnvlength(tumors,freecdir,subfolder)
write.table(dulakres,file="/fh/scratch/delete30/dai_j/freec/dulakcnvlengthcount.txt",row.names = F,sep="\t",quote=F)
dulaknumres=cnvnum(tumors,freecdir,subfolder)
write.table(dulaknumres,file="/fh/scratch/delete30/dai_j/freec/dulakcnvnumcount.txt",row.names = F,sep="\t",quote=F)


#dulaknormaldata cutoff
dulakres1=cnvlength(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(dulakres1,file="/fh/scratch/delete30/dai_j/freec/dulakcnvlengthcount_stringent.txt",row.names = F,sep="\t",quote=F)
dulaknumres1=cnvnum(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(dulaknumres1,file="/fh/scratch/delete30/dai_j/freec/dulakcnvnumcount_stringent.txt",row.names = F,sep="\t",quote=F)

#esccdata
tumors=paste0("T",c(1:6,8:18))
freecdir="/fh/scratch/delete30/dai_j/escc/freec"
subfolder="ploid2degree3force0"
esccres=cnvlength(tumors,freecdir,subfolder)
write.table(esccres,file="/fh/scratch/delete30/dai_j/escc/freec/escccnvlengthcount.txt",row.names = F,sep="\t",quote=F)
esccnumres=cnvnum(tumors,freecdir,subfolder)
write.table(esccnumres,file="/fh/scratch/delete30/dai_j/escc/freec/escccnvnumcount.txt",row.names = F,sep="\t",quote=F)

#escc normaldata cutoff
esccres1=cnvlength(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(esccres1,file="/fh/scratch/delete30/dai_j/escc/freec/escccnvlengthcount_stringent.txt",row.names = F,sep="\t",quote=F)
esccnumres1=cnvnum(tumors,freecdir,subfolder,pcutoff=cutoff)
write.table(esccnumres1,file="/fh/scratch/delete30/dai_j/escc/freec/escccnvnumcount_stringent.txt",row.names = F,sep="\t",quote=F)


#for gains:
pvalues=rep(NA,72)
for (i in 1:72)
{
  print(i)
  pvalues[i]=computep(dulakres[,i+1],henanres[,i+1])
}
idx=which(pvalues<0.05)
print(idx-24)
selcol=paste0("loss",c(5,16,17,18,20))
for (i in 1:length(selcol))
{
  print(selcol[i])
  print(computep(dulakres[,selcol[i]],henanres[,selcol[i]]))
}

#for whole length, including all chromosomes:
henangain=rowSums(henanres[,2:25])
dulakgain=rowSums(dulakres[,2:25])
print(computep(dulakgain,henangain))
henanloss=rowSums(henanres[,26:49])
dulakloss=rowSums(dulakres[,26:49])
print(computep(dulakloss,henanloss))
henanloh=rowSums(henanres[,50:73])
dulakloh=rowSums(dulakres[,50:73])
print(computep(dulakloh,henanloh))


#plotcnvs
normals=c('2A','4A','6A','8A','10A','12A')
tumors=c('1A','3A','5A','7A','9A','11A')
freecdir="/fh/scratch/delete30/dai_j/henan/freec"
subfolder="w2000"
cnvtype="loss"
selchrids=c(5,6,17,18,20)
tmp=NULL
for (i in 1:length(selchrids))
{
  tmp=paste(tmp,selchrids[i],sep="_")
}

outfig=paste0("CNV_plot_chr",tmp,".png")
plotcnv(normals,tumors,freecdir,subfolder,selchrids,cnvtype,outfig)


freecdir="/fh/scratch/delete30/dai_j/freec/out1"
subfolder="ploid2degree3force0"
outfig=paste0("/fh/scratch/delete30/dai_j/freec/DulakCNV_plot_chr",tmp,".png")
plotcnv(wgsnormals,wgstumors,freecdir,subfolder,selchrids,cnvtype,outfig)


#count overlap cnvs in each sample
normals=c('2A','4A','6A','8A','10A','12A')
tumors=c('1A','3A','5A','7A','9A','11A')
freecdir="/fh/scratch/delete30/dai_j/henan/freec"
subfolder="w2000"
#subfolder="ploids"

#first try the consensus table
#for loss
consensus_del=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/owes_rx0_conf0.95_armpeel0_brlen0.98_broad1/result_allcytoband_del.txt",header=T,sep="\t")
#consensus_amp=readboundaries(consensus_amp)
ycol=4:9
idx=rowSums(consensus_del[,ycol])>=3
consensus_del=consensus_del[idx,]
#A reftable needs to have columns chr,start,end,status
consensus_del=addstatus(consensus_del,"loss")
henan_consensus_del_count=overlapcnvcount(consensus_del,normals,tumors,freecdir,subfolder)
write.table(henan_consensus_del_count,file="henan_consensus_del_count.txt",col.names = T,sep="\t",quote=F,row.names=F)
#for gain
consensus_amp=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/owes_rx0_conf0.95_armpeel0_brlen0.98_broad1/result_allcytoband_amp.txt",header=T,sep="\t")
#consensus_amp=readboundaries(consensus_amp)
ycol=4:9
idx=rowSums(consensus_amp[,ycol])>=3
consensus_amp=consensus_amp[idx,]
#A reftable needs to have columns chr,start,end,status
consensus_amp=addstatus(consensus_amp,"gain")
henan_consensus_amp_count=overlapcnvcount(consensus_amp,normals,tumors,freecdir,subfolder)
write.table(henan_consensus_amp_count,file="henan_consensus_amp_count.txt",col.names = T,sep="\t",quote=F,row.names=F)

#second check with Dulak's regions
dulak_del=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/owes_rx0_conf0.95_armpeel0_brlen0.98_broad1/result_Dulak_our_del.txt",header=T,sep="\t")
dulak_del=readboundaries(dulak_del)
dulak_del=addstatus(dulak_del,"loss")
henan_dulak_del_count=overlapcnvcount(dulak_del,normals,tumors,freecdir,subfolder)
write.table(henan_dulak_del_count,file="henan_dulak_del_count.txt",col.names = T,sep="\t",quote=F,row.names=F)
dulak_amp=read.table("/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/owes_rx0_conf0.95_armpeel0_brlen0.98_broad1/result_Dulak_our_amp.txt",header=T,sep="\t")
dulak_amp=readboundaries(dulak_amp,opt=1) #times 10^6 for start and end
dulak_amp=addstatus(dulak_amp,"gain")
henan_dulak_amp_count=overlapcnvcount(dulak_amp,normals,tumors,freecdir,subfolder)
write.table(henan_dulak_amp_count,file="henan_dulak_amp_count.txt",col.names = T,sep="\t",quote=F,row.names=F)

#third check with Xiaohong's result Cancer prevention 2015
file="xiaohong.txt" #remember to remove spaces before (
xiaohong_ref=readtable(file)
henan_xiaohong_count=overlapcnvcount(xiaohong_ref,normals,tumors,freecdir,subfolder)
write.table(henan_xiaohong_count,file="henan_xiaohong_count.txt",col.names=T,row.names=F,sep="\t",quote=F)
