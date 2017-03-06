#! /usr/bin/env Rscript
#SBATCH -t 0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

gettumor=function(folder)
{
  filelist=list.files(folder)
  tumorfile=grep("*.pileup.gz_CNVs",filelist,value=T)
  tumor=strsplit(tumorfile,'.',fixed=T)
  tumor=tumor[[1]][1]
}

args <- commandArgs(trailingOnly = TRUE)
library(rtracklayer)
library(GenomicRanges)
preface=as.character(args[1]) #freecdir/wgstumor
print(preface)

ratiofile=paste0(preface,".pileup.gz_ratio.txt")
cnvfile=paste0(preface,".pileup.gz_CNVs")
output=paste0(preface,".pvalue.txt")

dataTable <-read.table(ratiofile, header=TRUE)
ratio<-data.frame(dataTable)

#dataTable <- read.table(args[4], header=FALSE)
dataTable <- read.table(cnvfile, header=FALSE, fill=T)
cnvs<- data.frame(dataTable)

ratio$Ratio[which(ratio$Ratio==-1)]=NA

cnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  
ratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)

overlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)
normals <- setdiff(ratio.bed,cnvs.bed)
normals <- subsetByOverlaps(ratio.bed,normals)

#mu <- mean(score(normals),na.rm=TRUE)
#sigma<- sd(score(normals),na.rm=TRUE)

#hist(score(normals),n=500,xlim=c(0,2))
#hist(log(score(normals)),n=500,xlim=c(-1,1))

#shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])
#qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))
#qqline(score(normals)[which(!is.na(score(normals)))], col = 2)

#shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])
#qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))
#qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)

numberOfCol=length(cnvs)

for (i in c(1:length(cnvs[,1]))) {
  values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))
  #wilcox.test(values,mu=mu)
  W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)
                                if(class(resultw)=="try-error") return(list("statistic"=NA,"parameter"=NA,"p.value"=NA,"null.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultw}
  KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)
                                 if(class(resultks)=="try-error") return(list("statistic"=NA,"p.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultks}
  #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)
  #  if(class(resultks)=="try-error") NA) else resultks
  cnvs[i,numberOfCol+1]=W(values,normals)$p.value
  cnvs[i,numberOfCol+2]=KS(values,normals)$p.value
}

if (numberOfCol==5) {
  names(cnvs)=c("chr","start","end","copy number","status","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==7) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}
if (numberOfCol==9) {
  names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","somatic/germline","precentageOfGermline","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
}

write.table(cnvs, file=output,sep="\t",quote=F,row.names=F)