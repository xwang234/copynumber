#! /usr/bin/env Rscript
#SBATCH -t 1-1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

library(GenomicRanges)
library("biomaRt")
mart=useMart("ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #GRCh37
findgenes=function(chr,start,end)
{
  attributes <- c("hgnc_symbol")
  filters <- c("chromosome_name","start","end")
  if (grepl('chr',chr))
    chr=substr(chr,4,nchar(as.character(chr)))
  
  values=list(chromosome=chr,start=start,end=end)
  genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
  genes=unique(genes[,1])
  idxkeep=which(genes!="")
  genes=genes[idxkeep]
  
  genes1=NULL
  if (length(genes)>0)
  {
    if (length(genes)==1)
    {
      genes1=genes
    }else
    {
      for (j in 1:length(genes))
      {
        if (j==1)
        {
          genes1=genes[j]
        }
        else
        {
          genes1=paste(genes1,genes[j],sep=",")
        }
      }
    }
  }
  if (is.null(genes1)) genes1=NA
  return(genes1)
}


setwd("/fh/scratch/delete30/dai_j/henan/absolute")
args <- commandArgs(trailingOnly = TRUE)
numpair=as.integer(args[1])
#numpair=1
normals=paste0((c(3,11,13,15,17,25,29,33,37,41)+1),"A")
tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
normal=normals[numpair]
tumor=tumors[numpair]
print(tumor)
library(ABSOLUTE1)
cbsfile=paste0("/fh/scratch/delete30/dai_j/henan/gistic/",tumor,".cnv.freecseg")
cbstable=read.table(cbsfile,header=F)
colnames(cbstable)=c("ID","Chromosome","Start","End","Num_Probes","Segment_Mean")
cbstable=cbstable[is.finite(cbstable$Segment_Mean),]
cbstable=cbstable[!is.na(cbstable$Segment_Mean),]
cnratio=2^cbstable$Segment_Mean
len=cbstable$End-cbstable$Start
weight=len/sum(as.numeric(len))
cnratio1=cnratio*weight
cnratio=cnratio/sum(cnratio1)
cbstable$Segment_Mean=log2(cnratio)

seg.dat.fn=paste0(tumor,".cnv.freec.txt")
write.table(cbstable,file=seg.dat.fn,col.names=T,row.names=F,sep="\t",quote=F)

maf.fn=paste0(tumor,".maf.txt")

#if (! file.exists(maf.fn))
#{
  chrs=c(1:22,"X","Y")
  mutectfile=paste0("/fh/scratch/delete30/dai_j/henan/mutect/",tumor,".Mutect_out_keep.txt")
  mutecttable=read.table(mutectfile,header=T,stringsAsFactors = F,sep="\t",quote="")
  colnames(mutecttable)[1]="Chromosome"
  mutecttable[,1]=gsub("chr","",mutecttable[,1],fixed = T)
  mutecttable=mutecttable[mutecttable[,1] %in% chrs,]
  colnames(mutecttable)[2]="Start_position"
  colnames(mutecttable)[which(colnames(mutecttable)=="dbsnp_site")]="dbSNP_Val_Status"
  if (! "Hugo_Symbol" %in% colnames(mutecttable)) mutecttable$Hugo_Symbol=rep("",nrow(mutecttable))
  if (! "Tumor_Sample_Barcode" %in% colnames(mutecttable)) mutecttable$Tumor_Sample_Barcode=rep(tumor,nrow(mutecttable))
  for (i in 1:nrow(mutecttable))
  {
    mutecttable$Hugo_Symbol[i]=findgenes(mutecttable[i,1],mutecttable[i,2],mutecttable[i,2])
  }
  write.table(mutecttable,maf.fn,col.names = T,row.names = F,sep="\t",quote=F)
#}  




#sigma.p=0.0
sigma.p=0.05
#max.sigma.h=0.05
max.sigma.h=0.5
min.ploidy=0.1
max.ploidy=5
primary.disease="EAC"
platform="Illumina_WES"

sample.name=tumor
results.dir=paste0("./",tumor)
max.as.seg.count=5000
max.non.clonal=10
max.neg.genome=10
copy_num_type="total"
output.fn.base=tumor
verbose=TRUE
RunAbsolute(seg.dat.fn=seg.dat.fn, sigma.p=sigma.p, max.sigma.h=max.sigma.h, min.ploidy=min.ploidy, max.ploidy=max.ploidy, 
                primary.disease=primary.disease, platform=platform, sample.name=sample.name, results.dir=results.dir,
                max.as.seg.count=max.as.seg.count, max.non.clonal=max.non.clonal, max.neg.genome=max.neg.genome,
                copy_num_type=copy_num_type, maf.fn = maf.fn, min.mut.af = 0.1, 
                output.fn.base=output.fn.base, verbose=verbose)

