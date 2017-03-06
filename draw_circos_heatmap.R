#!/usr/bin/env Rscript
library(GenomicRanges)
library("biomaRt")
library(RCircos)
# data(UCSC.HG19.Human.CytoBandIdeogram)
# chr.exclude <- NULL
# chr.exclude <- c("chrX", "chrY")
# cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
# #how many tracks inside:
# tracks.inside <- 3
# tracks.outside <- 0
# RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
# 
# # rcircos.params <- RCircos.Get.Plot.Parameters()
# # rcircos.cyto <- RCircos.Get.Plot.Ideogram()
# # rcircos.position <- RCircos.Get.Plot.Positions()
# #RCircos.List.Parameters()
# rcircos.params <- RCircos.Get.Plot.Parameters()
# #rcircos.params$base.per.unit <- 30000
# rcircos.params$track.in.start = 2
# RCircos.Reset.Plot.Parameters(rcircos.params)
# 
# readsegsfromgistic=function(segfile)
# {
#   gistictable=read.table(segfile,header=F,sep="\t")
#   idxNA=is.na(gistictable[1,])
#   gistictable=gistictable[,!idxNA]
#   res=data.frame(matrix(NA,nrow=ncol(gistictable)-1,ncol=8))
#   colnames(res)=c('chr','start','end','cytoband','qvalue','residqvalue','numgenes','genes')
#   for (i in 2:ncol(gistictable))
#   {
#     tmp=as.character(gistictable[4,i])
#     tmp1=unlist(strsplit(tmp,":"))
#     tmp2=unlist(strsplit(tmp1[2],"-"))
#     res[i-1,1]=tmp1[1]
#     res[i-1,2]=tmp2[1]
#     res[i-1,3]=tmp2[2]
#     res[i-1,4]=as.character(gistictable[1,i])
#     res[i-1,5]=as.character(gistictable[2,i])
#     res[i-1,6]=as.character(gistictable[3,i])
#     genes=as.character(gistictable[5,i])
#     n=1
#     for (j in 6:nrow(gistictable))
#     {
#       gene=as.character(gistictable[j,i])
#       if (gene!="")
#       {
#         n=n+1
#         genes=paste0(genes,',',gene)
#       }
#     }
#     res[i-1,7]=n
#     res[i-1,8]=genes
#   }
#   res[,1]=gsub("chr","",res[,1],fixed=T)
#   res[,2]=as.integer(res[,2])
#   res[,3]=as.integer(res[,3])
#   chrs=c(1:22,"X","Y")
#   res1=NULL
#   for (i in 1:length(chrs))
#   {
#     tmp=res[res[,1]==chrs[i],]
#     res1=rbind(res1,tmp)
#   }
#   res1=cbind.data.frame(res1,addnum=rep(1,nrow(res1)))
# }
# 
# #col column determines color
# RCircos.Heatmap.Plot1=function(heatmap.data, track.num, side="in") 
# {
#   RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
#   RCircos.Pos <- RCircos.Get.Plot.Positions()
#   RCircos.Par <- RCircos.Get.Plot.Parameters()
#   heatmap.data <- RCircos.Get.Plot.Data(heatmap.data, "plot")
#   #ColorRamp <- RCircos.Get.Heatmap.ColorScales(RCircos.Par$heatmap.color)
#   #columns <- 5:(ncol(heatmap.data) - 1)
#   #min.value <- min(as.matrix(heatmap.data[, columns]))
#   #max.value <- max(as.matrix(heatmap.data[, columns]))
#   #min.value=max.value=1
#   #ColorLevel <- seq(min.value, max.value, length = length(ColorRamp))
#   heatmap.locations <- as.numeric(heatmap.data[, ncol(heatmap.data)])
#   start <- heatmap.locations - RCircos.Par$heatmap.width/2
#   end <- heatmap.locations + RCircos.Par$heatmap.width/2
#   data.chroms <- as.character(heatmap.data[, 1])
#   chromosomes <- unique(data.chroms)
#   cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
#   for (a.chr in 1:length(chromosomes)) {
#     cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
#     locations <- as.numeric(RCircos.Cyto$Location[cyto.rows])
#     chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]]
#     chr.end <- max(locations)
#     data.rows <- which(data.chroms == chromosomes[a.chr])
#     start[data.rows[start[data.rows] < chr.start]] <- chr.start
#     end[data.rows[end[data.rows] > chr.end]] <- chr.end
#   }
#   locations <- RCircos.Track.Positions(side, track.num)
#   out.pos <- locations[1]
#   in.pos <- locations[2]
#   chroms <- unique(RCircos.Cyto$Chromosome)
#   for (a.chr in 1:length(chroms)) {
#     the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
#                             ]
#     the.start <- the.chr$Location[1] - the.chr$Unit[1] + 
#       1
#     the.end <- the.chr$Location[nrow(the.chr)]
#     polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
#                    RCircos.Pos[the.end:the.start, 1] * in.pos)
#     polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
#                    RCircos.Pos[the.end:the.start, 2] * in.pos)
#     polygon(polygon.x, polygon.y, col = "white")
#   }
#   #heatmap.value <- as.numeric(heatmap.data[, data.col])
#   for (a.point in 1:nrow(heatmap.data)) {
#     #the.level <- which(ColorLevel >= heatmap.value[a.point])
#     cell.color <- heatmap.data$col[a.point]
#     the.start <- start[a.point]
#     the.end <- end[a.point]
#     polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
#                    RCircos.Pos[the.end:the.start, 1] * in.pos)
#     polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
#                    RCircos.Pos[the.end:the.start, 2] * in.pos)
#     polygon(polygon.x, polygon.y, col = cell.color, border = cell.color,lwd=2)
#   }
# }
# 
# gisticdir1='/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx0_conf0.95_armpeel0_brlen0.98_broad1'
# segfile1=paste0(gisticdir1,"/","amp_genes.conf_95.txt")
# 
# gisticdir2='/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx0_conf0.95_armpeel0_brlen0.98_broad1'
# segfile2=paste0(gisticdir2,"/","amp_genes.conf_95.txt")
# 
# gisticdir3='/fh/scratch/delete30/dai_j/escc/gistic/escc_ploid2degree3force0_cnv_rx0_conf0.95_armpeel0_brlen0.98_broad1'
# segfile3=paste0(gisticdir3,"/","amp_genes.conf_95.txt")
# 
# segdata1=readsegsfromgistic(segfile1)
# segdata2=readsegsfromgistic(segfile2)
# segdata3=readsegsfromgistic(segfile3)
# 
# 
# layout(matrix(data=c(1,1), nrow=1, ncol=1), widths=4, heights=4)
# RCircos.Set.Plot.Area()
# RCircos.Chromosome.Ideogram.Plot()
# side <- "in"
# RCircos.Heatmap.Plot1(segdata1, track.num=1, side,col="red")
# RCircos.Heatmap.Plot1(segdata2, track.num=2, side,col="red")
# RCircos.Heatmap.Plot1(segdata3, track.num=3, side,col="red")
# 
# 
# segfile1=paste0(gisticdir1,"/","del_genes.conf_95.txt")
# segfile2=paste0(gisticdir2,"/","del_genes.conf_95.txt")
# segfile3=paste0(gisticdir3,"/","del_genes.conf_95.txt")
# segdata1=readsegsfromgistic(segfile1)
# segdata2=readsegsfromgistic(segfile2)
# segdata3=readsegsfromgistic(segfile3)
# RCircos.Set.Plot.Area()
# RCircos.Chromosome.Ideogram.Plot()
# RCircos.Heatmap.Plot1(segdata1, track.num=1, side)
# RCircos.Heatmap.Plot1(segdata2, track.num=2, side)
# RCircos.Heatmap.Plot1(segdata3, track.num=3, side)

##--------------------the updated implementation
findcytoband=function(chr,start,end)
{
  cytobands=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
  cytobands=cytobands[cytobands$chrom==chr,]
  gr_cytobands=GRanges(seqnames=cytobands$chrom,ranges=IRanges(start=cytobands$start,end=cytobands$end),name=cytobands$name)
  gr_segment=GRanges(seqnames=chr,ranges=IRanges(start=start,end=end))
  res=NA
  maxlen=0
  for (i in 1:length(gr_cytobands))
  {
    olap=intersect(gr_cytobands[i,],gr_segment)
    if (length(olap)>0)
    {
      len=width(olap)
      if (len>maxlen)
      {
        res=as.character(mcols(gr_cytobands[i,])["name"][[1]])
        maxlen=len
      }
    }
  }
  if (!is.na(res))
  {
    res=paste0(gsub("chr","",chr),res)
  }
  return(res)
}

getwholetable1=function(alllessionsfile,opt="AMP",opt1="peak",opt2="binary")
{
  
  alllessiontable=read.table(alllessionsfile,header=T,sep="\t",stringsAsFactors=F)
  #remove the last NA column introduced by \t
  idx=which(colSums(is.na(alllessiontable))==nrow(alllessiontable))
  alllessiontable=alllessiontable[,-idx]
  idx=which(colnames(alllessiontable)=="Amplitude.Threshold")
  idx=idx+1 # tumor samples start
  idxkeep=which(rowSums(abs(alllessiontable[1:(nrow(alllessiontable)/2),idx:ncol(alllessiontable)]),na.rm = T)>1)
  if (opt2=="binary")
  {
    alllessiontable=alllessiontable[1:(nrow(alllessiontable)/2),]
  }else
  {
    alllessiontable=alllessiontable[(nrow(alllessiontable)/2+1):nrow(alllessiontable),]
  }
  alllessiontable=alllessiontable[idxkeep,]
  alllessiontable$Descriptor=gsub(" ","",alllessiontable$Descriptor,fixed=T)
  idx1=grepl("Amplification",alllessiontable$Unique.Name)
  if (opt=="AMP")
  {
    alllessiontable=alllessiontable[idx1,]
  }else
  {
    alllessiontable=alllessiontable[!idx1,]
  }

  res=data.frame(matrix(NA,nrow=nrow(alllessiontable),ncol=5))
  colnames(res)=c("cytoband","chr","start","end","qvalue")
  res$chr=as.character(res$chr)
  res$start=as.integer(res$start)
  res$end=as.integer(res$end)
  res$cytoband=alllessiontable$Descriptor
  res$qvalue=alllessiontable$Residual.q.values.after.removing.segments.shared.with.higher.peaks
  res=cbind.data.frame(res,alllessiontable[,idx:ncol(alllessiontable)])
  for (i in 1:nrow(res))
  {
    #tmp=unlist(strsplit(alllessiontable$Wide.Peak.Limits[i],":",fixed=T))
    if (opt1=="peak")
    {
      tmp=unlist(strsplit(alllessiontable$Peak.Limits[i],":",fixed=T))
    }else
    {
      tmp=unlist(strsplit(alllessiontable$Region.Limits[i],":",fixed=T))
    }
    
    res$chr[i]=tmp[1]
    tmp=unlist(strsplit(tmp[2],"(",fixed=T))
    tmp=unlist(strsplit(tmp[1],"-",fixed=T))
    res$start[i]=as.integer(tmp[1])
    res$end[i]=as.integer(tmp[2])
  }
#   gr_res=GRanges(seqnames=res$chr,ranges=IRanges(start=res$start,end=res$end),cytoband=res$cytoband)
#   gr_res=reduce(gr_res)
#   res=data.frame(matrix(NA,nrow=length(gr_res),ncol=4))
#   colnames(res)=c("cytoband","chr","start","end")
#   res$chr=seqnames(gr_res)
#   res$start=start(gr_res)
#   res$end=end(gr_res)
#   res$cytoband=as.character(res$cytoband)
#   res$chr=as.character(res$chr)
#   for (i in 1:nrow(res))
#   {
#     res$cytoband[i]=findcytoband(chr=res$chr[i],start=res$start[i],end=res$end[i])
#   }
  #add rownames
  idxdup=duplicated(res$cytoband)
  rownames(res)[!idxdup]=res$cytoband[!idxdup]
  dupsign=0
  if (sum(idxdup)==1)
  {
    rownames(res)[idxdup]=paste0(res$cytoband[idxdup],"(1)")
  }
  if (sum(idxdup)>1)
  {
    idxdup=which(idxdup==T)
    rownames(res)[idxdup[1]]=paste0(res$cytoband[idxdup[1]],"(1)")
    dupsign=1
    for (i in 2:length(idxdup))
    {
      if(res$cytoband[idxdup[i]]==res$cytoband[idxdup[i-1]])
      {
        dupsign=dupsign+1
        rownames(res)[idxdup[i]]=paste0(res$cytoband[idxdup[i]],"(",dupsign,")")
      }else
      {
        dupsign=1
        rownames(res)[idxdup[i]]=paste0(res$cytoband[idxdup[i]],"(1)")
      }
    }
  }
  
  return(res)
}


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


overlapsegdata=function(segdata1,segdata2,segdata3)
{
  segdata=rbind.data.frame(segdata1[,1:5],segdata2[,1:5],segdata3[1:5])
  chrs=paste0("chr",c(1:22,"X","Y"))
  segdata_=NULL
  for (i in 1:length(chrs))
  {
    tmp=segdata[segdata$chr==chrs[i],]
    tmp=tmp[order(tmp$start),]
    segdata_=rbind(segdata_,tmp)
  }
  
  cytobands=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
  gr_cytobands=GRanges(seqnames=cytobands$chrom,ranges=IRanges(start=cytobands$start,end=cytobands$end),name=cytobands$name)
  
  allsegs=data.frame(matrix(NA,nrow=0,ncol=6))
  colnames(allsegs)=c("chr","start","end","cytoband","gene")
  allsegs$gene=as.character(allsegs$gene)
  uniqchrs=unique(segdata_$chr)
  for (j in 1:length(uniqchrs))
  {
    idx=which(segdata_$chr==uniqchrs[j])
    if (length(idx)>0)
    {
      tmptable=segdata_[idx,c("chr","start","end")]
      ends=unique(c(tmptable$start,tmptable$end))
      ends=ends[order(ends)]
      for (k in 1:(length(ends)-1))
      {
        gr_seg=GRanges(seqnames=uniqchrs[j],ranges=IRanges(start=ends[k],end=ends[k+1]))
        olap=subsetByOverlaps(gr_cytobands,gr_seg)
        if (length(olap)>0)
        {
          idx1=which.max(width(olap))
          tmpcytoband=as.character(mcols(olap)["name"][,1])[idx1]
        }
        tmpgene=findgenes(chr=uniqchrs[j],start=ends[k],end=ends[k+1])
        tmpseg=data.frame(chr=uniqchrs[j],start=ends[k],end=ends[k+1],cytoband=tmpcytoband,gene=tmpgene)
        allsegs=rbind.data.frame(allsegs,tmpseg)
      }
    }
  }
  allsegscolnames=c(colnames(allsegs),colnames(segdata1)[6:ncol(segdata1)],colnames(segdata2)[6:ncol(segdata2)],colnames(segdata3)[6:ncol(segdata3)])
  allsegs=cbind.data.frame(allsegs,data.frame(matrix(NA,nrow=nrow(allsegs),ncol=ncol(segdata1)-5+ncol(segdata2)-5+ncol(segdata3)-5)))
  colnames(allsegs)=allsegscolnames
  
  gr_allsegs=GRanges(seqnames=allsegs$chr,ranges=IRanges(start=allsegs$start,end=allsegs$end))
  gr_segdata1=GRanges(seqnames=segdata1$chr,ranges=IRanges(start=segdata1$start,end=segdata1$end))
  gr_segdata2=GRanges(seqnames=segdata2$chr,ranges=IRanges(start=segdata2$start,end=segdata2$end))
  gr_segdata3=GRanges(seqnames=segdata3$chr,ranges=IRanges(start=segdata3$start,end=segdata3$end))
  for (i in 1:length(gr_allsegs))
  {
    olap1=subsetByOverlaps(gr_segdata1,gr_allsegs[i,],minoverlap=10)
    if (length(olap1)>0)
    {
      idx=which(segdata1$chr==seqnames(olap1)[1] & segdata1$start==start(olap1)[1] & segdata1$end==end(olap1)[1])
      allsegs[i,colnames(segdata1)[6:ncol(segdata1)]]=segdata1[idx,6:ncol(segdata1)]
    }
    olap2=subsetByOverlaps(gr_segdata2,gr_allsegs[i,],minoverlap=10)
    if (length(olap2)>0)
    {
      idx=which(segdata2$chr==seqnames(olap2)[1] & segdata2$start==start(olap2)[1] & segdata2$end==end(olap2)[1])
      allsegs[i,colnames(segdata2)[6:ncol(segdata2)]]=segdata2[idx,6:ncol(segdata2)]
    }
    olap3=subsetByOverlaps(gr_segdata3,gr_allsegs[i,],minoverlap=10)
    if (length(olap3)>0)
    {
      idx=which(segdata3$chr==seqnames(olap3)[1] & segdata3$start==start(olap3)[1] & segdata3$end==end(olap3)[1])[1]
      allsegs[i,colnames(segdata3)[6:ncol(segdata3)]]=segdata3[idx,6:ncol(segdata3)]
    }
  }
  noidx=is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]]))
  allsegs=allsegs[!noidx,]
  idx1=which(! is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx2=which(is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & ! is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx3=which(is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & ! is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx12=which(! is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & ! is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx13=which(! is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & ! is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx23=which(is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & ! is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & ! is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  idx123=which(! is.na(rowSums(allsegs[,colnames(segdata1)[6:ncol(segdata1)]])) & ! is.na(rowSums(allsegs[,colnames(segdata2)[6:ncol(segdata2)]])) & ! is.na(rowSums(allsegs[,colnames(segdata3)[6:ncol(segdata3)]])))
  segs1=allsegs[idx1,]
  segs2=allsegs[idx2,]
  segs3=allsegs[idx3,]
  segs12=allsegs[idx12,]
  segs13=allsegs[idx13,]
  segs23=allsegs[idx23,]
  segs123=allsegs[idx123,]
  result=list(allsegs=allsegs,segs1=segs1,segs2=segs2,segs3=segs3,segs12=segs12,segs13=segs13,segs23=segs23,segs123=segs123)
}

countoverlap2segs=function(segdata1,segdata2)
{
  cytobands=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
  gr_cytobands=GRanges(seqnames=cytobands$chrom,ranges=IRanges(start=cytobands$start,end=cytobands$end),name=cytobands$name)
  segdata=rbind.data.frame(segdata1[,1:4],segdata2[,1:4])
  chrs=paste0("chr",c(1:22,"X","Y"))
  segdata_=NULL
  for (i in 1:length(chrs))
  {
    tmp=segdata[segdata$chr==chrs[i],]
    tmp=tmp[order(tmp$start),]
    segdata_=rbind(segdata_,tmp)
  }
  gr_segdata1=GRanges(seqnames=segdata1$chr,ranges=IRanges(start=segdata1$start,end=segdata1$end),cytoband=segdata1$cytoband)
  gr_segdata1=reduce(gr_segdata1)
  gr_segdata2=GRanges(seqnames=segdata2$chr,ranges=IRanges(start=segdata2$start,end=segdata2$end),cytoband=segdata2$cytoband)
  gr_segdata2=reduce(gr_segdata2)
  gr_segdata=GRanges(seqnames=segdata_$chr,ranges=IRanges(start=segdata_$start,end=segdata_$end),cytoband=segdata_$cytoband)
  gr_segdata=reduce(gr_segdata)
  overlapseg=data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(overlapseg)=c("cytoband","chr","start","end")
  overlapseg$cytoband=as.character(overlapseg$cytoband)
  overlapseg$chr=as.character(overlapseg$chr)
  for (i in 1:length(gr_segdata1))
  {
    olap=intersect(gr_segdata1[i,],gr_segdata2)
    if (length(olap)>0)
    {
      for (j in 1:length(olap))
      {
        chr=as.character(seqnames(olap[j,]))
        olap1=subsetByOverlaps(gr_cytobands,olap[j,])
        if (length(olap1)>0)
        {
          idx1=which.max(width(olap1))
          tmpcytoband=paste0(substr(chr,4,nchar(chr)),as.character(mcols(olap1)["name"][,1])[idx1])
          overlapseg=rbind(overlapseg,data.frame(cytoband=tmpcytoband,chr=chr,start=start(olap[j,]),end=end(olap[j,])))
        }
      }
    }
  }
  num_segdata1=num_segdata2=num_overlapseg=0
  print("overlapsegs:")
  for (i in 1:length(gr_segdata))
  {
    olap1=intersect(gr_segdata[i,],gr_segdata1)
    olap2=intersect(gr_segdata[i,],gr_segdata2)
    if (length(olap1)>0) num_segdata1=num_segdata1+1
    if (length(olap2)>0) num_segdata2=num_segdata2+1
    if (length(olap1)>0 & length(olap2)>0)
      {
        num_overlapseg=num_overlapseg+1
        print(gr_segdata[i,])
      }
  }
  
  width_segdata=as.numeric(format(sum(width(gr_segdata))/10^6,digits = 3,nsmall=1))
  width_segdata1=as.numeric(format(sum(width(gr_segdata1))/10^6,digits = 3,nsmall=1))
  width_segdata2=as.numeric(format(sum(width(gr_segdata2))/10^6,digits = 3,nsmall=1))
  gr_overlapseg=GRanges(seqnames=overlapseg$chr,ranges=IRanges(start=overlapseg$start,end=overlapseg$end),cytoband=overlapseg$cytoband)
  width_overlapseg=as.numeric(format(sum(width(reduce(gr_overlapseg)))/10^6,digits = 3,nsmall=1))
  result=list(allsegs=gr_segdata,overlapseg=overlapseg,width_segdata=width_segdata,width_segdata1=width_segdata1,width_segdata2=width_segdata2,width_overlapseg=width_overlapseg,
              num_segdata1=num_segdata1,num_segdata2=num_segdata2,num_overlapseg=num_overlapseg)
}

segonlyinseg1=function(segdata1,segdata2)
{
  gr_segdata1=GRanges(seqnames=segdata1$chr,ranges=IRanges(start=segdata1$start,end=segdata1$end),cytoband=segdata1$cytoband)
  #gr_segdata1=reduce(gr_segdata1)
  gr_segdata2=GRanges(seqnames=segdata2$chr,ranges=IRanges(start=segdata2$start,end=segdata2$end),cytoband=segdata2$cytoband)
  #gr_segdata2=reduce(gr_segdata2)
  seg1only=data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(seg1only)=c("cytoband","chr","start","end")
  for (i in 1:length(gr_segdata1))
  {
    olap=subsetByOverlaps(gr_segdata1[i,],gr_segdata2)
    if (length(olap)==0)
    {
      #chr=as.character(seqnames(gr_segdata1[i,]))
      #start=start(gr_segdata1[i,])
      #end=end(gr_segdata1[i,])
      #tmp=data.frame(cytoband=findcytoband(chr=chr,start=start,end=end),chr=chr,start=start,end=end)
      tmp=data.frame(cytoband=segdata1$cytoband[i],chr=segdata1$chr[i],start=segdata1$start[i],end=segdata1$end[i])
      seg1only=rbind.data.frame(seg1only,tmp)
    }
  }
  return(seg1only)
}

printoverlap2segs=function(overlapsegs)
{
  print(paste0("seg1:",overlapsegs$width_segdata1))
  print(paste0("seg2:",overlapsegs$width_segdata2))
  print(paste0("total width:",overlapsegs$width_segdata))
  print(paste0("overlap:",overlapsegs$width_overlapseg))
  print(paste0("overlap ratio:",overlapsegs$width_overlapseg/overlapsegs$width_segdata))
  print(unique(as.character(olapseg12$overlapseg$cytoband)))
}

RCircos.Get.Plot.Data1=function (genomic.data, plot.type) 
{
  genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type)
  #data.points <- rep(0, nrow(genomic.data))
  data.points=data.frame(matrix(0,nrow=nrow(genomic.data),ncol=2))
  for (a.row in 1:nrow(genomic.data)) {
    chromosome <- as.character(genomic.data[a.row, "chr"])
    data.points[a.row,1] <- RCircos.Data.Point(chromosome, genomic.data[a.row,"start"]) #for start
    data.points[a.row,2] <- RCircos.Data.Point(chromosome, genomic.data[a.row,"end"])
  }
  genomic.data["Location_start"] <- data.points[,1]
  genomic.data["Location_end"] <- data.points[,2]
  genomic.data <- genomic.data[order(genomic.data$Location_start), 
                               ]
  return(genomic.data)
}
#use the actual segment width
RCircos.Heatmap.Plot2=function(heatmap.data, track.num, side="in") 
{
  heatmap.data=heatmap.data[,c("chr","start","end","col")]
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  heatmap.data <- RCircos.Get.Plot.Data1(heatmap.data, "plot")
  #ColorRamp <- RCircos.Get.Heatmap.ColorScales(RCircos.Par$heatmap.color)
  #columns <- 5:(ncol(heatmap.data) - 1)
  #min.value <- min(as.matrix(heatmap.data[, columns]))
  #max.value <- max(as.matrix(heatmap.data[, columns]))
  #min.value=max.value=1
  #ColorLevel <- seq(min.value, max.value, length = length(ColorRamp))
  #heatmap.locations <- as.numeric(heatmap.data[, ncol(heatmap.data)])
  start=heatmap.data$Location_start
  end=heatmap.data$Location_end
  for (i in 1:length(start))
  {
    if (start[i]==end[i])
    {
      start[i]=start[i]-RCircos.Par$heatmap.width/2
      end[i]=end[i]+RCircos.Par$heatmap.width/2
    }
  }
  #start <- heatmap.locations - RCircos.Par$heatmap.width/2
  #end <- heatmap.locations + RCircos.Par$heatmap.width/2
  data.chroms <- as.character(heatmap.data[, 1])
  chromosomes <- unique(data.chroms)
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chromosomes)) {
    cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
    locations <- as.numeric(RCircos.Cyto$Location[cyto.rows])
    chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]]
    chr.end <- max(locations)
    data.rows <- which(data.chroms == chromosomes[a.chr])
    start[data.rows[start[data.rows] < chr.start]] <- chr.start
    end[data.rows[end[data.rows] > chr.end]] <- chr.end
  }
  locations <- RCircos.Track.Positions(side, track.num)
  out.pos <- locations[1]
  in.pos <- locations[2]
  chroms <- unique(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chroms)) {
    the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
                            ]
    the.start <- the.chr$Location[1] - the.chr$Unit[1] + 
      1
    the.end <- the.chr$Location[nrow(the.chr)]
    polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
                   RCircos.Pos[the.end:the.start, 1] * in.pos)
    polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
                   RCircos.Pos[the.end:the.start, 2] * in.pos)
    polygon(polygon.x, polygon.y, col = "white")
  }
  #heatmap.value <- as.numeric(heatmap.data[, data.col])
  for (a.point in 1:nrow(heatmap.data)) {
    #the.level <- which(ColorLevel >= heatmap.value[a.point])
    cell.color <- heatmap.data$col[a.point]
    the.start <- start[a.point]
    the.end <- end[a.point]
    polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
                   RCircos.Pos[the.end:the.start, 1] * in.pos)
    polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
                   RCircos.Pos[the.end:the.start, 2] * in.pos)
    polygon(polygon.x, polygon.y, col = cell.color, border = cell.color,lwd=2)
  }
}

printsegdata=function(segdata,output=NULL)
{
  if (! "gene" %in% colnames(segdata))
  {
    segdata=cbind.data.frame(segdata,gene=rep("",nrow(segdata)))
    segdata$gene=as.character(segdata$gene)
    for (i in 1:nrow(segdata))
    {
      segdata$gene[i]=findgenes(chr=segdata$chr[i],start=segdata$start[i],end=segdata$end[i])
    }
  }
  
  segdata$start[segdata$start<1000]=0
  segdata$start=format(segdata$start/10^6,nsmall=2,digits = 2)
  segdata$end[segdata$end<1000]=0
  segdata$end=format(segdata$end/10^6,nsmall=2,digits = 2)
  if (is.null(output))
  {
    print(segdata)
  }else
  {
    write.table(segdata,file=output,sep="\t",row.names=F,quote=F)  
  }
}

data(UCSC.HG19.Human.CytoBandIdeogram)
#chr.exclude <- NULL
chr.exclude <- c("chrX", "chrY")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
#how many tracks inside:
tracks.inside <- 3
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
#rcircos.params$base.per.unit <- 30000
rcircos.params$track.in.start = 2
RCircos.Reset.Plot.Parameters(rcircos.params)

gisticdir1="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile1=paste0(gisticdir1,"/all_lesions.conf_95.txt")
wgstumors1=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
             "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868","SRR1001635")
wgstumors2=paste0("X",c(3,11,13,15,17,25,29,33,37,41),"A")
gisticdir2="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile2=paste0(gisticdir2,"/all_lesions.conf_95.txt")
wgstumors3=paste0("T",c(1:6,8:18))
wgstumors3=paste0("T",c(1:4,6,8:18))
gisticdir3="/fh/scratch/delete30/dai_j/escc/gistic/escc_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile3=paste0(gisticdir3,"/all_lesions.conf_95.txt")

ampsegdata1=getwholetable1(alllessionsfile1,opt="AMP")
#printsegdata(ampsegdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_ampseg.txt")
ampsegdata2=getwholetable1(alllessionsfile2,opt="AMP")
#printsegdata(ampsegdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_ampseg.txt")
ampsegdata3=getwholetable1(alllessionsfile3,opt="AMP")
#printsegdata(ampsegdata3,output="/fh/scratch/delete30/dai_j/escc/gistic/escc_gistic_ampseg.txt")
ampolapseg12=countoverlap2segs(ampsegdata1,ampsegdata2)
#printsegdata(olapseg12$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_henan_gistic_ampoverlapseg.txt")
printoverlap2segs(ampolapseg12)
ampolapseg13=countoverlap2segs(ampsegdata1,ampsegdata3)
#printsegdata(olapseg13$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_escc_gistic_ampoverlapseg.txt")
printoverlap2segs(ampolapseg13)
ampolapseg23=countoverlap2segs(ampsegdata2,ampsegdata3)
#printsegdata(olapseg23$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/henan_escc_gistic_ampoverlapseg.txt")
printoverlap2segs(ampolapseg23)
olapseg123=countoverlap2segs(olapseg12$overlapseg,ampsegdata3)
#printsegdata(olapseg123$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_henan_escc_gistic_ampoverlapseg.txt")
printoverlap2segs(olapseg123)
#unique segs:
uniqampseg1=segonlyinseg1(ampsegdata1,ampsegdata2)
printsegdata(uniqampseg1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_ampdistinctseg.txt")
uniqampseg2=segonlyinseg1(ampsegdata2,ampsegdata1)
printsegdata(uniqampseg2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_ampdistinctseg.txt")
ampseg2=segonlyinseg1(ampseg2,ampsegdata3)
ampseg3=segonlyinseg1(ampsegdata3,ampsegdata2)
printsegdata(ampseg3,output="/fh/scratch/delete30/dai_j/escc/gistic/escc_gistic_ampdistinctseg.txt")
ampseg3=segonlyinseg1(ampseg3,ampsegdata1)



library(VennDiagram)
grid.newpage()
draw.triple.venn(area1 =olapseg12$width_segdata1, area2 = olapseg12$width_segdata2, area3 = olapseg13$width_segdata2, n12 = olapseg12$width_overlapseg, 
                 n23 = olapseg23$width_overlapseg, n13 = olapseg13$width_overlapseg, n123 = olapseg123$width_overlapseg, category = c("US-EA", "CH-EA", "CH-ESCC"), lty = "blank", 
                 fill = c("blue", "green", "red"),cex=rep(2,7),cat.cex = rep(2, 3),cat.col=c("blue", "green", "red"),euler.d=T)
#ampres=overlapsegdata(segdata1,segdata2,segdata3)
#load("tmp.RData") #where the results (ampres,delres) were saved

addcolor4segdata=function(segdata,segres,col1,col2)
{
  gr_segdata=GRanges(seqnames=segdata$chr,ranges=IRanges(start=segdata$start,end=segdata$end))
  commseg=rbind.data.frame(segres$segs12,segres$segs13,segres$segs23)
  gr_commseg=GRanges(seqnames=commseg$chr,ranges=IRanges(start=commseg$start,end=commseg$end))
  if ( !"col" %in% colnames(segdata))
    segdata=cbind.data.frame(segdata,data.frame(col=rep(col1,nrow(segdata)))) #add c color column
  segdata$col=as.character(segdata$col)
  for (i in 1:nrow(segdata))
  {
    olap=subsetByOverlaps(gr_commseg,gr_segdata[i,])
    if (length(olap)>0)
    {
      segdata$col[i]=col2
    }
  }
  return(segdata)
}
# segdata1=addcolor4segdata(segdata=segdata1,segres=ampres,col1="blue",col2="red")
# segdata1=segdata1[,-1]
# segdata2=addcolor4segdata(segdata=segdata2,segres=ampres,col1="blue",col2="red")
# segdata2=segdata2[,-1]
# segdata3=addcolor4segdata(segdata=segdata3,segres=ampres,col1="blue",col2="red")
# segdata3=segdata3[,-1]
ampsegdata1=cbind.data.frame(ampsegdata1,col=rep("red",nrow(ampsegdata1)))
ampsegdata1$col=as.character(ampsegdata1$col)
ampsegdata1=ampsegdata1[,-1]
ampsegdata2=cbind.data.frame(ampsegdata2,col=rep("red",nrow(ampsegdata2)))
ampsegdata2$col=as.character(ampsegdata2$col)
ampsegdata2=ampsegdata2[,-1]
ampsegdata3=cbind.data.frame(ampsegdata3,col=rep("red",nrow(ampsegdata3)))
ampsegdata3$col=as.character(ampsegdata3$col)
ampsegdata3=ampsegdata3[,-1]


RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(ampsegdata1, track.num=1)
RCircos.Heatmap.Plot2(ampsegdata2, track.num=2)
RCircos.Heatmap.Plot2(ampsegdata3, track.num=3)


opt="DEL"
delsegdata1=getwholetable1(alllessionsfile1,opt="DEL")
#printsegdata(segdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_delseg.txt")
delsegdata2=getwholetable1(alllessionsfile2,opt="DEL")
#printsegdata(segdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_delseg.txt")
delsegdata3=getwholetable1(alllessionsfile3,opt="DEL")
#printsegdata(segdata3,output="/fh/scratch/delete30/dai_j/escc/gistic/escc_gistic_delseg.txt")
delolapseg12=countoverlap2segs(delsegdata1,delsegdata2)
#printsegdata(olapseg12$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_henan_gistic_deloverlapseg.txt")
printoverlap2segs(olapseg12)
delolapseg13=countoverlap2segs(delsegdata1,delsegdata3)
#printsegdata(olapseg13$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_escc_gistic_deloverlapseg.txt")
printoverlap2segs(olapseg13)
delolapseg23=countoverlap2segs(delsegdata2,delsegdata3)
#printsegdata(olapseg23$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/henan_escc_gistic_deloverlapseg.txt")
printoverlap2segs(olapseg23)
delolapseg123=countoverlap2segs(olapseg12$overlapseg,segdata3)
#printsegdata(olapseg123$overlapseg,output="/fh/scratch/delete30/dai_j/gistic/dulak_henan_escc_gistic_deloverlapseg.txt")
printoverlap2segs(olapseg123)
#unique segs:
uniqdelseg1=segonlyinseg1(delsegdata1,delsegdata2)
printsegdata(uniqdelseg1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_deldistinctseg.txt")
uniqdelseg2=segonlyinseg1(delsegdata2,delsegdata1)
printsegdata(delseg2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_deldistinctseg.txt")
delseg2=segonlyinseg1(delseg2,segdata3)
delseg3=segonlyinseg1(segdata3,segdata2)
printsegdata(delseg3,output="/fh/scratch/delete30/dai_j/escc/gistic/escc_gistic_deldistinctseg.txt")
delseg3=segonlyinseg1(delseg3,segdata1)

grid.newpage()
draw.triple.venn(area1 =olapseg12$width_segdata1, area2 = olapseg12$width_segdata2, area3 = olapseg13$width_segdata2, n12 = olapseg12$width_overlapseg, 
                 n23 = olapseg23$width_overlapseg, n13 = olapseg13$width_overlapseg, n123 = olapseg123$width_overlapseg, category = c("US-EA", "CH-EA", "CH-ESCC"), lty = "blank", 
                 fill = c("blue", "green", "red"),cex=rep(1.4,7),cat.cex = rep(1.5, 3),cat.col=c("blue", "green", "red"),euler.d=T)

#delres=overlapsegdata(segdata1,segdata2,segdata3)
# segdata1=addcolor4segdata(segdata=segdata1,segres=delres,col1="blue",col2="red")
# segdata1=segdata1[,-1]
# segdata2=addcolor4segdata(segdata=segdata2,segres=delres,col1="blue",col2="red")
# segdata2=segdata2[,-1]
# segdata3=addcolor4segdata(segdata=segdata3,segres=delres,col1="blue",col2="red")
# segdata3=segdata3[,-1]
segdata1=cbind.data.frame(segdata1,col=rep("blue",nrow(segdata1)))
segdata1$col=as.character(segdata1$col)
segdata1=segdata1[,-1]
segdata2=cbind.data.frame(segdata2,col=rep("blue",nrow(segdata2)))
segdata2$col=as.character(segdata2$col)
segdata2=segdata2[,-1]
segdata3=cbind.data.frame(segdata3,col=rep("blue",nrow(segdata3)))
segdata3$col=as.character(segdata3$col)
segdata3=segdata3[,-1]

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(segdata1, track.num=1)
RCircos.Heatmap.Plot2(segdata2, track.num=2)
RCircos.Heatmap.Plot2(segdata3, track.num=3)

#10samples from dulakdata
gisticdir1="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_10samples_0_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile1=paste0(gisticdir1,"/all_lesions.conf_95.txt")
gisticdir2="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_10samples_1_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile2=paste0(gisticdir2,"/all_lesions.conf_95.txt")
gisticdir3="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_10samples_2_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile3=paste0(gisticdir3,"/all_lesions.conf_95.txt")


opt="AMP"
opt="DEL"
segdata1=getwholetable1(alllessionsfile1,opt=opt)
segdata2=getwholetable1(alllessionsfile2,opt=opt)
segdata3=getwholetable1(alllessionsfile3,opt=opt)
olapseg12=countoverlap2segs(segdata1,segdata2)
printoverlap2segs(olapseg12)
olapseg13=countoverlap2segs(segdata1,segdata3)
printoverlap2segs(olapseg13)
olapseg23=countoverlap2segs(segdata2,segdata3)
printoverlap2segs(olapseg23)
olapseg123=countoverlap2segs(olapseg12$overlapseg,segdata3)
printoverlap2segs(olapseg123)
grid.newpage()
draw.triple.venn(area1 =olapseg12$width_segdata1, area2 = olapseg12$width_segdata2, area3 = olapseg13$width_segdata2, n12 = olapseg12$width_overlapseg, 
                 n23 = olapseg23$width_overlapseg, n13 = olapseg13$width_overlapseg, n123 = olapseg123$width_overlapseg, category = c("A", "B", "C"), lty = "blank", 
                 fill = c("blue", "green", "red"),cex=rep(1.4,7),cat.cex = rep(1.5, 3),cat.col=c("blue", "green", "red"))


segdata1=cbind.data.frame(segdata1,col=rep("red",nrow(segdata1)))
segdata1$col=as.character(segdata1$col)
segdata1=segdata1[,-1]
segdata2=cbind.data.frame(segdata2,col=rep("red",nrow(segdata2)))
segdata2$col=as.character(segdata2$col)
segdata2=segdata2[,-1]
segdata3=cbind.data.frame(segdata3,col=rep("red",nrow(segdata3)))
segdata3$col=as.character(segdata3$col)
segdata3=segdata3[,-1]


RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(segdata1, track.num=1)
RCircos.Heatmap.Plot2(segdata2, track.num=2)
RCircos.Heatmap.Plot2(segdata3, track.num=3)

segdata1=cbind.data.frame(segdata1,col=rep("blue",nrow(segdata1)))
segdata1$col=as.character(segdata1$col)
segdata1=segdata1[,-1]
segdata2=cbind.data.frame(segdata2,col=rep("blue",nrow(segdata2)))
segdata2$col=as.character(segdata2$col)
segdata2=segdata2[,-1]
segdata3=cbind.data.frame(segdata3,col=rep("blue",nrow(segdata3)))
segdata3$col=as.character(segdata3$col)
segdata3=segdata3[,-1]


RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(segdata1, track.num=1)
RCircos.Heatmap.Plot2(segdata2, track.num=2)
RCircos.Heatmap.Plot2(segdata3, track.num=3)

#2dataset
gisticdir="/fh/scratch/delete30/dai_j/henan/gistic2data/dulak_henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile=paste0(gisticdir,"/all_lesions.conf_95.txt")

ampsegdata=getwholetable1(alllessionsfile,opt="AMP")
delsegdata=getwholetable1(alllessionsfile,opt="DEL")
colnames(ampsegdata)[6:21]=paste0("US-EA",1:16)
colnames(ampsegdata)[22:ncol(ampsegdata)]=paste0("CH-EA",1:10)
colnames(delsegdata)=colnames(ampsegdata)
library("ComplexHeatmap")
library(circlize)
Heatmap(ampsegdata[,6:ncol(ampsegdata)],cluster_rows = T,name="amp event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1,2), labels = c("none", "gain","amplification")),
        row_names_gp = gpar(fontsize = 10)
        )
Heatmap(delsegdata[,6:ncol(delsegdata)],cluster_rows = T,name="del event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1), labels = c("none", "loss")),
        row_names_gp = gpar(fontsize = 10)
)


#compare two datasets with one dataset

ampsegdata1d=cbind.data.frame(ampsegdata1,col=rep("red",nrow(ampsegdata1)),stringsAsFactors = FALSE)
ampsegdata1d=ampsegdata1d[,-1]
ampsegdata2d=cbind.data.frame(ampsegdata2,col=rep("red",nrow(ampsegdata2)),stringsAsFactors = FALSE)
ampsegdata2d=ampsegdata2d[,-1]
ampsegdatad=cbind.data.frame(ampsegdata,col=rep("red",nrow(ampsegdata)),stringsAsFactors=FALSE)
ampsegdatad=ampsegdatad[,-1]

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(ampsegdata1d, track.num=1)
RCircos.Heatmap.Plot2(ampsegdata2d, track.num=2)
#RCircos.Heatmap.Plot2(ampsegdata_d, track.num=3)


segdata1d=cbind.data.frame(segdata1,col=rep("blue",nrow(segdata1)),stringsAsFactors = FALSE)
segdata1d=segdata1d[,-1]
segdata2d=cbind.data.frame(segdata2,col=rep("blue",nrow(segdata2)),stringsAsFactors = FALSE)
segdata2d=segdata2d[,-1]
segdatad=cbind.data.frame(segdata,col=rep("blue",nrow(segdata)),stringsAsFactors=FALSE)
segdatad=segdatad[,-1]

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(segdata1d, track.num=1)
RCircos.Heatmap.Plot2(segdata2d, track.num=2)
#RCircos.Heatmap.Plot2(segdata_d, track.num=3)

#only work on combined data
get_proportion=function(segdata=ampsegdata,idx1=9,idx2=24,opt1="peak")
{
  res1=res2=data.frame(matrix(NA,nrow=nrow(segdata),ncol=5))
  colnames(res1)=colnames(res2)=c("cytoband","chr","start","end","proportion")
  res1$cytoband=res2$cytoband=segdata$Descriptor
  
  for (i in 1:nrow(segdata))
  {
    if (opt1=="peak")
    {
      tmp=unlist(strsplit(alllessiontable$Peak.Limits[i],":",fixed=T))
    }else
    {
      tmp=unlist(strsplit(alllessiontable$Region.Limits[i],":",fixed=T))
    }

    res1$chr[i]=res2$chr[i]=tmp[1]
    tmp=unlist(strsplit(tmp[2],"(",fixed=T))
    tmp=unlist(strsplit(tmp[1],"-",fixed=T))
    res1$start[i]=res2$start[i]=as.integer(tmp[1])
    res1$end[i]=res2$end[i]=as.integer(tmp[2])
    res1$proportion[i]=sum(segdata[i,idx1:idx2]>0)/(idx2-idx1+1)
    res2$proportion[i]=sum(segdata[i,(idx2+1):ncol(segdata)]>0)/(ncol(segdata)-idx2)
  }
  return(result=list(res1=res1,res2=res2))
}
uniq_uniqampseg1=segonlyinseg1(uniqampseg1,ampsegdata)
uniq_uniqampseg2=segonlyinseg1(uniqampseg2,ampsegdata)
uniq_ampseg=segonlyinseg1(ampsegdata,rbind(ampsegdata1[,1:4],ampsegdata2[,1:4]))

uniq_uniqdelseg1=segonlyinseg1(uniqdelseg1,delsegdata)
uniq_uniqdelseg2=segonlyinseg1(uniqdelseg2,delsegdata)
uniq_delseg=segonlyinseg1(delsegdata,rbind(delsegdata1[,1:4],delsegdata2[,1:4]))

uniq_uniq_uniqdelseg1=segonlyinseg1(uniqdelseg1,uniq_uniqdelseg1)
idx=which(!uniq_uniqdelseg1 %in% uniq_uniqdelseg1$cytoband)

Heatmap(delsegdata[,6:ncol(delsegdata)],cluster_rows = T,name="del event",show_column_names=T,
        cluster_columns = F,show_column_dend=F,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1), labels = c("none", "loss")),
        row_names_gp = gpar(fontsize = 10)
)
