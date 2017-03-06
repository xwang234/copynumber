#!/usr/bin/env Rscript

getwholetable=function(alllessiontable,wgstumors)
{
  alllessiontable$Descriptor=gsub(" ","",alllessiontable$Descriptor,fixed=T)
  uniq_cytoband=unique(alllessiontable$Descriptor)
  cytobandtable=data.frame(matrix(" ",nrow=length(uniq_cytoband),ncol=length(wgstumors)))
  rownames(cytobandtable)=uniq_cytoband
  colnames(cytobandtable)=wgstumors
  for (i in 1:ncol(cytobandtable)) cytobandtable[,i]=as.character(cytobandtable[,i])
  c_gain=0.1
  c_amp=0.9
  c_loss=-0.1
  c_del=-1.2
  for (i in 1:nrow(cytobandtable))
  {
    cytoband=rownames(cytobandtable)[i]
    idx=which(alllessiontable$Descriptor==cytoband)
    if (length(idx)>1) #multiple peaks 
    {
      idx=idx[which.min(alllessiontable$q.values[idx])]
    }
    for (j in 10:(9+length(wgstumors)))
    {
      if (alllessiontable[idx,j]<c_del)
      {
        cytobandtable[i,j-9]="Deletion"
      }else
      {
        if (alllessiontable[idx,j]<c_loss)
        {
          cytobandtable[i,j-9]="Loss"
        }else
        {
          if (alllessiontable[idx,j]>c_amp)
          {
            cytobandtable[i,j-9]="Amplification"
          }else
          {
            if (alllessiontable[idx,j]>c_gain)
            {
              cytobandtable[i,j-9]="Gain"
            }
          }
        }
      }
    }
  }
  return(cytobandtable)
}

getampdeltable=function(alllessiontable,wgstumors)
{
  amplessiontable=alllessiontable[grepl("Amplification",alllessiontable$Unique.Name),]
  ampcytotable=getwholetable(amplessiontable,wgstumors)
  for (i in 1:nrow(ampcytotable))
  {
    for (j in 1:ncol(ampcytotable))
    {
      ampcytotable[i,j]=gsub("Loss"," ",ampcytotable[i,j])
      ampcytotable[i,j]=gsub("Deletion"," ",ampcytotable[i,j])
    }
  }
  dellessiontable=alllessiontable[! grepl("Amplification",alllessiontable$Unique.Name),]
  delcytotable=getwholetable(dellessiontable,wgstumors)
  for (i in 1:nrow(delcytotable))
  {
    for (j in 1:ncol(delcytotable))
    {
      delcytotable[i,j]=gsub("Gain"," ",delcytotable[i,j])
      delcytotable[i,j]=gsub("Amplification"," ",delcytotable[i,j])
    }
  }
  result=list(ampcytotable=ampcytotable,delcytotable=delcytotable)
  return(result)
}


library(ComplexHeatmap)
wgstumors=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
            "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868","SRR1001635")
#gisticdir="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/nwgs_rx1_conf0.95_armpeel0_brlen0.98_broad1"
gisticdir="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
name="US-EA"
allsegfile="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv.combinedfreecseg.txt"


wgstumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
gisticdir="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
name="CH-EA"
allsegfile="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv.combinedfreecseg.txt"

wgstumors=paste0("T",c(1:6,8:18))
wgstumors=paste0("T",c(1:4,6,8:18))
gisticdir="/fh/scratch/delete30/dai_j/escc/gistic/escc_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
name="CH-ESCC"
allsegfile="/fh/scratch/delete30/dai_j/escc/gistic/escc_ploid2degree3force0_cnv.combinedfreecseg.txt"


alllessionsfile=paste0(gisticdir,"/all_lesions.conf_95.txt")
alllessiontable=read.table(alllessionsfile,header=T,sep="\t",stringsAsFactors=F)
alllessiontable1=alllessiontable[1:(nrow(alllessiontable)/2),]
idxkeep=which(rowSums(alllessiontable1[,10:(10+length(wgstumors)-1)],na.rm = T)>1)
alllessiontable=alllessiontable[(nrow(alllessiontable)/2+1):nrow(alllessiontable),]
alllessiontable=alllessiontable[idxkeep,]
#cytobandtable=getwholetable(alllessiontable)
twotables=getampdeltable(alllessiontable,wgstumors)

#count number of alterations
numdelalt=numampalt=numalt=rep(0,length(wgstumors))
allseg=read.table(file=allsegfile,header=F,sep="\t")
for (i in 1:length(wgstumors))
{
  tmptable=allseg[allseg[,1]==wgstumors[i],]
  numalt[i]=sum(tmptable[,6] < -0.1 | tmptable[,6] > 0.1)
  numdelalt[i]=sum(tmptable[,6] < -0.1)
  numampalt[i]=sum(tmptable[,6] > 0.1)
}

copynumber_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Amplification=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Gain=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "orange", col = NA))
  },
  Loss=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "skyblue", col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA))
  }
)

col = c("Amplification" = "red", "Gain"="orange",
        "Loss"="skyblue", "Deletion" = "blue")


# oncoPrint(as.matrix(cytobandtable), get_type = function(x) strsplit(x, ";")[[1]],
#           row_order = NULL,column_order = NULL,
#           remove_empty_columns = TRUE,
#           alter_fun = copynumber_fun, col = col, 
#           column_title = name,
#           column_title_gp = gpar(fontsize = 22),
#           show_column_names = FALSE,
#           show_pct = FALSE,
#           axis_gp = gpar(fontsize = 16),# size of axis
#           row_names_gp = gpar(fontsize = 16),  # set size for row names
#           pct_gp = gpar(fontsize = 16), # set size for percentage labels
#           row_barplot_width = unit(4, "cm"), #size barplot
#           heatmap_legend_param = list(title = "Copynumber", at = c("Amplification", "Gain", "Loss","Deletion"), 
#                                       labels = c("Amplification", "Gain", "Loss","Deletion")))

ha = HeatmapAnnotation(total_alterations= anno_barplot(numampalt,axis=T,axis_gp = gpar(fontsize = 12),axis_side="right",border=F),
                       show_annotation_name = T,annotation_name_offset = unit(2, "cm"),gap = unit(3, "mm"))

oncoPrint(as.matrix(twotables$ampcytotable), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = TRUE,
          alter_fun = copynumber_fun, col = col, 
          column_title =name,
          column_title_gp = gpar(fontsize = 18),
          show_column_names = FALSE,
          show_pct = FALSE,
          axis_gp = gpar(fontsize = 16),# size of axis
          row_names_gp = gpar(fontsize = 16),  # set size for row names
          pct_gp = gpar(fontsize = 16), # set size for percentage labels
          row_barplot_width = unit(4, "cm"), #size barplot
          bottom_annotation=ha,
          bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "Copynumber", at = c("Amplification", "Gain", "Loss","Deletion"), 
                                      labels = c("Amplification", "Gain", "Loss","Deletion")))

ha = HeatmapAnnotation(total_alterations= anno_barplot(numdelalt,axis=T,axis_gp = gpar(fontsize = 12),axis_side="right",border=F),
                       show_annotation_name = T,annotation_name_offset = unit(2, "cm"),gap = unit(3, "mm"))

oncoPrint(as.matrix(twotables$delcytotable), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = TRUE,
          alter_fun = copynumber_fun, col = col, 
          column_title = name,
          column_title_gp = gpar(fontsize = 18),
          show_column_names = FALSE,
          show_pct = FALSE,
          axis_gp = gpar(fontsize = 16),# size of axis
          row_names_gp = gpar(fontsize = 10),  # set size for row names
          pct_gp = gpar(fontsize = 16), # set size for percentage labels
          row_barplot_width = unit(4, "cm"), #size barplot
          bottom_annotation=ha,
          bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "Copynumber", at = c("Amplification", "Gain", "Loss","Deletion"), 
                                      labels = c("Amplification", "Gain", "Loss","Deletion")))