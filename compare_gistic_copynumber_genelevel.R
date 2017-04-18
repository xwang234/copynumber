#/usr/bin/env Rscript

gisticdir1="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
file1=paste0(gisticdir1,"/all_thresholded.by_genes.txt")
wgstumors1=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
             "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868")
wgstumors2=paste0("X",c(3,11,13,15,17,25,29,33,37,41),"A")
gisticdir2="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
file2=paste0(gisticdir2,"/all_thresholded.by_genes.txt")
wgstumors3=paste0("T",c(1:6,8:18))
wgstumors3=paste0("T",c(1:4,6,8:18))
gisticdir3="/fh/scratch/delete30/dai_j/escc/gistic/escc_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
file3=paste0(gisticdir3,"/all_thresholded.by_genes.txt")

compute_pvalues_gistic_gene=function(file1,file2)
{
  data1=read.table(file1,header=T,sep="\t",stringsAsFactors=F)
  data2=read.table(file2,header=T,sep="\t",stringsAsFactors=F)
  data12=merge(data1,data2,by="Gene.Symbol")
  data1=data12[,1:ncol(data1)]
  data2[,1:3]=data12[,1:3]
  data2[,4:ncol(data2)]=data12[,(ncol(data1)+3):ncol(data12)]
  res=data.frame(matrix(NA,ncol=9,nrow=nrow(data1)))
  res[,1:3]=data1[,1:3]
  colnames(res)=c(colnames(data1)[1:3],"num_amp1","num_amp2","p_amp","num_del1","num_del2","p_del")
  for (i in 1:nrow(data1))
  {
    numamp1=sum(data1[i,4:ncol(data1)]>=1)
    numamp2=sum(data2[i,4:ncol(data2)]>=1)
    mat=matrix(c(numamp1,ncol(data1)-3-numamp1,numamp2,ncol(data2)-3-numamp2),nrow=2,byrow = T)
    res$num_amp1[i]=numamp1
    res$num_amp2[i]=numamp2
    res$p_amp[i]=fisher.test(mat)$p.value
    numdel1=sum(data1[i,4:ncol(data1)] <= -1)
    numdel2=sum(data2[i,4:ncol(data2)] <= -1)
    mat=matrix(c(numdel1,ncol(data1)-3-numdel1,numdel2,ncol(data2)-3-numdel2),nrow=2,byrow = T)
    res$num_del1[i]=numdel1
    res$num_del2[i]=numdel2
    res$p_del[i]=fisher.test(mat)$p.value
  }
  return(res)
}

res_dulak_henan=compute_pvalues_gistic_gene(file1,file2)
res_dulak_escc=compute_pvalues_gistic_gene(file1,file3)
res_henan_escc=compute_pvalues_gistic_gene(file2,file3)

