#! /usr/bin/env bash

data="escc"
if [[ $data=="escc" ]]
then

  bamdir=/fh/scratch/delete30/dai_j/escc
  ids={{1..6},{8..18}}
  normals=()
  tumors=()
  for id in ${ids[@]}
  do
    $normals[$i]=N$id
    $tumors[$i]=T$id
  done
  for ((i=0;i<${#normals[@]};i++))
  do
     normalbam=$bamdir/${normals[$i]}.merged.deduprealigned.bam
     tumorbam=$bamdir/${tumors[$i]}.merged.deduprealigned.bam
    for chr in {0..23}
    do
        echo sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/freec_chr.sh $normalbam $tumorbam $chr
        #sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/freec_chr.sh $normalbam $tumorbam $chr
        sleep 1
    done
  done
fi
