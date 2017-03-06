#!/usr/bin/env bash

data="henan_normals"
echo $data
sleep 10
if [[ $data == "henan" ]]
then
  #chr=16
  freecdir=/fh/scratch/delete30/dai_j/henan/freec
  tumors=()
  tumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
  for ((i=0;i<${#tumors[@]};i++))
  do
    preface=$freecdir/${tumors[$i]}/ploid2degree3force0/${tumors[$i]}
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/copynumber/cnv_significance.R $preface
    sleep 1
  done
fi

if [[ $data == "henan_normals" ]]
then
  freecdir=/fh/scratch/delete30/dai_j/henan/freec
  tumors=()
  tumors=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
  for ((i=0;i<${#tumors[@]};i++))
  do
    preface=$freecdir/${tumors[$i]}/ploid2degree3force0/${tumors[$i]}
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/copynumber/cnv_significance.R $preface
    sleep 1
  done
fi
 
if [[ $data == "escc" ]]
then
  #chr=16
  freecdir=/fh/scratch/delete30/dai_j/escc/freec
  ids=( {{1..6},{8..18}} )
  tumors=()
  for ((i=0;i<${#ids[@]};i++))
  do
    tumors[$i]=T${ids[$i]}
  done
  for ((i=0;i<${#tumors[@]};i++))
  do
    preface=$freecdir/${tumors[$i]}/ploid2degree3force0/${tumors[$i]}
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/copynumber/cnv_significance.R $preface
    sleep 1
  done
fi

if [[ $data == "dulak" ]]
then
  #chr=16
  freecdir=/fh/scratch/delete30/dai_j/freec
  tumors=()
  tumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
  for ((i=0;i<${#tumors[@]};i++))
  do
    preface=$freecdir/${tumors[$i]}/ploid2degree3force0/${tumors[$i]}
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/copynumber/cnv_significance.R $preface
    sleep 1
  done
fi
