#! /usr/bin/env bash
#SBATCH -t 5-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org
#SBATCH --mem 30000
#the output is in $outdir/$expname

#the chromosome id
#chr=${1?"the chromsomeid"}
numpair=${1?"the samplepairid"}
ploidy=${2?"ploidy"}
degree=${3?"degree"}

reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta
#data=dulak
#data=henan_normals
if [[ $data == "dulak" ]]
then
  normal=(SRR1002719 SRR999433 SRR999687 SRR1000378 SRR1002786 SRR999559 SRR1001730 SRR10018461 SRR1002703 SRR999599 SRR1002792 SRR1001839 SRR9994281 SRR1002710 SRR9995631 SRR1021476)
  tumor=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
#  	       0          1         2          3         4           5         6         7          8         9          10         11        12         13       14          15   
  fromdir=/fh/scratch/delete30/dai_j/freec
  bamdir=/fh/scratch/delete30/dai_j/dulak
  ext=dedup.realigned.recal.bam
fi


if [[ $data == "henan" ]]
then
  normal=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
  tumor=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
#  	  0  1   2   3   4   5   6   7   8   9   
  fromdir=/fh/scratch/delete30/dai_j/henan/freec
  bamdir=/fh/scratch/delete30/dai_j/henan
  ext=merged.deduprealigned.bam
fi

if [[ $data == "henan_normals" ]]
then
  normal=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
  tumor=(14A 26A 4A 42A 34A 12A 38A 18A 30A 16A)
#  	  0  1   2   3   4   5   6   7   8   9   
  fromdir=/fh/scratch/delete30/dai_j/henan/freec
  bamdir=/fh/scratch/delete30/dai_j/henan
  ext=merged.deduprealigned.bam
fi

if [[ $data == "escc" ]]
then
  normal=(N1 N2 N3 N4 N5 N6 N8 N9 N10 N11 N12 N13 N14 N15 N16 N17 N18)
  tumor=(T1 T2 T3 T4 T5 T6 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18)
  fromdir=/fh/scratch/delete30/dai_j/escc/freec
  bamdir=/fh/scratch/delete30/dai_j/escc
  ext=merged.deduprealigned.bam
fi
freecdir=/fh/fast/dai_j/CancerGenomics/Tools/FREEC
normalsample=${normal[$numpair]}
tumorsample=${tumor[$numpair]}

echo $normalsample
echo $tumorsample

outputdir=$fromdir/$tumorsample/ploid${ploidy}degree${degree}force0
if [[ ! -d $fromdir/$tumorsample/ ]];then mkdir $fromdir/$tumorsample/;fi
if [[ ! -d $outputdir ]];then mkdir $outputdir;fi

format=pileup
orientation=FR
normalpileup=$fromdir/$normalsample.pileup.gz
tumorpileup=$fromdir/$tumorsample.pileup.gz

normalbam=$bamdir/$normalsample.$ext
tumorbam=$bamdir/$tumorsample.$ext
if [[ ! -s $normalpileup ]]
then 
  normalpileup=$fromdir/$normalsample.pileup
  echo "start to generate pileup file..." 
  samtools mpileup -q 1 -f $reference $normalbam >$normalpileup
  gzip $normalpileup
fi
if [[ ! -s $tumorpileup ]] 
then 
  tumorpileup=$fromdir/$tumorsample.pileup
  samtools mpileup -q 1 -f $reference $tumorbam >$tumorpileup
  gzip $tumorpileup
fi

normalpileup=$fromdir/$normalsample.pileup.gz
tumorpileup=$fromdir/$tumorsample.pileup.gz

cat $freecdir/freecwgsconfig1.txt | sed -e "s|myOutputDirectory|$outputdir|g" -e "s|TUMORARACHNE|$tumorpileup|g" -e "s|NORMALARACHNE|$normalpileup|g" -e "s|myFormat|$format|g" -e"s|myOrientation|$orientation|g" -e "s|myPloidy|$ploidy|g" -e "s|mydegree|$degree|g" - >$outputdir/${tumorsample}.config.$format.txt

$freecdir/freec -conf $outputdir/${tumorsample}.config.$format.txt 2>&1 | tee $outputdir/$tumorsample.log
