#! /usr/bin/env bash
#SBATCH -t 0-8
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org
#SBATCH --mem 30000

#run in the freec folder

normalbam=${1?"normalbam"}
tumorbam=${2?"tumorbam"}
chr=${3?"chrid,0-23"}
normalname=$(basename $normalbam)
normalname=${normalname%%.*}
tumorname=$(basename $tumorbam)
tumorname=${tumorname%%.*}

echo $normalbam
echo $tumorbam

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
chr=chr${chrs[$chr]}

#output folder:
if [ ! -d $tumorname ];then mkdir $tumorname;fi
if [ ! -d $tumorname/$chr ];then mkdir $tumorname/$chr;fi
outputdir=$tumorname/$chr

freecdir=/fh/fast/dai_j/CancerGenomics/Tools/FREEC/FREEC-10.2/src
reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta

format=pileup
orientation=0 #for sorted
degree=3

echo "create mpileup..."
samtools mpileup -B -q 1 -f $reference -r $chr $normalbam >${SCRATCH}/${normalname}_$chr.mpileup
samtools mpileup -B -q 1 -f $reference -r $chr $tumorbam >${SCRATCH}/${tumorname}_$chr.mpileup

gzip ${SCRATCH}/${normalname}_$chr.mpileup
gzip ${SCRATCH}/${tumorname}_$chr.mpileup


cat $freecdir/freecwgsconfig4.txt | sed -e "s|myOutputDirectory|$outputdir|g" -e "s|NORMALARACHNE|${SCRATCH}/${normalname}_$chr.mpileup.gz|g" -e "s|TUMORARACHNE|${SCRATCH}/${tumorname}_$chr.mpileup.gz|g" -e "s|myFormat|$format|g" -e"s|myOrientation|$orientation|g" -e "s|mydegree|$degree|g" - >${tumorname}_$chr.config.$format.txt
$freecdir/freec -conf ${tumorname}_$chr.config.$format.txt 2>&1 | tee ${tumorname}_$chr.log

