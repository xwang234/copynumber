#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org
normalbam=${1?"normalbamfile"}
tumorbam=${2?"tumorbamfile"}
rundir=${3?"rundir"}

echo $tumorbam
echo $SLURM_CPUS_ON_NODE

if [ ! -d $rundir ]; then mkdir $rundir; fi
#generate configure file
/fh/fast/dai_j/CancerGenomics/Tools/manta-1.0.3.centos5_x86_64/bin/configManta.py \
--normalBam $normalbam \
--tumorBam $tumorbam \
--referenceFasta /fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta \
--runDir $rundir

sleep 1
$rundir/runWorkflow.py -m local -j $SLURM_CPUS_ON_NODE

