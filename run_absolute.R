#! /usr/bin/env Rscript
setwd("/fh/scratch/delete30/dai_j/henan/absolute")
for (numpair in 1:10)
{
  cmd=paste0("sbatch ./absolute1.R ",numpair)
  system(cmd)
  cmd="sleep 1"
  system(cmd)
}
