#!/usr/bin/bash
#$ -cwd
#$ -pe smp 18
#$ -l mem_free=1T
#$ -l h_rt=100:00:00
#$ -m eab
#$ -M kathryn.mccauley@ucsf.edu

module load CBI r

Rscript --version

Rscript $1 $NSLOTS

qstat -j $JOB_ID
