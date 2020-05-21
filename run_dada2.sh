#!/usr/bin/bash
#$ -cwd
#$ -pe smp 10
#$ -l mem_free=10G
#$ -l h_rt=100:00:00
#$ -m eab
#$ -M kathryn.mccauley@ucsf.edu

module load CBI r/3.6.2

Rscript --version

Rscript $1 $NSLOTS

qstat -j $JOB_ID
