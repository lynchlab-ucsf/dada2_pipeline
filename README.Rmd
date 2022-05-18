## Scripts for Running DADA2 on NextSeq and MiSeq

These scripts are used to run DADA2 based on various types NextSeq or MiSeq runs. Based on the protocol used, each has various nuances requiring distinct scripts. Defunct scripts are included in a local `archive` dir, not tracked here.

*DADA2_NextSeqProcessing_30Oct.R* - The version of the script used for our original NextSeq runs where we used QIIMEv1 to generate the forward and reverse files in to distinct submission directories. Still used sometimes, so this script remains active.

*DADA2_NextSeq_SampSheet.R* - Using the Illumina SampleSheet format creates differently-structured names, thus requiring modifications to the original script.

*DADA2_MiSeq_16s.R* - Much like the original NextSeq, a script based on QIIMEv1's initial output to generate 16S rRNA SV table. Similar to 2x150, except not needing a sampsheet.

*DADA2_MiSeq_ITS2.R* - Script specifically for ITS2 sequencing, with links to the ITS2 UNITE database (for taxonomic assignment) and other important checks specific to ITS2.

*DADA2_MiSeq_SampSheet_V1V3.R* - Script to process a MiSeq run with a much longer read region (here V1-V3, though could probably be applied to anything up to 300 bp.

*DADA2_MiSeq_2x150_SampSheet.R* - As noted above, a script that combines a 2x150 MiSeq run with the use of a Sample Sheet csv file.

## Additional Scripts

This repo also includes scripts for processing data after running through DADA2, including:

*procrustes.R* - Command line interface (CLI) that allows you to compare two phyloseq objects with procrustes. Basic use is : `Rscript procrustes.R myphy1.rds myphy2.rds`. This script is helpful after multiply-rarefying to determine if different rarefying depths result in different data structures.

*alpha_rarefaction.R* - An early script to run alpha rarefaction. It's a little rough, but gives you a sense of the samples you retain at each sample depth, and the change in diversity across different depths.

*Multiply_Rarefy_phy.R* - The multiply-rarefy script modified for use with phyloseq objects. Uses the CLI, with basic use being `Rscript Multiply_Rarefy_phy.R inphy.rds 25000 outphy.rds` to rarefy a dataset to 25k reads per sample.

*DADA2_PostProcessing_phy.R* - Performs "routine" things after running dada2. This includes tree generation and low-abundance taxon filtering.

*NegControlRemoval_phy.R* - Custom negative control signal removal script that works on phyloseq objects from the CLI. Has decent documentation that prints out automatically when the script is run.

*NegControlRemoval_func.R* - Modification of the phyloseq version except it can be sourced into your R environment and used in your active R analysis script. Especially helpful if you need to check parameters and really understand the negative control signal (ie, needing to run filtering on a plate-by-plate basis, etc).
