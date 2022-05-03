#!/usr/bin/Rscript

## DADA2 Implementation 
## Author: Katie McCauley, based on the tutorial by Ben Callahan
## Most-recent update by Preston Tasoff January 6, 2021

## This script's main modification is to truncate the 300bp MiSeq read length so that it is shorter than the V4 region of 253 bp to allow for sufficient overlap without read-through into the amplicon. This is the main change from the 16S pipeline from the NextSeq pipeline


# 05/13/19: Add in code at the beginning to install first install pacman (if missing) and then to install other difficult packages without user input.
# 02/14/20: Modifying to run fully on each run instead of merging several runs. Should, theoretically, output run-specific OTU tables that can then be merged as needed.
# 05/21/20: Modifying to run after NextSeq Processing script
# 10/30/20: Modified to include generation of phyloseq object after processing, along with sequences. Removed tictoc package.
# 01/06/21: Preston added truncation to avoid the read-through of the 16S amplicon.

rm(list=ls())
args <- commandArgs(TRUE)
# This now checks for basic/difficult packages that need to be installed and installs them if not already done, but it isn't working for the most-recent update.
#if(!"pacman" %in% installed.packages()) install.packages("pacman", repos="https://cran.revolutionanalytics.com/")
#if(!"devtools" %in% installed.packages()) {
#install.packages("devtools", repos="https://cran.revolutionanalytics.com/")
#}
#if(!"dada2" %in% installed.packages()) {
#devtools::install_github("benjjneb/dada2", ref="v1.12") # change the ref argument to get other versions
#}
#if(!"DECIPHER" %in% installed.packages()) {
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DECIPHER")
#}

.libPaths()

pacman::p_load(dada2, parallel, purrr, DECIPHER, dplyr, tibble, phyloseq, ShortRead)

.libPaths()
## Add directory locations so that DADA2 code can be kept in a clean space that is distinct from all output files
## These directory locations shouldn't need to change, but ar available if needed
out_path <- "/wynton/group/lynch/MiSeq_data/MiSeq_Processed/"
dada2_path <- "/wynton/group/lynch/kmccauley/dada2_files/"

#source(paste0(dada2_path, "DADA2_modfunc.R")) # I am modifying some of the DADA2 functions. To see specific details, view the "DADA2_modfunc.R" script and search for "(KM)". I have noted my modifications and my reasoning for the change(s). Keep this in the seq_path directory.
threads <- as.numeric(args[1]) #Number of cores you want to use throughout the process. Modified here to use the number of cores specified in the bash script.
#You can also designate a "save path" (so that you can save your final files in a separate place from your sequence files)
threshold <- 50 # Drop a sample if it contains less than this threshold of reads (shouldn't be lower than 20).

## For the nextseq Processing, I want to be able to run without a CSV file -- basically just pull from whatever forward and reverse files exist in the submission directory
csv_file <- NULL

mapping_file <- as.character(args[3])

###################################################
#####       START DADA2 CODE            ###########
###################################################

 ns_dir <- as.character(args[2])
 save_path <- paste0(out_path, ns_dir, "/dada2_output")
map_file <- read.table(paste0(out_path, ns_dir, "/",mapping_file), check.names=F, sep="\t", comment="", header=TRUE)
rownames(map_file) <- map_file[,1]

print("Save path:")
save_path

## Developing the ability for a skip pattern to be implemented:
if(file.exists(paste0(save_path, "/dada2_result_object.RData"))) { ## If this has already been generated, then it gets loaded in.
   load(paste0(save_path, "/dada2_result_object.RData"))
}

if(!exists("initial_step")) {
 fnFs <- list()
 fnRs <- list()
 filtFs <- list()
 filtRs <- list()
 i <- 1
 dir.create(save_path, showWarnings=FALSE)
 fnFs[[i]] <- list.files(paste0(out_path, ns_dir,"/submission/R1"), pattern="fastq", full.names=TRUE)
 fnRs[[i]] <- list.files(paste0(out_path, ns_dir,"/submission/R2"), pattern="fastq", full.names=TRUE)
 filt_path <- file.path(save_path, "filtered")
 sample.names <- gsub("_R1.fastq.gz","",basename(fnFs[[i]]))
 filtFs[[i]] <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
 filtRs[[i]] <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))
initial_step <- 1
}

i <- 1
if(!exists("filterTrimComplete")) {
print("Begin the Filter and Trim Step")
length(filtFs[[i]])
length(fnFs[[i]])
out <- filterAndTrim(fwd=unlist(fnFs[[i]]), filt=unlist(filtFs[[i]]), rev=unlist(fnRs[[i]]), filt.rev=unlist(filtRs[[i]]), 
              maxEE=2, truncQ=2, rm.phix=TRUE, minLen=50,
              compress=TRUE, multithread=threads, verbose=TRUE, matchIDs = TRUE)

prop.out <- out[,"reads.out"]/out[,"reads.in"]
out <- cbind(out, prop.out)
print(out)

keep <- out[,"reads.out"] > threshold

keep.names <- gsub("R1.fastq.gz", "", names(keep[keep])) ## Removing the underscore to denote the end of the sample name in those that get dropped (for instance if S3 gets dropped, but there's also an S30, it would get dropped previously) This should mitigate that.
if(length(keep.names)<nrow(out)) {

filtFs[[i]] <- filtFs[[i]][grepl(paste(keep.names, collapse="|"), filtFs[[i]])]
filtRs[[i]] <- filtRs[[i]][grepl(paste(keep.names, collapse="|"), filtRs[[i]])]

print(paste("Samples were dropped due to having fewer than", threshold, "reads per sample:"))
print(paste(rownames(out)[!keep], collapse=","))
}
save.image(paste0(save_path, "/dada2_result_object.RData"))
filterTrimComplete <- 1
} else {
print("Filter and Trim Already Done")
}


save.image(paste0(save_path, "/dada2_result_object.RData"))
if(!exists("error_done")) {
errF <- list()
errR <- list()
errF[[i]] <- learnErrors(filtFs[[i]], multithread=threads, randomize=TRUE, nbases=1e8)
errR[[i]] <- learnErrors(filtRs[[i]], multithread=threads, randomize=TRUE, nbases=1e8)
error_done <- 1
} else {
print("Error Checking Complete")
}
save.image(paste0(save_path, "/dada2_result_object.RData"))
## Run DADA2
if(!exists("dada_done")) {
print("DADA2 Step")
dadaFs <- list()
dadaRs <- list()
dadaFs[[i]] <- dada(filtFs[[i]], err=errF[[i]], pool=FALSE, multithread=threads)
dadaRs[[i]] <- dada(filtRs[[i]], err=errR[[i]], pool=FALSE, multithread=threads)
dada_done <- 1
} else {
print("DADA denoising already complete")
}

 # This will be the end of the run-specific processing for now.
save.image(paste0(save_path, "/dada2_result_object.RData"))

## NOW Merge Paired Ends
if(!exists("merge_reads")) {
print("Merge Paired Reads")
library(purrr) #bringing in this package because it "undoes" the list-of-lists that I had, which are no longer needed because I have already taken care of the run-specific processing
mergers <- mergePairs(flatten(dadaFs), unlist(filtFs), flatten(dadaRs), unlist(filtRs), minOverlap=25, verbose=TRUE)

#I think this is useful information to print, so I will keep it here.
head(mergers[[1]])
str(mergers[[1]])
dim(mergers[[2]])
merge_reads <- 1
save.image(paste0(save_path, "/dada2_result_object.RData"))

} else {
print("Read Merging already Complete")
}

#Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#A data check..
print("Initial Sequence Lengths Pre-Chimera Removal")
table(nchar(getSequences(seqtab)))
# Distribution of sequence reads:
seq.length <- sort(table(nchar(getSequences(seqtab))))
## Remove Chimeras:

if(!exists("chimera_check")) {
seqtab.nochim.init <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)
dim(seqtab.nochim.init)
print("Proportion of non-chimeric reads")
sum(seqtab.nochim.init)/sum(seqtab)
chimera_check <- 1
} else {print("Chimera Check already complete")}

print("Initial Sequence Lengths")
table(nchar(getSequences(seqtab.nochim.init))) # Moving this to *after* the bimera removal, since it's possible that the odd-sized reads are considered bimeras.

seqtab.nochim <- seqtab.nochim.init[, nchar(getSequences(seqtab.nochim.init)) > 100]

write.table(x=getSequences(seqtab.nochim), file=paste0(save_path,"/", ns_dir,"_seqs.txt"), sep="\t", col.names=F, row.names=F, quote=F)
write.table(x=t(seqtab.nochim), file=paste0(save_path,"/",ns_dir,"_otutable.txt"), sep="\t", col.names=T, row.names=T, quote=F)

print("New Sequence Length Distribution")
table(nchar(getSequences(seqtab.nochim)))

save.image(paste0(save_path, "/dada2_result_object.RData"))

getN <- function(x) sum(getUniques(x))
sumdadaFs <- as.data.frame(sapply(flatten(dadaFs), getN)) %>%
	rownames_to_column(var="sampleid") %>%
	mutate(sampleid = gsub( "_R1_filt.fastq.gz","",sampleid)) %>%
	rename('denoisedF' = 'sapply(flatten(dadaFs), getN)')
sumdadaRs <- as.data.frame(sapply(flatten(dadaRs), getN)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R2_filt.fastq.gz","",sampleid)) %>%
        rename('denoisedR' = 'sapply(flatten(dadaRs), getN)')
sum.merge <- as.data.frame(sapply(mergers, getN)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R1_filt.fastq.gz","",sampleid)) %>%
        rename('merged' = 'sapply(mergers, getN)')
sum.seqtab <- as.data.frame(rowSums(seqtab.nochim)) %>%
        rownames_to_column(var="sampleid") %>%
        mutate(sampleid = gsub("_R1_filt.fastq.gz","",sampleid)) %>%
        rename('nonchim' = 'rowSums(seqtab.nochim)')
sum.out <- as.data.frame(out) %>%
	rownames_to_column(var="sampleid") %>%
	mutate(sampleid = gsub("_R1.fastq.gz","", sampleid))
track <- sum.out %>%
	full_join(sumdadaFs) %>%
	full_join(sumdadaRs) %>%
	full_join(sum.merge) %>%
	full_join(sum.seqtab)

write.csv(track, paste0(save_path, "/TrackedReadsThruDADA2.csv"))
#rm(list=c("dadaFs","dadaRs","mergers"))

## Taxonomy

#Need to figure out where Silva will end up going. Right now, I have copied it into the directory for this run, but that's not sustainable.
#Discussion of whether this should be DECIPHER or the usual taxonomy assignment. I will generate both and leave it up to the user.
if(!exists("taxa.print")) {
print("Running the assignTaxonomy taxonomy Assignment")
taxa <- assignTaxonomy(seqtab.nochim, paste0(dada2_path, "/sh_general_release_dynamic_04.02.2020.fasta"), multithread=threads, tryRC=TRUE, minBoot=80,outputBootstraps=TRUE)
write.csv(taxa$boot, paste0(save_path, "/TaxonomyBootstraps.csv"), row.names=F)
taxa.print <- taxa$tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
###table(is.na(taxa.print[,"Genus"]))
} else {print("DADA2 Taxonomy already generated")}
write.table(x=taxa.print, file=paste0(save_path,"/", ns_dir , "_tax_table.txt"), sep="\t", col.names=T, row.names=T, quote=F)

save.image(paste0(save_path, "/dada2_result_object.RData"))

## Handoff to phyloseq:
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs

orig <- seqtab.nochim
rownames(orig) <- sub("_R1_filt.fastq.gz","", rownames(orig))
dada2_phyobj <- phyloseq(otu_table(orig, taxa_are_rows=FALSE), tax_table(taxa$tax), sample_data(map_file), DNAStringSet(seqs))
dada2_phyobj
saveRDS(dada2_phyobj, paste0(save_path, "/dada2_phy_obj_raw.rds"))
