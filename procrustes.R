library(vegan)
library(phyloseq)

args <- commandArgs(TRUE)
print("Reading in Phyloseq Objects")
ps1 <- readRDS(as.character(args[1]))
ps2 <- readRDS(as.character(args[2]))

#Ensure the same samples are in the objects
print("Subsetting the objects to the same set of samples")
ps1.s <- subset_samples(ps1, sample_names(ps1) %in% sample_names(ps2))
ps2.s <- subset_samples(ps2, sample_names(ps2) %in% sample_names(ps1))

print("Performing ordination")
pcoa1 <- ordinate(ps1.s, method="PCoA", distance="bray")
pcoa2 <- ordinate(ps2.s, method="PCoA", distance="bray")

print("Performing procrustes")
proc <- protest(pcoa1$vectors, pcoa2$vectors)
print(proc)
