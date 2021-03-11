#! /usr/bin/Rscript

## Will need to set this up so that it's in the correct directory. Will likely involve a bash script...
library(phyloseq)
phy <- readRDS("dada2_phy_obj_raw.rds")

#phy_filtB <- subset_taxa(phy, Kingdom %in% "Bacteria") ## Removing after discussion with Sue
phy_filtP <- prune_taxa(taxa_sums(phy) > 0.000001*sum(taxa_sums(phy)), phy)
phy_filtP

## Begin Tree Generation
library(phangorn)
library(msa)
library(DECIPHER)

print("Running Alignment")
alignment <- AlignSeqs(DNAStringSet(refseq(phy_filtP)), anchor=NA)
print("Running phangorn")
phang.align <- as.phyDat(as(alignment, "matrix"), type="DNA")
print("Making distance matrix")
dm <- dist.ml(phang.align)
print("Making tree from distance matrix")
treeNJ <- NJ(dm)
print("Fitting")
fit <- pml(treeNJ, data=phang.align)
write.tree(fit$tree, "dada2_init_tree.tre")
phy_filtT <- merge_phyloseq(phy_filtP, phy_tree(fit$tree))
saveRDS(phy_filtT, file="dada2_phy_init_tree.rds")

fitGTR <- update(fit, k=4, inv=0)
print("Optimizing Fit")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
          rearrangement = "stochastic", control = pml.control(trace = 0))
write.tree(fitGTR$tree, "dada2_optim_tree.tre")
phy_filtF <- merge_phyloseq(phy_filtP, phy_tree(fitGTR$tree))
saveRDS(phy_filtF, file="dada2_phy_pruned_wtree.rds")

