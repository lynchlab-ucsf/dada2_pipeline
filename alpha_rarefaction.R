library(vegan)
library(phyloseq)

args <- commandArgs(TRUE)

phy <- readRDS(as.character(args[1]))

depth <- samples.kept <- c()
for(i in 1:25) {
  depth[i] <- floor(max(sample_sums(phy))/25)*i
  samples.kept[i] <- sum(sample_sums(phy) > depth[i])
}

pdf("Alpha_Rarefaction_Curve.pdf",height=8, width=11)
plt <- rarecurve(t(otu_table(phy)), step=floor(max(sample_sums(phy))/25),cex=0.5)
abline(v=depth)
dev.off()

print(cbind(depth, samples.kept))
