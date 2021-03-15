library(vegan)
library(phyloseq)
library(dplyr)

args <- commandArgs(TRUE)

phy <- readRDS(as.character(args[1]))

set.max <- as.numeric(args[2])

uselevels <- data.frame(read.depth=sample_sums(phy)) %>%
  mutate(groups=factor(ntile(read.depth, 25))) %>% 
  group_by(groups) %>%
  summarize(min(read.depth))

depth <- samples.kept <- c()
for(i in 1:nrow(uselevels[,2])) {
  depth[i] <- uselevels[i,2]
  samples.kept[i] <- sum(sample_sums(phy) > depth[i])
}

pdf("Alpha_Rarefaction_Curve.pdf",height=8, width=11)
plt <- rarecurve(t(otu_table(phy)), step=floor(max(sample_sums(phy))/25),cex=0.5)
abline(v=depth)
dev.off()

print(cbind(depth, samples.kept))
