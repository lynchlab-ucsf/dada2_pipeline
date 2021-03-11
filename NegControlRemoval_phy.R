# Code to remove negative controls.
# Katie McCauley
# Is it possible to develop this code so that it plays well with the pipeline so far.
#Maybe use the input CSV to determine which samples belong to each run and do run-specific negative control filtering?

library(phyloseq)
library(ggplot2)
setwd(".")
phy_filt_tree <- readRDS("dada2_phy_pruned_wtree.rds")
myNTCs=c("NTC","EMPTY")
sample_data(phy_filt_tree)$SampleName <- sample_names(phy_filt_tree)

print("Generating First PCoA Plot")
ord <- ordinate(phy_filt_tree, method="PCoA", distance="bray")
pcoa.plt <- plot_ordination(phy_filt_tree, ord, color="SampleType", justDF=TRUE)
bray.pcoa <- ggplot(pcoa.plt, aes(x=Axis.1, y=Axis.2, color=SampleType)) + 
  geom_point() + 
  geom_text(data=subset(pcoa.plt, SampleType %in% myNTCs), aes(label=SampleName))
print("Generating Second PCoA Plot")
ord <- ordinate(phy_filt_tree, method="PCoA", distance="canberra")
pcoa.plt <- plot_ordination(phy_filt_tree, ord, color="SampleType", justDF=TRUE)
can.pcoa <- ggplot(pcoa.plt, aes(x=Axis.1, y=Axis.2, color=SampleType)) +                         
  geom_point() + 
  geom_text(data=subset(pcoa.plt, SampleType %in% myNTCs), aes(label=SampleName))

if(file.exists("ExcludeNegatives.txt")) {
have.neg.list <- 1
neg.list <- read.table("ExcludeNegatives.txt", header=F, sep="\t", comment="")
print("Negative Control Exclusion list found. Removing the following negative controls from consideration as negative controls.")
print(neg.list$V1)
phy_filt_tree <- subset_samples(phy_filt_tree, !sample_names(phy_filt_tree) %in% neg.list$V1)
}

if(sum(phy_filt_tree@sam_data$SampleType %in% myNTCs) == 0) print("NO NEGATIVE CONTROLS FOUND. IF THIS IS NOT EXPECTED, please ensure that a 'SampleType' variable exists in your mapping file and it has values of either 'NTC' or 'EMPTY' to identify negative control samples. More details can be found in the 16s Pipeline document.")

neg.dat <- subset_samples(phy_filt_tree, SampleType %in% myNTCs)
samp.dat <- subset_samples(phy_filt_tree, !SampleType %in% myNTCs)

print("Negative Controls Found:")
print(sample_names(neg.dat))
print("If this is unexpected, this script uses the SampleType variable in the mapping file to identify Negative Controls")

negsA <- taxa_sums(otu_table(neg.dat)>0)/nsamples(neg.dat)
sampsA <- taxa_sums(otu_table(samp.dat)>0)/nsamples(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(samp.dat)@.Data, by.x="#OTU ID",by.y=0)
dat2$Genus[is.na(dat2$Genus)] <- "NA"

library(ggplot2)
pre.cleaning <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.15), aes(label=Genus)) +
  labs(tag="A")
#ggsave("Pre_Cleaning_Figure.pdf", pre.cleaning, device="pdf", height=6, width=7)

print("Filtering taxa in greater than 15% of negative controls and less than 15% of samples")

outright.drop <- dat2$`#OTU ID`[dat2$negs > 0.15 & dat2$samps < 0.15]

#ord = outright.drop

all.datA <- subset_taxa(phy_filt_tree, !taxa_names(phy_filt_tree) %in% outright.drop)
####

#otu.sub
neg.dat2 <- subset_taxa(neg.dat, taxa_sums(neg.dat) >0 & !taxa_names(neg.dat) %in% outright.drop)

mean.in.NTC <- ceiling(apply(t(otu_table(neg.dat2)), 1, mean))

print("Performing subtraction of the mean of remaining negative control reads from samples")
otu_sub <- otu_table(all.datA)
for(i in names(mean.in.NTC)) {
  otu_sub[, i] <- otu_sub[, i] - mean.in.NTC[i]
  otu_sub[, i][otu_sub[, i] < 0] <- 0
}
otu_submean <- all.datA
otu_table(otu_submean) <- otu_table(t(otu_sub), taxa_are_rows = TRUE)

neg.dat <- subset_samples(otu_submean, SampleType %in% myNTCs)
samp.dat <- subset_samples(otu_submean, !SampleType %in% myNTCs)

negsA <- taxa_sums(otu_table(neg.dat)>0)/nsamples(neg.dat)
sampsA <- taxa_sums(otu_table(samp.dat)>0)/nsamples(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(samp.dat)@.Data, by.x="#OTU ID",by.y=0)
dat2$Genus[is.na(dat2$Genus)] <- "NA"

library(ggplot2)
post.clean <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.15), aes(label=Genus)) +
  labs(tag="A")

library(gridExtra)
if(exists("have.neg.list")) {
extname <- "_negex"
} else {
extname <- ""
}
ggsave(paste0("Combined_Cleaning_Figure",extname,".png"), grid.arrange(bray.pcoa, can.pcoa, pre.cleaning, post.clean, nrow=2), device="png", height=11, width=14)
samp.dat <- subset_taxa(samp.dat, taxa_sums(samp.dat) > 0)
samp.dat <- subset_samples(samp.dat, sample_sums(samp.dat) >0)
saveRDS(samp.dat, paste0("phyloseq_noneg",extname,".rds"))
