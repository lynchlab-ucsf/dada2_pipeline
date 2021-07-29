# Code to remove signal arising from negative controls.
# Katie McCauley
## Added command line options. 
writeLines(
"Code to remove signal arising from negative controls.
Written by Katie McCauley
Last Updated: 5/17/21

When run as part of the complete_16s_pipeline, the default options are already chosen for you. If you need to go back and change defaults based on your results, below is the syntax. Unfortunately, it is difficult to apply named arguments in Rscript functions, so this will need to suffice. See example usage below.

## Option 1: Phyloseq object name/location
## Option 2: Outright remove if sequence is greater than this proportion in negative controls
## Option 3: Outright remove if sequence is less than this proportion in samples
## Option 4: To subtract mean or max in negative controls from samples

## Example: Rscript NegControlRemoval_phy.R dada2_phy_obj_raw.rds 0.15 0.15 mean
")

library(phyloseq)
library(ggplot2)
setwd(".")
argv <- commandArgs(TRUE)
phy_filt_tree <- readRDS(argv[1])
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

print(paste0("Filtering taxa in greater than ", argv[2],"% of negative controls and less than ", argv[3],"% of samples."))
neg.prop <- as.numeric(argv[2])
samp.prop <- as.numeric(argv[3])

outright.drop <- dat2$`#OTU ID`[dat2$negs > neg.prop & dat2$samps < samp.prop]

#ord = outright.drop

all.datA <- subset_taxa(phy_filt_tree, !taxa_names(phy_filt_tree) %in% outright.drop)
####

#otu.sub
neg.dat2 <- subset_taxa(neg.dat, taxa_sums(neg.dat) >0 & !taxa_names(neg.dat) %in% outright.drop)

sub.in.NTC <- ceiling(apply(t(otu_table(neg.dat2)), 1, argv[4]))

print(paste0("Performing subtraction of the ", argv[4]," of remaining negative control reads from samples"))
otu_sub <- otu_table(all.datA)
for(i in names(sub.in.NTC)) {
  otu_sub[, i] <- otu_sub[, i] - sub.in.NTC[i]
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
