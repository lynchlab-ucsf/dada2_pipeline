# Code to remove signal arising from negative controls.
# Katie McCauley
## Added command line options.

ntc_filter <- function(phy, ntcs = c("NTC","EMPTY"), pct_neg=0.15, pct_samp=0.15, stat="mean", exclude_ntcs=NULL) {
require(phyloseq)
require(ggplot2)
require(gridExtra)

sample_data(phy)$SampleName <- sample_names(phy)

print("Generating First PCoA Plot")
ord <- ordinate(phy, method="PCoA", distance="bray")
pcoa.plt <- plot_ordination(phy, ord, color="SampleType", justDF=TRUE)
bray.pcoa <- ggplot(pcoa.plt, aes(x=Axis.1, y=Axis.2, color=SampleType)) + 
  geom_point() + 
  geom_text(data=subset(pcoa.plt, SampleType %in% ntcs), aes(label=SampleName)) +
  ggtitle("Bray Curtis")
print("Generating Second PCoA Plot")
ord <- ordinate(phy, method="PCoA", distance="canberra")
pcoa.plt <- plot_ordination(phy, ord, color="SampleType", justDF=TRUE)
can.pcoa <- ggplot(pcoa.plt, aes(x=Axis.1, y=Axis.2, color=SampleType)) +                         
  geom_point() + 
  geom_text(data=subset(pcoa.plt, SampleType %in% ntcs), aes(label=SampleName)) +
  ggtitle("Canberra")

if(!is.null(exclude_ntcs)) {
have.neg.list <- 1
neg.list <- exclude_ntcs
print("Negative Control Exclusion list found. Removing the following negative controls from consideration as negative controls.")
print(neg.list$V1)
phy <- subset_samples(phy, !sample_names(phy) %in% neg.list$V1)
}

if(sum(phy@sam_data$SampleType %in% ntcs) == 0) print("NO NEGATIVE CONTROLS FOUND. IF THIS IS NOT EXPECTED, please ensure that a 'SampleType' variable exists in your mapping file and it has values of either 'NTC' or 'EMPTY' to identify negative control samples. More details can be found in the 16s Pipeline document.")

ntc_samps <- phy@sam_data$SampleType %in% ntcs

neg.dat <- prune_samples(ntc_samps, phy)
samp.dat <- prune_samples(!ntc_samps, phy)

print("Negative Controls Found:")
print(sample_names(neg.dat))
print("If this is unexpected, this script uses the SampleType variable in the mapping file to identify Negative Controls")

negsA <- taxa_sums(otu_table(neg.dat)>0)/nsamples(neg.dat)
sampsA <- taxa_sums(otu_table(samp.dat)>0)/nsamples(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(samp.dat)@.Data, by.x="#OTU ID",by.y=0)
dat2$Genus[is.na(dat2$Genus)] <- "NA"


pre.cleaning <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.15 & samps < 0.15), aes(label=Genus), size=2)

print(paste0("Filtering taxa in greater than or equal to ", pct_neg*100,"% of negative controls and less than or equal to ", pct_samp*100,"% of samples."))
neg.prop <- as.numeric(pct_neg)
samp.prop <- as.numeric(pct_samp)

outright.drop <- dat2$`#OTU ID`[dat2$negs >= neg.prop & dat2$samps <= samp.prop]
print(dat2[dat2$`#OTU ID` %in% outright.drop, !names(dat2) %in% c("Species")])

all.datA <- prune_taxa(!taxa_names(phy) %in% outright.drop, phy)

neg.dat2 <- prune_taxa(!taxa_names(phy) %in% outright.drop, neg.dat)

sub.in.NTC <- ceiling(apply(t(otu_table(neg.dat2)), 1, stat))

print(paste0("Performing subtraction of the ", stat," of remaining negative control reads from samples"))
otu_sub <- otu_table(all.datA)
for(i in names(sub.in.NTC)) {
  otu_sub[, i] <- otu_sub[, i] - sub.in.NTC[i]
  otu_sub[, i][otu_sub[, i] < 0] <- 0
}
otu_submean <- all.datA
otu_table(otu_submean) <- otu_table(t(otu_sub), taxa_are_rows = TRUE)

neg.dat <- prune_samples(ntc_samps, otu_submean)
samp.dat <- prune_samples(!ntc_samps, otu_submean)

negsA <- taxa_sums(otu_table(neg.dat)>0)/nsamples(neg.dat)
sampsA <- taxa_sums(otu_table(samp.dat)>0)/nsamples(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(samp.dat)@.Data, by.x="#OTU ID",by.y=0)
dat2$Genus[is.na(dat2$Genus)] <- "NA"

post.clean <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.15 & samps < 0.15), aes(label=Genus), size=2)

plt <- grid.arrange(bray.pcoa, can.pcoa, pre.cleaning, post.clean)
samp.dat <- prune_taxa(taxa_sums(samp.dat) > 0, samp.dat)
samp.dat <- prune_samples(sample_sums(samp.dat) >0, samp.dat)
plot(plt)

return(samp.dat)
}
