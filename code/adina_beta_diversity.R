library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
 
setwd('~/Documents/StripsDiversityWebsite/')
worle <- readRDS("data/Worle_curated.RDS")

# Looking at only soil / no manure treatment for the biodiversity paper
soil <- subset_samples(worle, matrix == c("soil"))
head(sample_data(soil))
colnames(sample_data(soil))
#[1] "unique_id"        "experiment"       "matrix"           "treatment"       
#[5] "plot"             "sample_day"       "depth"            "in_plot_location"
#[9] "block"            "strip"            "manure_treatment" "soil_type"       
sample_data(soil)$plot <- factor(sample_data(soil)$plot)
soil_nm <- subset_samples(soil, treatment == "no_manure_strip")
times <- unique(data.frame(sample_data(soil)$sample_day))

#stats based on each time point
# I removed border since it doesn't seem to be consistently strippy or croppy
# but it is often in the mix betwweenthe two, you can add it to the pcoa and look
# when the soil is dryer, it seems like the soil type makes a bigger difference, and also we see more plot by plot differences
# This can be repeated for depths 1 and 2, there are some interseting differences for depths
# depth 2 i think is more similar overall
day_list <- unique(sample_data(soil)$sample_day)
for (i in 1:length(day_list)){
  baseline <- subset_samples(soil, sample_day == day_list[i] & treatment == "no_manure_strip")
  baseline2 <- subset_samples(baseline, depth == "1" & soil_type != "border")
  bray <- phyloseq::distance(physeq = baseline2, method = "bray")
  print(day_list[i])
  print(adonis(bray~soil_type+as.factor(plot), data = data.frame(sample_data(baseline2)), permutations = 999, methods = "bray"))
  #print(adonis(bray~in_plot_location+as.factor(plot), data = data.frame(sample_data(baseline2)), permutations = 999, methods = "bray"))
  #no significance by in plot location, this may be an interesting result
  }


baseline <- subset_samples(soil, treatment == "no_manure_strip")
baseline2 <- subset_samples(baseline, depth == "1" & soil_type != "border")
baseline2b <- subset_samples(baseline2, in_plot_location %in% c("s1", "s2", "s8", "s9"))
#ord <- ordinate(baseline2, "PCoA", "bray") 
ord <- ordinate(baseline2b, "PCoA", "bray") 
p1 = plot_ordination(baseline2b, ord, color="soil_type") 
p1+ geom_point(size=2) +facet_grid(sample_day~plot)
#cca makes things very pretty but not sure if its approriate

#Strips are more heterogeneous than crop
bray <- phyloseq::distance(baseline2, "bray")
sampledf <- data.frame(sample_data(baseline2))
beta <- betadisper(bray, sampledf$soil_type)
pmod <- permutest(beta, permutations = 999, pairwise=TRUE)
pmod

# Betadiversity measures
library(data.table)
ord <- ordinate(baseline2, "PCoA", "bray") 
ord.df <- data.table(ord$vectors, keep.rownames = TRUE)
setnames(ord.df, "rn", "SampleID")
ord.df$SampleID
ggplot(ord.df, aes())
sdt = data.table(as(sample_data(baseline2), "data.frame"),
                  keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
c <- subset(sdt, soil_type == "crop")
combined<- merge(ord.df, sdt, by="SampleID")
head(combined)

combined2 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1=mean(Axis.1), SE = sd(Axis.1)/sqrt(length(Axis.1)))
limits<-aes(ymin=MEAN.1-SE, ymax=MEAN.1+SE)
ggplot(combined2, aes(x=sample_day, y=MEAN.1, group = soil_type, color=soil_type))+geom_point()+geom_line()+geom_errorbar(limits, width=0)
#shows more distance from baseline in PC1 of strips compared to rop
# 10/19/20 Ask Adina if MEAN.1 should have called Axis.2, it had Axis.1
combined2 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1=mean(Axis.2), SE = sd(Axis.2)/sqrt(length(Axis.2)))
limits<-aes(ymin=MEAN.1-SE, ymax=MEAN.1+SE)
ggplot(combined2, aes(x=sample_day, y=MEAN.1, group = soil_type, color=soil_type))+geom_point()+geom_line()+geom_errorbar(limits, width=0)
#shows more distance from baseline in PC2 of strips compared to rop

#Unifracing
#I dont have the tree but I dont think its going to be sueful to show here.



# deseq code but i didn't clean this up



baseline <- subset_samples(soil, sample_day == "Baseline" & treatment == "no_manure_strip")
baseline2 <- subset_samples(baseline, depth == "1")

diagdds = phyloseq_to_deseq2(baseline2, ~ soil_type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)

alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(baseline2)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

s <- prune_taxa(c("TACGTAGGGGGCGAGCGTTGTCCGGAGTTACTGGGCGTAAAGGGCATGCAGGCGGTCGGACGCGTTTCGGGTGACAGCGGGTGGCTTAACTACCCGAGTGGTCGAAAGACGGTTTGACTTGAGGACGCGAGAGGGACACGGAATTCCGGGTGGAGTGGTGAAATGCGTAGAGATCCGGAGGAACACCGAAGGCGAAGGCAGTGTCCTGGCGCGTAACTGACGCTGAGGTGCGAAAGCCAGGGGAGCGAACGGG"), baseline2)
plot_bar(s)

library("devtools")
library(microbiome)
pseq.rel <- microbiome::transform(baseline2, "compositional")
strips <- subset_samples(baseline2, soil_type == "strip")
strips = filter_taxa(strips, function(x) sum(x) > 0, TRUE)

crops <- subset_samples(baseline2, soil_type == "crop")
crops = filter_taxa(crops, function(x) sum(x) > 0, TRUE)

strips.core <- core(strips, detection = 0, prevalence = 0.5)
crops.core <- core(crops, detection = 0, prevalence = 0.5)
strips.core.rel <-  microbiome::transform(strips.core, "compositional")
crops.core.rel <-  microbiome::transform(crops.core, "compositional")

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(strips.core.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")
plot_core(crops.core.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

rownames(otu_table(crops.core.rel)) %in% rownames(otu_table(strips.core.rel))
rownames(otu_table(strips.core.rel)) %in% rownames(otu_table(crops.core.rel))
