library("DESeq2")

worle <- readRDS("data/Worle_curated.RDS")
soil <- subset_samples(worle, matrix == c("soil"))
baseline <- subset_samples(soil, sample_day == "Baseline" & treatment == "no_manure_strip")
baseline2 <- subset_samples(baseline, depth == "1") # If changed to depth 2, no seqs are sig nif dif

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
plot_bar(s, fill = "soil_type")

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
