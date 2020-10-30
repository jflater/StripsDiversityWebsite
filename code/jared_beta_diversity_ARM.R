library(vegan)
library(plyr)
library(tidyverse)
library(phyloseq)
library(phylosmith)
setwd('~/Documents/StripsDiversityWebsite/')
worle <- readRDS("data/Arm_curated.RDS")

# Two baseline phyloseq, one for each depth 
baseline1 <- subset_samples(worle, treatment == "no_manure_strip" & depth == "1" & soil_type != "border" & matrix == "soil") %>%
  filter_taxa(function(x) sum(x) >= 1, T)
baseline2 <- subset_samples(worle, treatment == "no_manure_strip" & depth == "2" & soil_type != "border" & matrix == "soil")  %>%
  filter_taxa(function(x) sum(x) >= 1, T)

pcoa_phyloseq(baseline1, treatment = c("soil_type", "sample_day"), x = 1, y = 2, method = 'bray', circle = 0.95, colors = 'default', labels = NULL)
pcoa_phyloseq(baseline2, treatment = c("soil_type", "sample_day"), x = 1, y = 2, method = 'bray', circle = 0.95, colors = 'default', labels = NULL)
nmds_phyloseq(baseline1, treatment = c("soil_type", "sample_day"), method = 'bray', circle = 0.95, colors = 'default', labels = NULL)
nmds_phyloseq(baseline2, treatment = c("soil_type", "sample_day"), method = 'bray', circle = 0.95, colors = 'default', labels = NULL)

# Plot PCoA ordinations
# D1
ord1 <- ordinate(baseline1, "PCoA", "bray")
p1d1 <- plot_ordination(baseline1, ord1, color = "soil_type") +
  geom_point(size = 2) +
  facet_grid(sample_day ~ plot)
p1d1

# D2
ord2 <- ordinate(baseline2, "PCoA", "bray")
p2d2 <- plot_ordination(baseline2, ord2, color = "soil_type") +
  geom_point(size = 2) +
  facet_grid(sample_day ~ plot)
p2d2

#stats
#d1
day_list <- unique(sample_data(worle)$sample_day)
for (i in 1:length(day_list)) {
  baseline.1 <- subset_samples(worle, sample_day == day_list[i] & treatment == "no_manure_strip" & depth == "1" & soil_type != "border" & matrix == "soil") %>%
    filter_taxa(function(x) sum(x) >= 1, T) 
  bray <- phyloseq::distance(physeq = baseline.1, method = "bray")
  print(day_list[i])
  print(adonis(bray ~ soil_type + as.factor(plot), data = data.frame(sample_data(baseline.1)), permutations = 999, methods = "bray"))
}
#d2
for (i in 1:length(day_list)) {
  baseline.2 <- subset_samples(worle, sample_day == day_list[i] & treatment == "no_manure_strip" & depth == "2" & soil_type != "border" & matrix == "soil") %>%
    filter_taxa(function(x) sum(x) >= 1, T) 
  bray <- phyloseq::distance(physeq = baseline.2, method = "bray")
  print(day_list[i])
  print(adonis(bray ~ soil_type + as.factor(plot), data = data.frame(sample_data(baseline.2)), permutations = 999, methods = "bray"))
}

# Betadiversity measures d1
baseline1 <- subset_samples(worle, treatment == "no_manure_strip" & depth == "1" & soil_type != "border" & matrix == "soil") %>%
  filter_taxa(function(x) sum(x) >= 1, T)
library(data.table)
ord <- ordinate(baseline1, "PCoA", "bray") 
ord.df <- data.table(ord$vectors, keep.rownames = TRUE)
setnames(ord.df, "rn", "SampleID")
ord.df$SampleID
ggplot(ord.df, aes())
sdt = data.table(as(sample_data(baseline1), "data.frame"),
                 keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
c <- subset(sdt, soil_type == "crop")
combined <- merge(ord.df, sdt, by = "SampleID")
head(combined)

combined2 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1 = mean(Axis.1), SE = sd(Axis.1)/sqrt(length(Axis.1)))
limits <- aes(ymin = MEAN.1 - SE, ymax = MEAN.1 + SE)
ggplot(combined2, aes(x = sample_day, y = MEAN.1, group = soil_type, color = soil_type)) + geom_point() + geom_line() + geom_errorbar(limits, width = 0)
ggsave("d1")
#shows more distance from baseline in PC1 of strips compared to crop
# 10/19/20 Ask Adina if MEAN.1 should have called Axis.2, it had Axis.1

combined3 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1 = mean(Axis.2), SE = sd(Axis.2)/sqrt(length(Axis.2)))
limits <- aes(ymin = MEAN.1 - SE, ymax = MEAN.1 + SE)
ggplot(combined3, aes(x = sample_day, y = MEAN.1, group = soil_type, color = soil_type)) + geom_point() + geom_line() + geom_errorbar(limits, width = 0)
#shows more distance from baseline in PC2 of strips compared to crop

# Betadiversity measures d2
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
combined <- merge(ord.df, sdt, by = "SampleID")
head(combined)

combined2 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1 = mean(Axis.1), SE = sd(Axis.1)/sqrt(length(Axis.1)))
limits <- aes(ymin = MEAN.1 - SE, ymax = MEAN.1 + SE)
ggplot(combined2, aes(x = sample_day, y = MEAN.1, group = soil_type, color = soil_type)) + geom_point() + geom_line() + geom_errorbar(limits, width = 0)
#shows more distance from baseline in PC1 of strips compared to crop
# 10/19/20 Ask Adina if MEAN.1 should have called Axis.2, it had Axis.1

combined3 <- plyr::ddply(combined, .(sample_day, soil_type),summarise, MEAN.1 = mean(Axis.2), SE = sd(Axis.2)/sqrt(length(Axis.2)))
limits <- aes(ymin = MEAN.1 - SE, ymax = MEAN.1 + SE)
ggplot(combined3, aes(x = sample_day, y = MEAN.1, group = soil_type, color = soil_type)) + geom_point() + geom_line() + geom_errorbar(limits, width = 0)
#shows more distance from baseline in PC2 of strips compared to crop