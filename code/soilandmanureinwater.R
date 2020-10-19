# we need list of ASVs from manure and from crop portion of each plot
# We will then look to see if strips filter these out of the water samples
library(tidyverse)
library(phyloseq)
library(venn)
library(ggpubr)
library(phylosmith)
setwd("~/Documents/StripsDiversityWebsite/")
worle <- readRDS("data/Worle_curated.RDS") 

taxa_names(worle) <- paste0("ASV", seq(ntaxa(worle)))

worle_manure <- subset_samples(worle, matrix == c("manure")) 
worle_soil <- subset_samples(worle, matrix == c("soil") & soil_type == c("crop") & sample_day == c("Baseline"))
worle_strips <- subset_samples(worle, matrix == c("soil") & soil_type == c("strip") & sample_day == c("Baseline"))

worle_manure_mintax1 <- worle_manure %>%
  filter_taxa(function(x) sum(x) >= 1, T) 
worle_crop_soil_mintax1 <- worle_soil %>%
  filter_taxa(function(x) sum(x) >= 1, T)
worle_strip_soil_mintax1 <- worle_strips %>%
  filter_taxa(function(x) sum(x) >= 1, T)

worle_manure_asvs <- taxa_names(worle_manure_mintax1)
worle_crop_soil_asvs <- taxa_names(worle_crop_soil_mintax1)
worle_strip_soil_asvs <- taxa_names(worle_strip_soil_mintax1)

worle_vvv_diag <- venn(list("Manure_ASVs" = worle_manure_asvs, "Crop_ASVs" = worle_crop_soil_asvs, "Strip_ASVs" = worle_strip_soil_asvs))
worle_vvv_diag
worle_manure_persitors <- attr(worle_vvv_diag, "intersections")$Manure_ASVs
worle_crop_persistors <- attr(worle_vvv_diag, "intersections")$Crop_ASVs
worle_strip_persistors <- attr(worle_vvv_diag, 'intersection')$Strip_ASVs
worle_cropstrip_persistors <- attr(worle_vvv_diag, 'intersections')$`Crop_ASVs:Strip_ASVs`


worle.water.phy <- subset_samples(worle, matrix == "water" & unique_id != "Comp-from-P8-7-10-26-17") %>%
  filter_taxa(function(x) sum(x) >= 1, T) 

tax_association <- tax_table(worle.water.phy) %>%
  data.frame() %>%
  rownames_to_column("ASV") %>%
  mutate(ASV_Association = ifelse(ASV %in% worle_manure_persitors, "Manure_associated",
                                         ifelse(ASV %in% worle_crop_persistors, "Crop_associated", 
                                                ifelse(ASV %in% worle_strip_persistors, "Strip_associated","No_unique_association")))) %>%
  column_to_rownames("ASV") %>%
  as.matrix()

#taxa_names(worle.water.phy)
#rownames(tax_association)
tax_table(worle.water.phy) <- tax_association

####
# Testing relative abundance plot with full water dataset
phylogeny_profile(worle.water.phy, classification = 'ASV_Association', treatment = c("treatment"), merge = TRUE, relative_abundance = TRUE)



# # worle.water.phy.merged <- merge_samples(worle.water.phy, "trtandwater")
# # sample_data(worle.water.phy.merged)$trtandwater <- levels(sample_data(worle.water.phy.merged)$trtandwater)
# # sample_data(worle.water.phy.merged)
# 
# worle.manure.water <- prune_taxa(worle_manure_persitors, worle.water.phy) 
# worle.manure.water
# worle.soil.water <- prune_taxa(worle_crop_persistors, worle.water.phy) 
# worle.soil.water
# 
# phylogeny_profile(worle.manure.water, classification = 'Phylum', treatment = c("treatment"), merge = TRUE, relative_abundance = F) + 
#   ggtitle("Manure ASVs in water")
# phylogeny_profile(worle.soil.water, classification = 'Phylum', treatment = c("treatment"), merge = TRUE, relative_abundance = F) +
#   ggtitle("Soil ASVs in water")
# 
# worle.manure.water.melted <- worle.manure.water %>%
#   psmelt() %>%
#   separate(Sample, c('trt', 'A'), sep = " ")
# worle.soil.water.melted <- worle.soil.water %>%
#   psmelt() %>%
#   separate(Sample, c('trt', 'A'), sep = " ")
# 
# manure <- ggplot(data = worle.manure.water.melted, aes(x = water_sample, y = Abundance)) + 
#   geom_bar(aes(fill = Phylum), stat = "identity") +
#   facet_grid(~trt) +
#   theme_pubclean() + 
#   ggtitle("Manure ASVs in water")
# 
# soil <- ggplot(data = worle.soil.water.melted, aes(x = water_sample, y = Abundance)) + 
#   geom_bar(aes(fill = Phylum), stat = "identity") +
#   facet_grid(~trt) +
#   theme_pubclean() +
#   ggtitle("Soil ASVs in water")
# 
# # Below to make plots prettier
# phylalist <- data.frame(tax_table(worle),row.names = NULL) %>%
#   select(Phylum) %>%
#   unique() 
# phylalist$Phylum <- as.character(phylalist$Phylum)
# phylalist$Phylum[is.na(phylalist$Phylum)] <- "Unclassified"
# 
# # this package will generate a pallette based on number and desired colors
# library(colorspace)
# colors <- sequential_hcl(n_distinct(phylalist), palette = "viridis") %>%
#   sample() %>%
#   setNames(phylalist$Phylum)
# 
# addSmallLegend <- function(myPlot, pointSize = 1, textSize = 5, spaceLegend = 0.5) {
#   myPlot +
#     guides(shape = guide_legend(override.aes = list(size = pointSize)),
#            color = guide_legend(override.aes = list(size = pointSize))) +
#     theme(legend.title = element_text(size = textSize), 
#           legend.text  = element_text(size = textSize),
#           legend.key.size = unit(spaceLegend, "lines"))
# }
# addSmallLegend(manure) + scale_color_manual(aesthetics = "fill", values = colors)
# addSmallLegend(soil) + scale_color_manual(aesthetics = "fill", values = colors)
