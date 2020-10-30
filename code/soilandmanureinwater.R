# we need list of ASVs from manure and from crop portion of each plot
# We will then look to see if strips filter these out of the water samples
library(tidyverse)
library(phyloseq)
library(venn)
library(ggpubr)
library(phylosmith)
library(microbiome)
setwd("~/Documents/StripsDiversityWebsite/")
worle <- readRDS("data/Worle_curated.RDS") 

taxa_names(worle) <- paste0("ASV", seq(ntaxa(worle)))
worle1 <- subset_samples(worle, matrix %in% c("manure", "soil")) %>%
  filter_taxa(function(x) sum(x) >= 2, T)

max(sample_sums(worle1))
min(sample_sums(worle1))
nsamples(worle1)
ntaxa(worle1)
summarize_phyloseq(worle1)

tdt <- sample_sums(worle1) %>%
  enframe(name = "ID", value = "Sample_sum")

tdt <- as_tibble(sample_data(worle1)) %>%
  cbind(., sample_sums(worle1)) %>%
  dplyr::group_by(matrix) %>%
  dplyr::summarize(Mean_sample_sums = mean(`sample_sums(worle)`, na.rm = TRUE), n_samples = n()) %>%
  write_csv("data/SSUMSWORLE.csv")


worle_manure <- subset_samples(worle, matrix == c("manure")) 
# worle_soil <- subset_samples(worle, matrix == c("soil") & soil_type == c("crop") & sample_day == c("Baseline"))
# worle_strips <- subset_samples(worle, matrix == c("soil") & soil_type == c("strip") & sample_day == c("Baseline"))

worle_soil <- subset_samples(worle, matrix == c("soil") & soil_type == c("crop") & treatment == c("no_manure_strip"))
worle_strips <- subset_samples(worle, matrix == c("soil") & soil_type == c("strip") & treatment == c("no_manure_strip"))

worle_manure_mintax1 <- worle_manure %>%
  filter_taxa(function(x) sum(x) >= 2, T) 
ntaxa(worle_manure_mintax1)
phylogeny_profile(worle_manure_mintax1, classification = "Phylum", merge = TRUE, relative_abundance = TRUE)
worle_crop_soil_mintax1 <- worle_soil %>%
  filter_taxa(function(x) sum(x) >= 2, T)
worle_strip_soil_mintax1 <- worle_strips %>%
  filter_taxa(function(x) sum(x) >= 2, T)

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
  filter_taxa(function(x) sum(x) >= 2, T) 

tax_association <- tax_table(worle.water.phy) %>%
  data.frame() %>%
  rownames_to_column("ASV") %>%
  mutate(ASV_Association = ifelse(ASV %in% worle_manure_persitors, "Manure_associated",
                                         ifelse(ASV %in% worle_crop_persistors, "Crop_associated", 
                                                ifelse(ASV %in% worle_strip_persistors, "Strip_associated","No_unique_association")))) %>%
  column_to_rownames("ASV") %>%
  as.matrix()


tax_table(worle.water.phy) <- tax_association

####
# Testing relative abundance plot with full water dataset
plot <- phylogeny_profile(worle.water.phy, classification = 'ASV_Association', treatment = c("treatment"), merge = TRUE, relative_abundance = TRUE)
plot
ggsave("images/manureinwater.png", plot = plot, height = 3, width = 11)
data <- plot$data %>%
  group_by(plot, treatment, ASV_Association) %>%
  summarise_at(vars(-OTU, -Sample, -unique_id, -experiment, -matrix, -sample_day, -depth, -in_plot_location, -block, -strip, -manure_treatment, -soil_type), funs(mean(., na.rm=TRUE)))
data
write_csv(x = data, file = "data/manuretaxainwater.csv")
#### Manure ASVs in soil
worle.soil.phy <- subset_samples(worle, matrix == "soil" & treatment %in% c("manured_control", "manured_strip") & sample_day != c("Baseline")) %>%
  filter_taxa(function(x) sum(x) >= 2, T) 

tax_association <- tax_table(worle.soil.phy) %>%
  data.frame() %>%
  rownames_to_column("ASV") %>%
  mutate(ASV_Association = ifelse(ASV %in% worle_manure_persitors, "Manure_associated",
                                  ifelse(ASV %in% worle_crop_persistors, "Crop_associated", 
                                         ifelse(ASV %in% worle_strip_persistors, "Strip_associated","No_unique_association")))) %>%
  column_to_rownames("ASV") %>%
  as.matrix()


tax_table(worle.soil.phy) <- tax_association

####
# Testing relative abundance plot with full soil dataset
plot <- phylogeny_profile(worle.soil.phy, classification = 'ASV_Association', treatment = c("treatment"), merge = TRUE, relative_abundance = TRUE)
plot
data <- plot$data %>%
  group_by(treatment, soil_type, ASV_Association) %>%
  summarise_at(vars(-OTU, -Sample, -unique_id, -experiment, -matrix, -sample_day, -depth, -in_plot_location, -block, -strip, -manure_treatment, -plot), funs(mean(., na.rm=TRUE)))
data
write_csv(x = data, file = "data/manuretaxainsoil.csv")
