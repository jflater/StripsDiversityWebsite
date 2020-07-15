library(DECIPHER)
library(phangorn)
library(dada2)
library(phyloseq)
library(doParallel)

setwd("~/Documents/StripsDiversityWebsite/")
strips <- readRDS("data/all_samples_merged_STRIPS.RDS") 

Allbase <- subset_samples(strips, sample_day %in% "Baseline" & !soil_type %in% "border" & depth %in% "1") %>%
  filter_taxa(function(x) sum(x >= 3) > (0.2*length(x)), T) 

rm(strips)
taxa_names(Allbase) <- paste0("ASV", seq(ntaxa(Allbase)))

seqtab <- refseq(Allbase)
seqs <- getSequences(seqtab)
names(seqs) <- taxa_names(Allbase) # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, processors = NULL)

print("Detecting # of cores")
detectCores(all.tests = TRUE)
cl <- makeCluster(detectCores(all.tests = TRUE))
print("Register parallel")
registerDoParallel(cores = 6)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
print("Now doing function NJ()")
treeNJ <- NJ(dm) # Note, tip order != sequence order
print("Now doing function pml()")
fit = pml(treeNJ, data = phang.align, )
print("Saving fit")
saveRDS(object = fit, "data/fittree.RDS")

## negative edges length changed to 0!
print("Updating fit")
fitGTR <- update(fit, k = 4, inv = 0.2)

# July 2nd: breaking here
print("Optimizing model")
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


print("Stop parallel")
stopCluster(cl)
rm(cl)
registerDoSEQ() 

print("Saving tree object and just tree as RDS")
saveRDS(object = fitGTR, "data/d1mintax5treemodel.RDS") # hasn't been run as of July 2nd
saveRDS(object = fitGTR$tree, "data/d1mintax5tree.RDS")
