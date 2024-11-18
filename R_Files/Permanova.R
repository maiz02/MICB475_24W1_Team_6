## PERMANOVA ##
# To conduct a PERMANOVA in R, we need to load the "vegan" and "phyloseq" libraries.

library(vegan)
library(phyloseq)

load("R_Files/ms_rare.RData")

# Then, we can create a distance matrix with our metric of choice. Here are three examples:

dm_unifrac <- UniFrac(ms_rare, weighted=TRUE) # Weighted UniFrac
dm_braycurtis <- vegdist(t(otu_table(ms_rare)), method="bray") # Bray-curtis
dm_jaccard <- vegdist(t(otu_table(ms_rare)), method="jaccard") # Jaccard
dm_wunifrac <- 

# Then, we use the adonis2 function. Here, let's assume we have two predictor variables (predictor1 and predictor2) in the sample metadata table, dat:

dat <- as.data.frame(sample_data(phyloseq_object))
adonis2(dm_unifrac ~ predictor1*predictor2, data=dat)