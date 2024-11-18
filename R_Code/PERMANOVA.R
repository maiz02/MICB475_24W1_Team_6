library(vegan)
library(phyloseq)

dm_unifrac <- UniFrac(ms_rare, weighted=TRUE) # Weighted UniFrac
dm_braycurtis <- vegdist(t(otu_table(ms_rare)), method="bray") # Bray-curtis
dm_jaccard <- vegdist(t(otu_table(ms_rare)), method="jaccard") # Jaccard

df <- as.data.frame(as.matrix(sample_data(ms_rare)))

#allergies
adonis2(dm_unifrac ~ upf_status * allergies, data =df)
adonis2(dm_braycurtis ~ upf_status*allergies, data =df)
adonis2(dm_jaccard ~ upf_status*allergies, data =df)

#asthma
adonis2(dm_unifrac ~ upf_status * asthma, data =df)
adonis2(dm_braycurtis ~ upf_status*asthma, data =df)
adonis2(dm_jaccard ~ upf_status*asthma, data =df)
