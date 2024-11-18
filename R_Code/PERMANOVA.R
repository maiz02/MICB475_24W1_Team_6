library(vegan)
library(phyloseq)

dm_unifrac <- UniFrac(ms_rare, weighted=TRUE) # Weighted UniFrac
dm_braycurtis <- vegdist(t(otu_table(ms_rare)), method="bray") # Bray-curtis
dm_jaccard <- vegdist(t(otu_table(ms_rare)), method="jaccard") # Jaccard

df <- as.data.frame(as.matrix(sample_data(ms_rare)))

#allergies
permanova_wuni_allergies <- adonis2(dm_unifrac ~ upf_status * allergies, data =df)
permanova_bray_allergies <- adonis2(dm_braycurtis ~ upf_status*allergies, data =df)
permanova_jacc_allergies <- adonis2(dm_jaccard ~ upf_status*allergies, data =df)

#asthma
permanova_wuni_asthma <- adonis2(dm_unifrac ~ upf_status * asthma, data =df)
permanova_bray_asthma <- adonis2(dm_braycurtis ~ upf_status*asthma, data =df)
permanova_jacc_asthma <- adonis2(dm_jaccard ~ upf_status*asthma, data =df)

# Save the results as a CSV file
write.csv(permanova_wuni_allergies, "R_Files/permanova_wuni_allergies.csv", row.names = TRUE)
write.csv(permanova_bray_allergies, "R_Files/permanova_bray_allergies.csv", row.names = TRUE)
write.csv(permanova_jacc_allergies, "R_Files/permanova_jacc_allergies.csv", row.names = TRUE)
write.csv(permanova_wuni_asthma, "R_Files/permanova_wuni_asthma.csv", row.names = TRUE)
write.csv(permanova_bray_asthma, "R_Files/permanova_bray_asthma.csv", row.names = TRUE)
write.csv(permanova_jacc_asthma, "R_Files/permanova_wuni_asthma.csv", row.names = TRUE)

