# Load necessary libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load datasets for high UPF and low UPF
load("R_Code/upf_phyloseq_final_high.RData")
load("R_Code/upf_phyloseq_final_low.RData")

# Convert datasets to relative abundance
upf_high_RA <- transform_sample_counts(upf_phyloseq_final_high, fun = function(x) x / sum(x))
upf_low_RA <- transform_sample_counts(upf_phyloseq_final_low, fun = function(x) x / sum(x))

# Filter datasets by asthma status
upf_high_asthma <- subset_samples(upf_high_RA, `asthma` == 1)
upf_high_noasthma <- subset_samples(upf_high_RA, `asthma` == 0)
upf_low_asthma <- subset_samples(upf_low_RA, `asthma` == 1)
upf_low_noasthma <- subset_samples(upf_low_RA, `asthma` == 0)

# Identify core microbiome members for each group
# Adjust prevalence and detection thresholds as necessary
high_asthma_ASVs <- core_members(upf_high_asthma, detection = 0.001, prevalence = 0.10)
high_noasthma_ASVs <- core_members(upf_high_noasthma, detection = 0.001, prevalence = 0.10)
low_asthma_ASVs <- core_members(upf_low_asthma, detection = 0.001, prevalence = 0.10)
low_noasthma_ASVs <- core_members(upf_low_noasthma, detection = 0.001, prevalence = 0.10)

# Combine lists into a single data structure
core_ASVs_list <- list(
  `High UPF - Asthma` = high_asthma_ASVs,
  `High UPF - No Asthma` = high_noasthma_ASVs,
  `Low UPF - Asthma` = low_asthma_ASVs,
  `Low UPF - No Asthma` = low_noasthma_ASVs
)

# Create a Venn diagram with all four groups
combined_venn <- ggVennDiagram(
  x = core_ASVs_list,
  label_alpha = 1.0,
  label_size = 5,
  edge_size = 0.5
) + 
  coord_cartesian(clip = "off") + 
  coord_fixed(ratio = 0.75) # Adjust ratio as needed for aesthetics

# Save the combined Venn diagram
ggsave("R_files/venn_asthma_upf_combined.png", plot = combined_venn, width = 12, height = 7, dpi = 300)



