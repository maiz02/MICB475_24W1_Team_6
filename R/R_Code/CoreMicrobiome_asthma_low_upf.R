BiocManager::install("microbiome")
install.packages("ggVennDiagram")

# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data #### need the unrarefied (final) phyloseq object
load("R_Code/upf_phyloseq_final_low.RData")

#### "core" microbiome ####

# Convert to relative abundance
upf_low_RA <- transform_sample_counts(upf_phyloseq_final_low, fun=function(x) x/sum(x))

# Filter dataset by asthma
upf_low_asthma <- subset_samples(upf_low_RA, `asthma`==1)
upf_low_noasthma <- subset_samples(upf_low_RA, `asthma`==0)

# Identify core microbiome members for each group
# Adjust prevalence and detection thresholds as necessary
upf_low_asthma_list <- core_members(upf_low_asthma, detection=0.001, prevalence = 0.50)
upf_low_noasthma_list <- core_members(upf_low_noasthma, detection=0.001, prevalence = 0.50)

# Combine lists into a single data structure
upf_low_asthma_list_full <- list(Asthma = upf_low_asthma_list, No_Asthma = upf_low_noasthma_list)

## MAC USERS BEWARE ##
install.packages("sf")
library("sf")
## MAC USERS BEWARE ENDS HERE##

# Create a Venn diagram using all the ASVs shared and unique to asthma and non-asthma

low_venn <- ggVennDiagram(
  x = upf_low_asthma_list_full,
  label_alpha = 1.0,
  label_size = 7,
  edge_size = 0.5
) +
  coord_cartesian(clip = "off") +
  coord_fixed(ratio = 0.75) + 
  labs(caption = "Low UPF") + 
  theme(
    plot.caption = element_text(family = "Helvetica", hjust = 0.5, size = 14, face = "bold", margin = margin(t = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    text = element_text(face = "bold", family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"))+
  scale_fill_gradient(low = "mistyrose", high = "salmon")



ggsave("R_files/core_microbiome/venn_asthma_low_upf.png", plot = low_venn, width = 12, height = 7, dpi = 300)

