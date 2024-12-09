# if you didn't install microbiome or ggVennDiagram packages, run the following code
# these can take quite a few minutes to install, give it 5-15 minutes for microbiome
# ggVennDiagram can take up to 10-25 minutes to install
BiocManager::install("microbiome")
install.packages("ggVennDiagram")

#!/usr/bin/env Rscript
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

# What ASVs are found in more than 70% of samples in each antibiotic usage category?
# trying changing the prevalence to see what happens
upf_low_asthma_ASVs <- core_members(upf_low_asthma, detection=0, prevalence = 0.7)
upf_low_noasthma_ASVs <- core_members(upf_low_noasthma, detection=0, prevalence = 0.7)

# What are these ASVs? you can code it in two different ways to see the same things
tax_table(prune_taxa(upf_low_asthma_ASVs,upf_phyloseq_final_low))
tax_table(prune_taxa(upf_low_noasthma_ASVs,upf_phyloseq_final_low))

# can plot those ASVs' relative abundance
prune_taxa(upf_low_asthma_ASVs,upf_low_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`asthma`, scales ="free")

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
upf_low_asthma_list <- core_members(upf_low_asthma, detection=0.001, prevalence = 0.10)
upf_low_noasthma_list <- core_members(upf_low_noasthma, detection=0.001, prevalence = 0.10)

upf_low_asthma_list_full <- list(Asthma = upf_low_asthma_list, No_Asthma = upf_low_noasthma_list)

## MAC USERS BEWARE ##
# If you are working on a new Mac running one of their new OS, you may get an error when
# you try to run the ggVennDiagram function that says something along the lines of a 
# a package called 'sf' is not found so you need to install said package and a couple more
# packages to ensure that ggVennDiagram is going to work on your computer
install.packages("sf")
# BUT when it asks if you want to install from sources, answer NO!
# load the package before running the next command, not good practice to load midway but this
# for mac users only
library("sf")
## MAC USERS BEWARE ENDS HERE##

# Create a Venn diagram using all the ASVs shared and unique to asthma and non-asthma
low_venn <- ggVennDiagram(x = upf_low_asthma_list_full, label_alpha = 1.0,
                           label_size = 5,        
                           edge_size = 0.5) + coord_cartesian(clip = "off") + coord_fixed(ratio = 0.55)

ggsave("R_files/venn_asthma_low_upf.png", plot = low_venn, width = 12, height = 7, dpi = 300)

