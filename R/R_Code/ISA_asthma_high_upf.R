
install.packages("indicspecies")

# load packages
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("upf_phyloseq_final_high.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
high_upf_genus <- tax_glom(upf_phyloseq_final_high, "Genus", NArm = FALSE)
high_upf_genus_RA <- transform_sample_counts(high_upf_genus, fun=function(x) x/sum(x))

#ISA
isa_high_upf <- multipatt(t(otu_table(high_upf_genus_RA)), cluster = sample_data(high_upf_genus_RA)$`asthma`,control = how(nperm=20000))
summary(isa_high_upf)
taxtable <- tax_table(upf_phyloseq_final_high) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_high_upf$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View() 

# Results: 1 species associated to Group 1, stat value of 0.3 --> less than 0.8, not very indicative