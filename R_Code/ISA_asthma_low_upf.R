# if you didn't install the indicspecies package, run the following
install.packages("indicspecies")

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("upf_phyloseq_final_low.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
low_upf_genus <- tax_glom(upf_phyloseq_final_low, "Genus", NArm = FALSE)
low_upf_genus_RA <- transform_sample_counts(low_upf_genus, fun=function(x) x/sum(x))

#ISA
isa_low_upf <- multipatt(t(otu_table(low_upf_genus_RA)), cluster = sample_data(low_upf_genus_RA)$`asthma`,control = how(nperm=20000))
summary(isa_low_upf)
taxtable <- tax_table(upf_phyloseq_final_low) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_low_upf$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

