library(phyloseq)
library(tidyverse)
library(picante)

load("R_Code/upf_phyloseq_rare_high.RData")
load("R_Code/upf_phyloseq_rare_low.RData")

#### HIGH UPF ####
# Observed Features and Shannon Diversity
gg_richness_high <- plot_richness(upf_phyloseq_rare_high, x = "upf_asthma", 
                             measures = c("Observed","Shannon")) +
  xlab("High UPF and Asthma") +
  geom_boxplot()
gg_richness_high

ggsave(filename = "R_Files/gg_richness_high.png", 
       gg_richness_high,
       height=4, width=6)

# Faith's PD
phylo_dist_high <- pd(t(otu_table(upf_phyloseq_rare_high)), phy_tree(upf_phyloseq_rare_high),
                 include.root=F) 

sample_data(upf_phyloseq_rare_high)$PD <- phylo_dist_high$PD

pd_plot_high <- ggplot(sample_data(upf_phyloseq_rare_high), aes(upf_asthma, PD)) + 
  geom_boxplot() +
  xlab("High UPF and Asthma") +
  ylab("Phylogenetic Diversity")
pd_plot_high

ggsave(filename = "R_Files/pd_plot_high.png", 
       pd_plot_high,
       height=4, width=6)

#### LOW UPF ####
# Observed Features and Shannon Diversity
gg_richness_low <- plot_richness(upf_phyloseq_rare_low, x = "upf_asthma", 
                             measures = c("Observed","Shannon")) +
  xlab("Low UPF and Asthma") +
  geom_boxplot()
gg_richness_low

ggsave(filename = "R_Files/gg_richness_low.png", 
       gg_richness_low,
       height=4, width=6)

# Faith's PD
phylo_dist_low <- pd(t(otu_table(upf_phyloseq_rare_low)), phy_tree(upf_phyloseq_rare_low),
                 include.root=F) 

sample_data(upf_phyloseq_rare_low)$PD <- phylo_dist_low$PD

pd_plot_low <- ggplot(sample_data(upf_phyloseq_rare_low), aes(upf_asthma, PD)) + 
  geom_boxplot() +
  xlab("Low UPF and Asthma") +
  ylab("Phylogenetic Diversity")
pd_plot_low

ggsave(filename = "R_Files/pd_plot_low.png", 
       pd_plot_low,
       height=4, width=6)