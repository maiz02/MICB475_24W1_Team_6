library(phyloseq)
library(ape)
library(tidyverse)

load("R_Files/ms_rare.RData")

# all countries
bc_asthma_countries_dm <- distance(ms_rare, method="bray")
pcoa_asthma_countries_bc <- ordinate(ms_rare, method="PCoA", distance=bc_asthma_countries_dm)
bc_asthma_countries_pcoa <- plot_ordination(ms_rare, pcoa_asthma_countries_bc, color = "country", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="Country")
bc_asthma_countries_pcoa

bc_allergies_countries_dm <- distance(ms_rare, method="bray")
pcoa_allergies_countries_bc <- ordinate(ms_rare, method="PCoA", distance=bc_allergies_countries_dm)
bc_allergies_countries_pcoa <- plot_ordination(ms_rare, pcoa_allergies_countries_bc, color = "country", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="Country")
bc_allergies_countries_pcoa

jaccard_asthma_countries_dm <- distance(ms_rare, method="jaccard")
pcoa_asthma_countries_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_asthma_countries_dm)
jaccard_asthma_countries_pcoa <- plot_ordination(ms_rare, pcoa_asthma_countries_jaccard, color = "country", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="Country")
jaccard_asthma_countries_pcoa

jaccard_allergies_countries_dm <- distance(ms_rare, method="jaccard")
pcoa_allergies_countries_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_allergies_countries_dm)
jaccard_allergies_countries_pcoa <- plot_ordination(ms_rare, pcoa_allergies_countries_jaccard, color = "country", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="Country")
jaccard_allergies_countries_pcoa

# bray-curtis for asthma
bc_asthma_dm <- distance(ms_rare, method="bray")
pcoa_asthma_bc <- ordinate(ms_rare, method="PCoA", distance=bc_asthma_dm)
bc_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_bc, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
bc_asthma_pcoa

# bray-curtis for allergies
bc_allergies_dm <- distance(ms_rare, method="bray")
pcoa_allergies_bc <- ordinate(ms_rare, method="PCoA", distance=bc_allergies_dm)
bc_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_bc, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
bc_allergies_pcoa

# bray-curtis for both - NOT SURE HOW TO DO THIS ONE
sample_data(ms_rare)$asthma_allergies_yn <- interaction(sample_data(ms_rare)$asthma_yn, sample_data(ms_rare)$allergies_yn)
bc_both_dm <- distance(ms_rare, method="bray")
pcoa_both_bc <- ordinate(ms_rare, method="PCoA", distance=bc_both_dm)
bc_both_pcoa <- plot_ordination(ms_rare, pcoa_both_bc, 
                                     color = "upf_status", 
                                     shape = "asthma_allergies_yn") +
  labs(pch = "Asthma and Allergies Status", col = "UPF Status") +
  scale_shape_manual(values = c(16, 17, 18, 19))
bc_both_pcoa

# jaccard for asthma
jaccard_asthma_dm <- distance(ms_rare, method="jaccard")
pcoa_asthma_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_asthma_dm)
jaccard_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_jaccard, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
jaccard_asthma_pcoa

# jaccard for allergies
jaccard_allergies_dm <- distance(ms_rare, method="jaccard")
pcoa_allergies_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_allergies_dm)
jaccard_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_jaccard, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
jaccard_allergies_pcoa

# unweighted unifrac for asthma
uni_asthma_dm <- distance(ms_rare, method="unifrac")
pcoa_asthma_uni <- ordinate(ms_rare, method="PCoA", distance=uni_asthma_dm)
uni_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_uni, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
uni_asthma_pcoa

# unweighted unifrac for allergies
uni_allergies_dm <- distance(ms_rare, method="unifrac")
pcoa_allergies_uni <- ordinate(ms_rare, method="PCoA", distance=uni_allergies_dm)
uni_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_uni, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
uni_allergies_pcoa

# weighed unifrac for asthma
wuni_asthma_dm <- distance(ms_rare, method="wunifrac")
pcoa_asthma_wuni <- ordinate(ms_rare, method="PCoA", distance=wuni_asthma_dm)
wuni_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_wuni, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
wuni_asthma_pcoa

# weighted unifrac for allergies
wuni_allergies_dm <- distance(ms_rare, method="wunifrac")
pcoa_allergies_wuni <- ordinate(ms_rare, method="PCoA", distance=wuni_allergies_dm)
wuni_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_wuni, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
wuni_allergies_pcoa

# saving plots
ggsave("R_Files/bc_asthma_pcoa.png",
       bc_asthma_pcoa,
       height=4, width=5)

ggsave("R_Files/bc_allergies_pcoa.png",
       bc_allergies_pcoa,
       height=4, width=5)

ggsave("R_Files/jaccard_asthma_pcoa.png",
       jaccard_asthma_pcoa,
       height=4, width=5)

ggsave("R_Files/jaccard_allergies_pcoa.png",
       jaccard_allergies_pcoa,
       height=4, width=5)

ggsave("R_Files/uni_asthma_pcoa.png",
       uni_asthma_pcoa,
       height=4, width=5)

ggsave("R_Files/uni_allergies_pcoa.png",
       uni_allergies_pcoa,
       height=4, width=5)

ggsave("R_Files/wuni_asthma_pcoa.png",
       wuni_asthma_pcoa,
       height=4, width=5)

ggsave("R_Files/wuni_allergies_pcoa.png",
       wuni_allergies_pcoa,
       height=4, width=5)
