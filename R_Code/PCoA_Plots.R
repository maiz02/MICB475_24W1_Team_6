library(phyloseq)
library(ape)
library(tidyverse)

load("R_Files/ms_rare.RData")
set.seed(500)

#### BRAY-CURTIS ####
# bray-curtis for asthma
bc_asthma_dm <- distance(ms_rare, method="bray")
pcoa_asthma_bc <- ordinate(ms_rare, method="PCoA", distance=bc_asthma_dm)
bc_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_bc, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
bc_asthma_pcoa

# updated bray-curtis for asthma
updated_bc_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_bc, shape = "upf_status", color = "asthma_yn") +
  labs(col = "Presence of asthma", pch = "UPF Status") + 
  stat_ellipse(aes(group = asthma_yn, color = asthma_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No asthma", "Asthma")) 
updated_bc_asthma_pcoa

# bray-curtis for allergies
bc_allergies_dm <- distance(ms_rare, method="bray")
pcoa_allergies_bc <- ordinate(ms_rare, method="PCoA", distance=bc_allergies_dm)
bc_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_bc, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
bc_allergies_pcoa

# updated bray-curtis for allergies
updated_bc_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_bc, shape = "upf_status", color = "allergies_yn") +
  labs(col = "Presence of allergies", pch = "UPF Status") + 
  stat_ellipse(aes(group = allergies_yn, color = allergies_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No allergies", "Allergies")) 
updated_bc_allergies_pcoa

# bray-curtis for both
bc_both_dm <- distance(ms_rare, method="bray")
pcoa_both_bc <- ordinate(ms_rare, method="PCoA", distance=bc_both_dm)
sample_data(ms_rare)$upf_status <- factor(
  sample_data(ms_rare)$upf_status,
  levels = c("upf_low", "upf_high"))
bc_both_pcoa <- plot_ordination(ms_rare, pcoa_both_bc,
                                color = "asthma_yn", 
                                shape = "allergies_yn") +
  geom_point(aes(size = upf_status), alpha = 0.7) +
  scale_size_manual(values = c("upf_low" = 1, "upf_high" = 4)) +
  labs(title = "Bray-Curtis PCoA Plot",
       color = "Asthma status",
       shape = "Allergies status",
       size = "UPF status") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2),
         shape = guide_legend(order = 3)) +
  theme_minimal() +
  theme(legend.position = "right")
bc_both_pcoa

#### JACCARD ####
# jaccard for asthma
jaccard_asthma_dm <- distance(ms_rare, method="jaccard")
pcoa_asthma_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_asthma_dm)
jaccard_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_jaccard, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
jaccard_asthma_pcoa

# updated jaccard for asthma
updated_jaccard_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_jaccard, shape = "upf_status", color = "asthma_yn") +
  labs(col = "Presence of asthma", pch = "UPF Status") + 
  stat_ellipse(aes(group = asthma_yn, color = asthma_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No asthma", "Asthma")) 
updated_jaccard_asthma_pcoa

# jaccard for allergies
jaccard_allergies_dm <- distance(ms_rare, method="jaccard")
pcoa_allergies_jaccard <- ordinate(ms_rare, method="PCoA", distance=jaccard_allergies_dm)
jaccard_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_jaccard, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
jaccard_allergies_pcoa

# updated jaccard for allergies
updated_jaccard_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_jaccard, shape = "upf_status", color = "allergies_yn") +
  labs(col = "Presence of allergies", pch = "UPF Status") + 
  stat_ellipse(aes(group = allergies_yn, color = allergies_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No allergies", "Allergies")) 
updated_jaccard_allergies_pcoa

# jaccard for both
jaccard_both_dm <- distance(ms_rare, method="jaccard")
pcoa_both_jaccard <- ordinate(ms_rare, method="PCoA", distance = jaccard_both_dm)
sample_data(ms_rare)$upf_status <- factor(
  sample_data(ms_rare)$upf_status,
  levels = c("upf_low", "upf_high"))
jaccard_both_pcoa <- plot_ordination(ms_rare, pcoa_both_jaccard,
                                color = "asthma_yn", 
                                shape = "allergies_yn") +
  geom_point(aes(size = upf_status), alpha = 0.7) +
  scale_size_manual(values = c("upf_low" = 1, "upf_high" = 4)) +
  labs(title = "Jaccard PCoA Plot",
       color = "Asthma status",
       shape = "Allergies status",
       size = "UPF status") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2),
         shape = guide_legend(order = 3)) +
  theme_minimal() +
  theme(legend.position = "right")
jaccard_both_pcoa

#### UNWEIGHTED UNIFRAC ####
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

# unweighted for both
uni_both_dm <- distance(ms_rare, method="unifrac")
pcoa_both_uni <- ordinate(ms_rare, method="PCoA", distance = uni_both_dm)
sample_data(ms_rare)$upf_status <- factor(
  sample_data(ms_rare)$upf_status,
  levels = c("upf_low", "upf_high"))
uni_both_pcoa <- plot_ordination(ms_rare, pcoa_both_uni,
                                     color = "asthma_yn", 
                                     shape = "allergies_yn") +
  geom_point(aes(size = upf_status), alpha = 0.7) +
  scale_size_manual(values = c("upf_low" = 1, "upf_high" = 4)) +
  labs(title = "Unweighted UniFrac PCoA Plot",
       color = "Asthma status",
       shape = "Allergies status",
       size = "UPF status") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2),
         shape = guide_legend(order = 3)) +
  theme_minimal() +
  theme(legend.position = "right")
uni_both_pcoa

#### WEIGHTED UNIFRAC ####
# weighed unifrac for asthma
wuni_asthma_dm <- distance(ms_rare, method="wunifrac")
pcoa_asthma_wuni <- ordinate(ms_rare, method="PCoA", distance=wuni_asthma_dm)
wuni_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_wuni, color = "upf_status", shape = "asthma_yn") +
  labs(pch = "Presence of asthma", col="UPF Status")
wuni_asthma_pcoa

# updated wuni for asthma
updated_wuni_asthma_pcoa <- plot_ordination(ms_rare, pcoa_asthma_wuni, shape = "upf_status", color = "asthma_yn") +
  labs(col = "Presence of asthma", pch = "UPF Status") + 
  stat_ellipse(aes(group = asthma_yn, color = asthma_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No asthma", "Asthma")) +
  xlim(-0.1, 0.05) + ylim(-0.050, 0.015)
updated_wuni_asthma_pcoa

# weighted unifrac for allergies
wuni_allergies_dm <- distance(ms_rare, method="wunifrac")
pcoa_allergies_wuni <- ordinate(ms_rare, method="PCoA", distance=wuni_allergies_dm)
wuni_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_wuni, color = "upf_status", shape = "allergies_yn") +
  labs(pch = "Presence of allergies", col="UPF Status")
wuni_allergies_pcoa

# updated wuni for allergies
updated_wuni_allergies_pcoa <- plot_ordination(ms_rare, pcoa_allergies_wuni, shape = "upf_status", color = "allergies_yn") +
  labs(col = "Presence of allergies", pch = "UPF Status") + 
  stat_ellipse(aes(group = allergies_yn, color = allergies_yn), level = 0.95, linewidth = 0.8, show.legend = TRUE) + 
  guides(shape = guide_legend(override.aes = list(linetype = 0))) + 
  scale_shape_manual(values = c(16, 17),
                     labels = c("High UPF", "Low UPF")) + 
  scale_color_manual(values = c("#F8766D", "#00BFC4"),  
                     labels = c("No allergies", "Allergies")) +
  xlim(-0.1, 0.05) + ylim(-0.050, 0.015)
updated_wuni_allergies_pcoa

# weighted unifrac for both
wuni_both_dm <- distance(ms_rare, method="wunifrac")
pcoa_both_wuni <- ordinate(ms_rare, method="PCoA", distance = wuni_both_dm)
sample_data(ms_rare)$upf_status <- factor(
  sample_data(ms_rare)$upf_status,
  levels = c("upf_low", "upf_high"))
wuni_both_pcoa <- plot_ordination(ms_rare, pcoa_both_wuni,
                                  color = "asthma_yn", 
                                  shape = "allergies_yn") +
  geom_point(aes(size = upf_status), alpha = 0.7) +
  scale_size_manual(values = c("upf_low" = 1, "upf_high" = 4)) +
  labs(title = "Weighted UniFrac PCoA Plot",
       color = "Asthma status",
       shape = "Allergies status",
       size = "UPF status") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2),
         shape = guide_legend(order = 3)) +
  theme_minimal() +
  theme(legend.position = "right")
wuni_both_pcoa

# updated wuni both
sample_data(ms_rare)$grouping <- paste(sample_data(ms_rare)$upf_status,
                                       sample_data(ms_rare)$asthma_yn,
                                       sample_data(ms_rare)$allergies_yn,
                                       sep = ", ")

updated_wuni_both_pcoa <- plot_ordination(ms_rare, pcoa_both_wuni,
                                  color = "grouping") +
  geom_point(alpha = 0.7) +
  stat_ellipse(aes(color = grouping), level = 0.95, linewidth = 0.8) +
  scale_color_manual(
    values = c("upf_high, no, no" = "#E41A1C", 
                "upf_high, no, yes" = "#377EB8",
                "upf_high, yes, no" = "#4DAF4A",
                "upf_high, yes, yes" = "#FF7F00",
                "upf_low, no, no" = "#F781BF",
                "upf_low, no, yes" = "#984EA3",
                "upf_low, yes, no" = "#FFBF00",
                "upf_low, yes, yes" = "#00BFC4"),
    labels = c("upf_high, no, no" = "High UPF, No asthma, No allergies", 
               "upf_high, no, yes" = "High UPF, No asthma, Allergies", 
               "upf_high, yes, no" = "High UPF, Asthma, No allergies",
               "upf_high, yes, yes" = "High UPF, Asthma, Allergies",
               "upf_low, no, no" = "Low UPF, No asthma, No allergies",
               "upf_low, no, yes" = "Low UPF, No asthma, Allergies",
               "upf_low, yes, no" = "Low UPF, Asthma, No allergies",
               "upf_low, yes, yes" = "Low UPF, Asthma, Allergies")) +
  labs(title = "Weighted UniFrac PCoA Plot",
       color = "Groups (UPF, Asthma, Allergies)") +
  guides(color = guide_legend(order = 1)) +
  theme(legend.position = "right") +
  xlim(-0.060, 0.05) +
  ylim(-0.050, 0.055)
updated_wuni_both_pcoa

#### COUNTRIES ####
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

#### SAVING PLOTS ####
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

ggsave("R_files/bc_both_pcoa.png",
       bc_both_pcoa,
       height=4, width=5)

ggsave("R_files/jaccard_both_pcoa.png",
       jaccard_both_pcoa,
       height=4, width=5)

ggsave("R_files/uni_both_pcoa.png",
       uni_both_pcoa,
       height=4, width=5)

ggsave("R_files/wuni_both_pcoa.png",
       wuni_both_pcoa,
       height=4, width=5)

ggsave("R_files/updated_wuni_both_pcoa.png",
       updated_wuni_both_pcoa,
       height=4, width=5)

ggsave("R_files/updated_bc_asthma_pcoa.png",
       updated_bc_asthma_pcoa,
       height=4, width=5)

ggsave("R_files/updated_bc_allergies_pcoa.png",
       updated_bc_allergies_pcoa,
       height=4, width=5)

ggsave("R_files/updated_jaccard_asthma_pcoa.png",
       updated_jaccard_asthma_pcoa,
       height=4, width=5)

ggsave("R_files/updated_jaccard_allergies_pcoa.png",
       updated_jaccard_allergies_pcoa,
       height=4, width=5)

ggsave("R_files/updated_wuni_asthma_pcoa.png",
       updated_wuni_asthma_pcoa,
       height=4, width=5)

ggsave("R_files/updated_wuni_allergies_pcoa.png",
       updated_wuni_allergies_pcoa,
       height=4, width=5)