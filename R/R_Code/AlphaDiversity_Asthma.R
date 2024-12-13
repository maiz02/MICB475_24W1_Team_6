install.packages("patchwork")

library(phyloseq)
library(tidyverse)
library(picante)
library(patchwork)
library(dplyr)
library(ape)

load("R_Code/upf_phyloseq_rare_high.RData")
load("R_Code/upf_phyloseq_rare_low.RData")

#### HIGH UPF ####
# Observed Features and Shannon Diversity
gg_richness_high <- plot_richness(upf_phyloseq_rare_high, x = "upf_asthma", 
                                  measures = c("Observed","Shannon")) +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  +
  geom_boxplot(aes(fill = asthma_yn), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("no" = "steelblue", "yes" = "salmon"),
                    name = "Asthma status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none")
gg_richness_high
ggsave(filename = "R_Files/gg_richness_high.png", 
       gg_richness_high,
       height=4, width=6)

# Faith's PD
phylo_dist_high <- pd(t(otu_table(upf_phyloseq_rare_high)), phy_tree(upf_phyloseq_rare_high),
                 include.root=F) 

sample_data(upf_phyloseq_rare_high)$PD <- phylo_dist_high$PD

pd_plot_high <- ggplot(sample_data(upf_phyloseq_rare_high), aes(upf_asthma, PD)) + 
  geom_point(position = position_dodge(width = 1), size = 1.5, alpha = 1) +
  geom_boxplot(aes(fill = asthma_yn), alpha = 0.7, outlier.shape = NA) + 
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma")) +  # Label adjustment
  ylab("Phylogenetic Diversity") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c("no" = "steelblue", "yes" = "salmon"),
                    name = "Asthma status") +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  scale_colour_manual(name="", values=c(Com = "lightgrey") )

pd_plot_high<- pd_plot_high+facet_grid(. ~ "Faith's PD") +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=8, colour="black"))
pd_plot_high


ggsave(filename = "R_Files/pd_plot_high.png", 
       pd_plot_high,
       height=4, width=3)


combined_plot <- gg_richness_high + pd_plot_high + 
  plot_layout(ncol = 2, widths = c(2, 1))
combined_plot
ggsave(filename = "R_Files/combined_high_upf.png", 
       combined_plot, height = 4, width = 9)

## Stats
# Observed and Shannon
alphadiv_high <- estimate_richness(upf_phyloseq_rare_high)
samp_dat_high <- sample_data(upf_phyloseq_rare_high)
samp_dat_high_wdiv <- data.frame(samp_dat_high, alphadiv_high)

wilcox.test(Observed ~ asthma, data=samp_dat_high_wdiv, exact = FALSE)
wilcox.test(Shannon ~ asthma, data=samp_dat_high_wdiv, exact = FALSE)

# Faith's PD
phylo_dist_high <- pd(t(otu_table(upf_phyloseq_rare_high)), phy_tree(upf_phyloseq_rare_high),
                 include.root=F) 
sample_data(upf_phyloseq_rare_high)$PD <- phylo_dist_high$PD
pd_values <- sample_data(upf_phyloseq_rare_high)$PD
asthma_high <- sample_data(upf_phyloseq_rare_high)$asthma
wilcox.test(pd_values ~ asthma_high)

#### LOW UPF ####
# Observed Features and Shannon Diversity
gg_richness_low <- plot_richness(upf_phyloseq_rare_low, x = "upf_asthma", 
                                  measures = c("Observed","Shannon")) +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  +
  geom_boxplot(aes(fill = asthma_yn), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("no" = "steelblue", "yes" = "salmon"),
                    name = "Asthma status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none")
gg_richness_low




ggsave(filename = "R_Files/gg_richness_low.png", 
       gg_richness_low,
       height=4, width=6)


# Faith's PD
phylo_dist_low <- pd(t(otu_table(upf_phyloseq_rare_low)), phy_tree(upf_phyloseq_rare_low),
                      include.root=F) 

sample_data(upf_phyloseq_rare_low)$PD <- phylo_dist_low$PD

pd_plot_low <- ggplot(sample_data(upf_phyloseq_rare_low), aes(upf_asthma, PD)) + 
  geom_point(position = position_dodge(width = 1), size = 1.5, alpha = 1) +
  geom_boxplot(aes(fill = asthma_yn), alpha = 0.7, outlier.shape = NA) + 
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma")) +  # Label adjustment
  ylab("Phylogenetic Diversity") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c("no" = "steelblue", "yes" = "salmon"),
                    name = "Asthma status") +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  scale_colour_manual(name="", values=c(Com = "lightgrey") )

pd_plot_low<- pd_plot_low+facet_grid(. ~ "Faith's PD") +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=8, colour="black"))
pd_plot_low


ggsave(filename = "R_Files/pd_plot_low.png", 
       pd_plot_low,
       height=4, width=3)


combined_plot_low <- gg_richness_low + pd_plot_low + 
  plot_layout(ncol = 2, widths = c(2, 1))
combined_plot_low
ggsave(filename = "R_Files/combined_low_upf.png", 
       combined_plot, height = 4, width = 9)


## Stats
# Observed and Shannon
alphadiv_low <- estimate_richness(upf_phyloseq_rare_low)
samp_dat_low <- sample_data(upf_phyloseq_rare_low)
samp_dat_low_wdiv <- data.frame(samp_dat_low, alphadiv_low)

wilcox.test(Observed ~ asthma, data=samp_dat_low_wdiv, exact = FALSE)
wilcox.test(Shannon ~ asthma, data=samp_dat_low_wdiv, exact = FALSE)

# Faith's PD
phylo_dist_low <- pd(t(otu_table(upf_phyloseq_rare_low)), phy_tree(upf_phyloseq_rare_low), include.root = F)
sample_data(upf_phyloseq_rare_low)$PD <- phylo_dist_low$PD
pd_values_low <- sample_data(upf_phyloseq_rare_low)$PD
asthma_low <- sample_data(upf_phyloseq_rare_low)$asthma
wilcox.test(pd_values_low ~ asthma_low)


#### COMBINED UPF ####
# Load data
load("R_Code/upf_phyloseq_rare.RData")
# ^ this loads the upf rare high for some reason idk how to fix it!

# Plot Observed and Shannon, high vs low UPF
gg_richness_high <- plot_richness(upf_phyloseq_rare_high, x = "upf_status", 
                                  measures = c("Observed","Shannon")) +
  xlab(" ") +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))  
gg_richness_high

# Stats for combined UPF, Observed and Shannon
alphadiv <- estimate_richness(upf_phyloseq_rare)
samp_dat <- sample_data(upf_phyloseq_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

wilcox.test(Observed ~ upf_status, data=samp_dat_high_wdiv, exact = FALSE)
wilcox.test(Shannon ~ upf_status, data=samp_dat_high_wdiv, exact = FALSE)

# Plot Faith's PD, high vs low UPF
phylo_dist_high <- pd(t(otu_table(upf_phyloseq_rare_high)), phy_tree(upf_phyloseq_rare_high),
                      include.root=F) 

sample_data(upf_phyloseq_rare_high)$PD <- phylo_dist_high$PD

pd_plot_high <- ggplot(sample_data(upf_phyloseq_rare_high), aes(upf_status, PD)) + 
  geom_point(position = position_dodge(width = 1), size = 1.5, alpha = 1) +
  geom_boxplot() +
  xlab(" ") +
  ylab("Phylogenetic Diversity") +
  theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11))+
  scale_colour_manual(name="", values=c(Com = "lightgrey") )

pd_plot_high<- pd_plot_high+facet_grid(. ~ "Faith's PD") +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size=8, colour="black"))
pd_plot_high


