install.packages("patchwork")

library(phyloseq)
library(tidyverse)
library(picante)
library(patchwork)
library(dplyr)
library(ape)

load("R_Code/upf_phyloseq_final.RData")
load("R_Code/upf_phyloseq_rare.RData")

obs_asthma <- wilcox.test(Observed ~ asthma, data = samp_dat_wdiv, exact = FALSE)$p.value
obs_upf <- wilcox.test(Observed ~ upf_status, data = samp_dat_wdiv, exact = FALSE)$p.value
shan_asthma <- wilcox.test(Shannon ~ asthma, data = samp_dat_wdiv, exact = FALSE)$p.value
shan_upf <- wilcox.test(Shannon ~ upf_status, data = samp_dat_wdiv, exact = FALSE)$p.value

# Format p-values for display
format_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", signif(p, 3)))
  }
}

obs_asthma_label <- format_p(obs_asthma)
obs_upf_label <- format_p(obs_upf)
shan_asthma_label <- format_p(shan_asthma)
shan_upf_label <- format_p(shan_upf)


# Alpha Diversity Plot (gg_richness)
gg_richness <- plot_richness(upf_phyloseq_final, x = "upf_asthma", 
                             measures = c("Observed", "Shannon")) +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma","No Asthma", "Asthma")) + 
  geom_boxplot(aes(fill = upf_status), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none") 
gg_richness

obs_annotations <- data.frame(
  x = c(2.5, 2.5),
  y = c(370, 400), 
  label = c(obs_asthma_label, obs_upf_label),
  variable = "Observed"  
)

shan_annotations <- data.frame(
  x = c(2.5, 2.5),
  y = c(5.5, 6), 
  label = c(shan_asthma_label, shan_upf_label),
  variable = "Shannon"  
)

gg_richness <- gg_richness +
  geom_text(
    data = obs_annotations, 
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE 
  ) +
  geom_text(
    data = shan_annotations, 
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE  
  )

gg_richness


ggsave(filename = "R_Files/AlphaDiversity_Plots/gg_richness_high_low_upf.png", 
       gg_richness, height = 4, width = 6)


# Faith's PD Plot (pd_plot)
phylo_dist <- pd(t(otu_table(upf_phyloseq_rare)), phy_tree(upf_phyloseq_rare), include.root = FALSE)
sample_data(upf_phyloseq_rare)$PD <- phylo_dist$PD

pd_plot <- ggplot(sample_data(upf_phyloseq_rare), aes(x = upf_asthma, y = PD, fill = upf_status)) + 
  geom_point(position = position_dodge(width = 0.8), size = 1.5, alpha = 1) + 
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma","No Asthma", "Asthma")) +  
  ylab("Phylogenetic Diversity") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11), 
        axis.title.y = element_text(size = 11))

pd_plot <- pd_plot + 
  facet_grid(. ~ "Faith's PD") +  
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 8, colour = "black"))

pd_plot 

ggsave(filename = "R_Files/AlphaDiversity_Plots/pd_plot_high_low_upf.png", 
       pd_plot, height = 4, width = 3)



combined_plot <- gg_richness + pd_plot + 
  plot_layout(ncol = 2, widths = c(2, 1))  
combined_plot

ggsave(filename = "R_Files/AlphaDiversity_Plots/combined_high_low_upf.png", 
       combined_plot, height = 4, width = 10)




## Stats
# Observed and Shannon
samp_dat <- sample_data(upf_phyloseq_rare)
alphadiv <- estimate_richness(upf_phyloseq_rare)
samp_dat_wdiv <- cbind(as.data.frame(samp_dat), alphadiv)
wilcox.test(Observed ~ asthma, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Observed ~ upf_status, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ asthma, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ upf_status, data=samp_dat_wdiv, exact = FALSE)

# Faith's PD
phylo_dist <- pd(t(otu_table(upf_phyloseq_rare)), phy_tree(upf_phyloseq_rare),
                 include.root=F) 
sample_data(upf_phyloseq_rare)$PD <- phylo_dist$PD
pd_values <- sample_data(upf_phyloseq_rare)$PD
asthma <- sample_data(upf_phyloseq_rare)$asthma
wilcox.test(pd_values ~ asthma)
wilcox.test(pd_values ~ upf_status)