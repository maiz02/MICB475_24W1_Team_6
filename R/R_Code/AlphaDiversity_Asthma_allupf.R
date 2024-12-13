install.packages("patchwork")

library(phyloseq)
library(tidyverse)
library(picante)
library(patchwork)
library(dplyr)
library(ape)

load("R_Code/upf_phyloseq_final.RData")
load("R_Code/upf_phyloseq_rare.RData")


# Alpha Diversity Plot (gg_richness)
gg_richness_shannon <- plot_richness(upf_phyloseq_final, x = "upf_asthma", 
                                     measures = "Shannon") +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma","No Asthma", "Asthma")) + 
  geom_boxplot(aes(fill = upf_status), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none") +
  
  # Annotating only for Shannon plot
  annotate("text", x = 2.5, y = 6.3, label = "**", size = 5, hjust = 0.5) +
  annotate("text", x = 3, y = 6.8, label = "*", size = 5, hjust = 0.5) +
  annotate("text", x = 2.5, y = 7.3, label = "*", size = 5, hjust = 0.5) +
  annotate("text", x = 1.5, y = 5.6, label = "*", size = 5, hjust = 0.5,color = "red") +
  annotate("text", x = 3.5, y = 5.6, label = "**", size = 5, hjust = 0.5,color = "red") +
  
  # Adding branches (lines) for Shannon plot
  geom_segment(aes(x = 1, xend = 4, y = 6.2, yend = 6.2), color = "black", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 6.1, yend = 6.2), color = "black", size = 0.7) +
  geom_segment(aes(x = 4, xend = 4, y = 6.1, yend = 6.2), color = "black", size = 0.7) +
  
geom_segment(aes(x = 2, xend = 4, y = 6.7, yend = 6.7), color = "black", size = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = 6.6, yend = 6.7), color = "black", size = 0.7) +
  geom_segment(aes(x = 4, xend = 4, y = 6.6, yend = 6.7), color = "black", size = 0.7) +
  
geom_segment(aes(x = 2, xend = 3, y = 7.2, yend = 7.2), color = "black", size = 0.7)+
  geom_segment(aes(x = 2, xend = 2, y =7.1, yend = 7.2), color = "black", size = 0.7)+
  geom_segment(aes(x = 3, xend = 3, y = 7.1, yend = 7.2), color = "black", size = 0.7) +
  
  geom_segment(aes(x = 1, xend = 2, y = 5.5, yend = 5.5), color = "red", size = 0.7)+
  geom_segment(aes(x = 1, xend = 1, y = 5.4, yend = 5.5), color = "red", size = 0.7)+
  geom_segment(aes(x = 2, xend = 2, y = 5.4, yend = 5.5), color = "red", size = 0.7) +
  
  geom_segment(aes(x = 3, xend = 4, y = 5.5, yend = 5.5), color = "red", size = 0.7)+
  geom_segment(aes(x = 4, xend = 4, y = 5.4, yend = 5.5), color = "red", size = 0.7)+
  geom_segment(aes(x = 3, xend = 3, y = 5.4, yend = 5.5), color = "red", size = 0.7) 

gg_richness_shannon
# Save Shannon plot
ggsave(filename = "R_Files/AlphaDiversity_Plots/gg_richness_shannon_upf.png", 
       gg_richness_shannon, height = 4, width = 6)


# Observed Richness Plot
gg_richness_observed <- plot_richness(upf_phyloseq_final, x = "upf_asthma", 
                                      measures = "Observed") +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma","No Asthma", "Asthma")) + 
  geom_boxplot(aes(fill = upf_status), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none") +
  
  # Annotating only for Observed plot
  annotate("text", x = 2.5, y = 385, label = "*", size = 5, hjust = 0.5) +
  annotate("text", x = 2, y = 410, label = "*", size = 5, hjust = 0.5) +
  annotate("text", x = 2.5, y = 440, label = "***", size = 5, hjust = 0.5) +
  annotate("text", x = 1.5, y = 355, label = "*", size = 5, hjust = 0.5,color = "red") +
  annotate("text", x = 3.5, y = 355, label = "***", size = 5, hjust = 0.5,color = "red") +

  # Adding branches (lines) for Shannon plot
  geom_segment(aes(x = 1, xend = 4, y = 380, yend = 380), color = "black", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 375, yend = 380), color = "black", size = 0.7) +
  geom_segment(aes(x = 4, xend = 4, y = 375, yend = 380), color = "black", size = 0.7) +
  
  geom_segment(aes(x = 1, xend = 3, y = 405, yend = 405), color = "black", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 400, yend = 405), color = "black", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 400, yend = 405), color = "black", size = 0.7) +
  
  geom_segment(aes(x = 2, xend = 3, y = 435, yend = 435), color = "black", size = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = 430, yend = 435), color = "black", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 430, yend = 435), color = "black", size = 0.7)+
  
  geom_segment(aes(x = 1, xend = 2, y = 350, yend = 350), color = "red", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 345, yend = 350), color = "red", size = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = 345, yend = 350), color = "red", size = 0.7)+
  
  geom_segment(aes(x = 3, xend = 4, y = 350, yend = 350), color = "red", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 345, yend = 350), color = "red", size = 0.7) +
  geom_segment(aes(x = 4, xend = 4, y = 345, yend = 350), color = "red", size = 0.7)
gg_richness_observed

# Save Observed plot
ggsave(filename = "R_Files/AlphaDiversity_Plots/gg_richness_observed_upf.png", 
       gg_richness_observed, height = 4, width = 6)


# Faith's PD Plot (pd_plot)
phylo_dist <- pd(t(otu_table(upf_phyloseq_rare)), phy_tree(upf_phyloseq_rare), include.root = FALSE)
sample_data(upf_phyloseq_rare)$PD <- phylo_dist$PD

pd_plot <- ggplot(sample_data(upf_phyloseq_rare), aes(x = upf_asthma, y = PD, fill = upf_status)) + 
  geom_point(position = position_dodge(width = 0.8), size = 1.5, alpha = 1) +  # Align points without changing color
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot alignment
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma","No Asthma", "Asthma")) +  
  ylab("Phylogenetic Diversity") +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11), 
        axis.title.y = element_text(size = 11)) + 
  facet_grid(. ~ "Faith's PD") +  # Adjust facet label
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 8, colour = "black"))+

  annotate("text", x = 2, y = 33, label = "*", size = 5, hjust = 0.5) +
  annotate("text", x = 2.5, y = 35.5, label = "***", size = 5, hjust = 0.5) +
  
  annotate("text", x = 1.5, y = 29, label = "*", size = 5, hjust = 0.5,color = "red") +
  annotate("text", x = 3.5, y = 29, label = "***", size = 5, hjust = 0.5,color = "red") +
  
  # Adding branches (lines) for Shannon plot
  geom_segment(aes(x = 1, xend = 3, y = 32.5, yend = 32.5), color = "black", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 32, yend = 32.5), color = "black", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 32, yend = 32.5), color = "black", size = 0.7) +
  
  geom_segment(aes(x = 2, xend = 3, y = 35, yend = 35), color = "black", size = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = 34.5, yend = 35), color = "black", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 34.5, yend = 35), color = "black", size = 0.7) +
  
  geom_segment(aes(x = 1, xend = 2, y = 28.5, yend = 28.5), color = "red", size = 0.7) +
  geom_segment(aes(x = 1, xend = 1, y = 28, yend = 28.5), color = "red", size = 0.7) +
  geom_segment(aes(x = 2, xend = 2, y = 28, yend = 28.5), color = "red", size = 0.7) +
  
  geom_segment(aes(x = 3, xend = 4, y = 28.5, yend = 28.5), color = "red", size = 0.7) +
  geom_segment(aes(x = 3, xend = 3, y = 28, yend = 28.5), color = "red", size = 0.7) +
  geom_segment(aes(x = 4, xend = 4, y = 28, yend = 28.5), color = "red", size = 0.7) 

pd_plot 
# Save pd_plot
ggsave(filename = "R_Files/AlphaDiversity_Plots/pd_plot_high_low_upf.png", 
       pd_plot, height = 4, width = 3)


# Combined Plot
combined_plot <- gg_richness_shannon +gg_richness_observed+ pd_plot + 
  plot_layout(ncol = 3, widths = c(1, 1))  
combined_plot

# Save combined plot
ggsave(filename = "R_Files/AlphaDiversity_Plots/combined_high_low_upf.png", 
       combined_plot, height = 4, width = 10)




## Stats
# Observed and Shannon
alphadiv <- estimate_richness(upf_phyloseq_rare)
samp_dat <- sample_data(upf_phyloseq_rare)
samp_dat_wdiv <- data.frame(upf_phyloseq_rare, alphadiv)

wilcox.test(Observed ~ asthma, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ asthma, data=samp_dat_wdiv, exact = FALSE)

# Faith's PD
phylo_dist <- pd(t(otu_table(upf_phyloseq_rare)), phy_tree(upf_phyloseq_rare),
                      include.root=F) 
sample_data(upf_phyloseq_rare)$PD <- phylo_dist$PD
pd_values <- sample_data(upf_phyloseq_rare)$PD
asthma <- sample_data(upf_phyloseq_rare)$asthma
wilcox.test(pd_values ~ asthma)

