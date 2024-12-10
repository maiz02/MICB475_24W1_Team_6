install.packages("patchwork")

library(phyloseq)
library(tidyverse)
library(picante)
library(patchwork)
library(dplyr)
library(ape)

load("R_Code/upf_phyloseq_final.RData")
load("R_Code/upf_phyloseq_rare.RData")

upf_phyloseq_rare <- upf_phyloseq_rare_high

###Two way ANOVA

mdl <- lm(Observed ~ upf_status*asthma, data=samp_dat_wdiv)
summary(aov(mdl))


##### figure generating part
samp_dat_wdiv <- data.frame(sample_data(upf_phyloseq_rare), estimate_richness(upf_phyloseq_rare))
format_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", signif(p, 3)))
  }
}



obs_asthma <- format_p(wilcox.test(Observed ~ upf_status, data = samp_dat_wdiv, exact = FALSE)$p.value)
shan_asthma <- format_p(wilcox.test(Shannon ~ upf_status, data = samp_dat_wdiv, exact = FALSE)$p.value)


# Alpha Diversity Plot (gg_richness)
gg_richness <- plot_richness(upf_phyloseq_final, x = "upf_asthma", 
                             measures = c("Observed", "Shannon")) +
  xlab(" ") +
  scale_x_discrete(labels = c("No Asthma", "Asthma", "No Asthma", "Asthma")) + 
  geom_boxplot(aes(fill = upf_status), alpha = 0.7, outlier.shape = NA) +  
  scale_fill_manual(values = c("upf_high" = "steelblue", "upf_low" = "salmon"),
                    name = "UPF Status") + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none") 

# Print the final plot
# Define annotation data for p-values and brackets
brackets_asthma <- data.frame(
  x = c(1, 4),                # Start and end positions for the bracket
  xend = c(4, 1),             # End positions for the bracket
  y = c(350, 8),            # Vertical position for the brackets
  variable = c("Observed", "Shannon")  # Facets to target
)

annotations_asthma <- data.frame(
  x = c(2.5, 2.5),            # Midpoint for p-value text
  y = c(360, 8.3),            # Position for p-value text
  label = c(obs_asthma, shan_asthma),
  variable = c("Observed", "Shannon")  # Facets to target
)

# Add brackets and annotations to the plot
gg_richness <- gg_richness +
  # Brackets for asthma comparison
  geom_segment(
    data = brackets_asthma,
    aes(x = x, xend = xend, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  # P-value annotations for asthma comparison
  geom_text(
    data = annotations_asthma,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE
  )

# Print the updated plot
gg_richness



ggsave(filename = "R_Files/AlphaDiversity_Plots/gg_richness_high_low_upf.png", 
       gg_richness, height = 4, width = 6)


# Faith's PD Plot (pd_plot)
phylo_dist_high <- pd(t(otu_table(upf_phyloseq_rare_high)), phy_tree(upf_phyloseq_rare_high), include.root = FALSE)
phylo_dist_low <- pd(t(otu_table(upf_phyloseq_rare_low)), phy_tree(upf_phyloseq_rare_low), include.root = FALSE)

# Add PD to sample_data in the phyloseq objects
sample_data(upf_phyloseq_rare_high)$PD <- phylo_dist_high$PD
sample_data(upf_phyloseq_rare_low)$PD <- phylo_dist_low$PD

# Convert updated sample_data to data frames
samp_dat_high <- as.data.frame(sample_data(upf_phyloseq_rare_high))
samp_dat_low <- as.data.frame(sample_data(upf_phyloseq_rare_low))

samp_dat_high$asthma <- as.factor(samp_dat_high$asthma)
samp_dat_low$asthma <- as.factor(samp_dat_low$asthma)

# Perform Wilcoxon test for PD
pd_asthma_high <- format_p(wilcox.test(PD ~ asthma, data = samp_dat_high, exact = FALSE)$p.value)
pd_asthma_low <- format_p(wilcox.test(PD ~ asthma, data = samp_dat_low, exact = FALSE)$p.value)


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

# Define annotation data for PD
pd_annotations <- data.frame(
  x = c(1.5, 3.5),
  y = c(5, 7),  # Adjust the y-values based on your data range for PD
  label = c(pd_asthma_high, pd_asthma_low),
  variable = "Phylogenetic Diversity"
)

# Define brackets for PD facet
brackets_pd <- data.frame(
  x = c(1, 3),
  xend = c(2, 4),
  y = c(6, 8),  # Adjust the y-values based on your data range for PD
  variable = "Phylogenetic Diversity"
)

# Modify pd_plot to include brackets and annotations
pd_plot <- pd_plot +
  geom_segment(data = brackets_pd, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
  geom_text(data = pd_annotations, aes(x = x, y = y, label = label), inherit.aes = FALSE)

# Print the final phylogenetic diversity plot with annotations
pd_plot


pd_plot 

ggsave(filename = "R_Files/AlphaDiversity_Plots/pd_plot_high_low_upf.png", 
       pd_plot, height = 4, width = 3)



combined_plot <- gg_richness + pd_plot + 
  plot_layout(ncol = 2, widths = c(2, 1))  
combined_plot

ggsave(filename = "R_Files/AlphaDiversity_Plots/combined_high_low_upf.png", 
       combined_plot, height = 4, width = 10)




## Stats - Wilcoxon test ##
# Observed and Shannon
samp_dat <- sample_data(upf_phyloseq_rare)
alphadiv <- estimate_richness(upf_phyloseq_rare, measures = c("Observed", "Shannon"))
samp_dat_wdiv <- cbind(as.data.frame(samp_dat), alphadiv)
wilcox.test(Observed ~ asthma, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Observed ~ upf_status, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ asthma, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ upf_status, data=samp_dat_wdiv, exact = FALSE)

# High with asthma + low no asthma
metadata <- sample_data(upf_phyloseq_rare)
metadata_df <- as.data.frame(metadata)
observed_features <- estimate_richness(upf_phyloseq_rare, measures = "Observed")
shannon_diversity <- estimate_richness(upf_phyloseq_rare, measures = "Shannon")
metadata_df$observed_features <- observed_features$Observed
metadata_df$shannon_diversity <- shannon_diversity$Shannon

group_high_observed_1 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_high, 1"]
group_low_observed_1 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_low, 0"]
wilcox_test_observed_1 <- wilcox.test(group_high_observed_1, group_low_observed_1, alternative = "two.sided")
wilcox_test_observed_1

group_high_shannon_1 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_high, 1"]
group_low_shannon_1 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_low, 0"]
wilcox_test_shannon_1 <- wilcox.test(group_high_shannon_1, group_low_shannon_1, alternative = "two.sided")
wilcox_test_shannon_1

# High with asthma + low with asthma
group_high_observed_2 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_high, 1"]
group_low_observed_2 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_low, 1"]
wilcox_test_observed_2 <- wilcox.test(group_high_observed_2, group_low_observed_2, alternative = "two.sided")
wilcox_test_observed_2

group_high_shannon_2 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_high, 1"]
group_low_shannon_2 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_low, 1"]
wilcox_test_shannon_2 <- wilcox.test(group_high_shannon_2, group_low_shannon_2, alternative = "two.sided")
wilcox_test_shannon_2

# High no asthma + low no asthma
group_high_observed_3 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_high, 0"]
group_low_observed_3 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_low, 0"]
wilcox_test_observed_3 <- wilcox.test(group_high_observed_3, group_low_observed_3, alternative = "two.sided")
wilcox_test_observed_3

group_high_shannon_3 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_high, 0"]
group_low_shannon_3 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_low, 0"]
wilcox_test_shannon_3 <- wilcox.test(group_high_shannon_3, group_low_shannon_3, alternative = "two.sided")
wilcox_test_shannon_3

# High no asthma + low with asthma
group_high_observed_4 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_high, 0"]
group_low_observed_4 <- metadata_df$observed_features[metadata_df$upf_asthma == "upf_low, 1"]
wilcox_test_observed_4 <- wilcox.test(group_high_observed_4, group_low_observed_4, alternative = "two.sided")
wilcox_test_observed_4

group_high_shannon_4 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_high, 0"]
group_low_shannon_4 <- metadata_df$shannon_diversity[metadata_df$upf_asthma == "upf_low, 1"]
wilcox_test_shannon_4 <- wilcox.test(group_high_shannon_4, group_low_shannon_4, alternative = "two.sided")
wilcox_test_shannon_4


# Faith's PD
phylo_dist <- pd(t(otu_table(upf_phyloseq_rare)), phy_tree(upf_phyloseq_rare),
                 include.root=F) 
sample_data(upf_phyloseq_rare)$PD <- phylo_dist$PD
pd_values <- sample_data(upf_phyloseq_rare)$PD
asthma_condition <- sample_data(upf_phyloseq_rare)$upf_asthma
wilcox.test(pd_values ~ asthma)
wilcox.test(pd_values ~ upf_status)

# High with asthma + low no asthma
group_high_pd_1 <- pd_values[asthma_condition == "upf_high, 1"]
group_low_pd_1 <- pd_values[asthma_condition == "upf_low, 0"]
wilcox_test_pd_1 <- wilcox.test(group_high_pd_1, group_low_pd_1, alternative = "two.sided")
wilcox_test_pd_1

# High with asthma + low with asthma
group_high_pd_2 <- pd_values[asthma_condition == "upf_high, 1"]
group_low_pd_2 <- pd_values[asthma_condition == "upf_low, 1"]
wilcox_test_pd_2 <- wilcox.test(group_high_pd_2, group_low_pd_2, alternative = "two.sided")
wilcox_test_pd_2

# High no asthma + low no asthma
group_high_pd_3 <- pd_values[asthma_condition == "upf_high, 0"]
group_low_pd_3 <- pd_values[asthma_condition == "upf_low, 0"]
wilcox_test_pd_3 <- wilcox.test(group_high_pd_3, group_low_pd_3, alternative = "two.sided")
wilcox_test_pd_3

# High no asthma + low with asthma
group_high_pd_4 <- pd_values[asthma_condition == "upf_high, 0"]
group_low_pd_4 <- pd_values[asthma_condition == "upf_low, 1"]
wilcox_test_pd_4 <- wilcox.test(group_high_pd_4, group_low_pd_4, alternative = "two.sided")
wilcox_test_pd_4
