library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

# creating rare file for PCoA plot
load("R files/mpt_final.RData")
ms_rare <- rarefy_even_depth(mpt_final, rngseed = 1, sample.size = 6000)
save(ms_rare, file="R files/ms_rare.RData")

metadata <- ms_rare@sam_data
metadata$asthma <- factor(metadata$asthma, levels = c(0, 1), labels = c("Without asthma", "With asthma"))


# all countries


# bray-curtis
bc_dm <- distance(ms_rare, method="bray")
pcoa_bc <- ordinate(ms_rare, method="PCoA", distance=bc_dm)
plot_ordination(ms_rare, pcoa_bc, color = "upf_status", 
                shape = "asthma") +
  scale_shape_manual(values = c(16, 17))


#### Beta diversity #####
bc_dm <- distance(mpt_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(mpt_rare, method="PCoA", distance=bc_dm)

plot_ordination(mpt_rare, pcoa_bc, color = "body.site", shape="subject")

gg_pcoa <- plot_ordination(mpt_rare, pcoa_bc, color = "body.site", shape="subject") +
  labs(pch="Subject #", col = "Body Site")
gg_pcoa

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)