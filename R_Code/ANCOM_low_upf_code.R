
library(tidyverse) # For all your basic data wrangling and plotting needs.
library(phyloseq) # Indispensable package for microbiome analyses. 
library(ggpubr)

set.seed(711)

#Load phyloseq
load("upf_phyloseq_final_low.RData")

#clean it up
hist(sample_sums(upf_phyloseq_final_low),breaks = 30) 
hist(log10(sample_sums(upf_phyloseq_final_low)),breaks = 30) 

table(below_1000 = sample_sums(upf_phyloseq_final_low)<=1000,
      expt = upf_phyloseq_final_low@sam_data$asthma_yn) 

upf_phyloseq_final_low = prune_samples(sample_sums(upf_phyloseq_final_low) >= 1000, upf_phyloseq_final_low)

table(below_1000 = sample_sums(upf_phyloseq_final_low)<=1000,
      expt = upf_phyloseq_final_low@sam_data$asthma_yn) 

#aggregate the data to the Family level for our analysis.
family = tax_glom(upf_phyloseq_final_low,'Family')
ntaxa(upf_phyloseq_final_low); ntaxa(family) 
#[1] 2056
#[1] 82


calculate_relative_abundance <- function(x) x / sum(x)
total_counts <- taxa_sums(family)
relative_abundance <- calculate_relative_abundance(total_counts) 
abundant <- relative_abundance > 0.001 
family <- prune_taxa(abundant, family) 
family 

#ANCOM
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

# Input 
ancom.family = ancombc(phyloseq = family, # Raw counts
                       formula = 'asthma_yn', # The explanatory variable
                       p_adj_method = "fdr",
                       prv_cut=0.10, # Max proportion of zeros allowed per taxon
                       lib_cut = 1000, # Can filter out samples below minimum seq depth here
                       group = 'asthma_yn', # If you're including structural zeros below, you need thisaa
                       struc_zero = T) # If true, any taxa present in only 1 category of 'group' will automatically be significant
str(ancom.family)
View(ancom.family)
#results = format_ancom_results(ancom.family,family,level_letter = 'f')

# Update the column names for each part of the results, because they're currently all the same.
colnames(ancom.family$res$lfc) = paste(colnames(ancom.family$res$lfc),'_beta',sep='')
colnames(ancom.family$res$se) = paste(colnames(ancom.family$res$se),'_se',sep='')
colnames(ancom.family$res$W) = paste(colnames(ancom.family$res$W),'_W',sep='')
colnames(ancom.family$res$p_val) = paste(colnames(ancom.family$res$p_val),'_p_val',sep='')
colnames(ancom.family$res$q_val) = paste(colnames(ancom.family$res$q_val),'_q_val',sep='')
colnames(ancom.family$res$diff_abn) = paste(colnames(ancom.family$res$diff_abn),'_diff_abn',sep='')

results = lapply(ancom.family$res,function(x) rownames_to_column(x,'Family')) %>% 
  lapply(as_tibble) %>% reduce(full_join)
results = results %>% dplyr::select(-contains('Intercept'))
srv = Vectorize(str_remove) 
colnames(results) = srv(colnames(results),'asthma_ynyes_') 
head(results)

library(microbiome)
family_tss = family %>% transform('compositional')
selected_taxa <- taxa_names(family)
family_tss <- prune_taxa(selected_taxa, family_tss)
family_tss_melt = family_tss %>% psmelt()

family_tss_melt = family_tss_melt %>% 
  dplyr::select(-c(OTU,Domain:Order)) %>%
  pivot_wider(names_from = Family, values_from = Abundance)

# Look at which p values have less than 0.05, and match the sample ID to the one seen in tax_table(family)
library(writexl)
write_xlsx(list('all_results' = results),
           'ancom_low_upf_results_family.xlsx')                
tax_table(family)


#put bug as whichever family is significant to make each graph. 

##Butyricicoccaceae
bug = "f__Butyricicoccaceae"
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   

#Coriobacteriaceae
bug = "f__Coriobacteriaceae"
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   


#Erysipelatoclostridiaceae
bug = "f__Erysipelatoclostridiaceae"  
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   


#Lachnospiraceae 
bug = "f__Lachnospiraceae"
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   



#Ruminococcaceae
bug = "f__Ruminococcaceae"
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   

#UCG-010
bug = "f__UCG-010"
  for (b in bug) {
  # Check if the column exists
  if (!b %in% colnames(family_tss_melt)) {
    warning(paste("Column", b, "not found in family_tss_melt. Skipping."))
    next
  }
  
  # Check if the column is numeric
  if (!is.numeric(family_tss_melt[[b]])) {
    family_tss_melt[[b]] <- as.numeric(family_tss_melt[[b]])
    if (anyNA(family_tss_melt[[b]])) {
      warning(paste("Column", b, "contains non-convertible values. Skipping."))
      next
    }
  }

  # Generate random colors
  colors <- c("blue", "red")

  # Create the plot
  p <- ggplot(family_tss_melt, aes(x = asthma_yn, y = family_tss_melt[[b]], fill = asthma_yn)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    xlab("Presence of Asthma") +
    ylab("% Ab")
  # Save the plot
  print(p)
  ggsave(paste0("tss_", b, ".jpeg"), plot = p, height = 5, width = 5)
}   


# Saving your results
# Save ANCOM results as a whole:
saveRDS(ancom.family,'R_Files/ANCOM/low_upf_ancom_results_family.rds')
# ancom.family = readRDS('R_Files/ANCOM/low_upf_ancom_results_family.rds') # to load .rds files

# Save your formatted results table as a .csv
write.csv(results.formatted,'ancom_results_family.csv',row.names = F) 
# We don't have any row names assigned, so it would have just created an irrelevant column of increasing numbers

# Save your formatted results table as an excel spreadsheet, where the significant results are on a separate sheet
library(writexl)
write_xlsx(list('all_results' = results),
           'R_Files/ANCOM/low_upf_ancom_results_family.xlsx')


