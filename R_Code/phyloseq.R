library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load data ####
metams <- "MS_Files/final_filtered_ms_metadata.tsv"
meta <- read_delim(metams, delim="\t")
meta_high<- meta %>%
  filter(upf_status == "upf_high")
meta_low<- meta %>%
  filter(upf_status == "upf_low")

otums <- "QIIME2_Files/export/feature-table.txt"
otu <- read_delim(file = otums, delim="\t", skip=1)

taxms <- "QIIME2_Files/export/taxonomy.tsv"
tax <- read_delim(taxms, delim="\t")

phylotreems <- "QIIME2_Files/export/tree.nwk"
phylotree <- read.tree(phylotreems)

#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
samp_df_high <- as.data.frame(meta_high[,-1:-2])
rownames(samp_df_high)<- meta_high$'sample-id'
SAMP_high <- sample_data(samp_df_high)
class(SAMP_high)

samp_df_low <- as.data.frame(meta_low[,-1:-2])
rownames(samp_df_low)<- meta_low$'sample-id'
SAMP_low <- sample_data(samp_df_low)
class(SAMP_low)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
upf_phyloseq_high<- phyloseq(OTU, SAMP_high, TAX, phylotree)
upf_phyloseq_low<- phyloseq(OTU, SAMP_low, TAX, phylotree)


#### Looking at phyloseq object #####
otu_table(upf_phyloseq_high)
sample_data(upf_phyloseq_high)
tax_table(upf_phyloseq_high)
phy_tree(upf_phyloseq_high)

otu_table(upf_phyloseq_low)
sample_data(upf_phyloseq_low)
tax_table(upf_phyloseq_low)
phy_tree(upf_phyloseq_low)

#### Accessor functions #### - UPF HIGH
sample_variables(upf_phyloseq_high)
get_variable(upf_phyloseq_high, c("upf_status","asthma"))


sample_names(upf_phyloseq_high)
nsamples(upf_phyloseq_high) #425
sample_sums(upf_phyloseq_high)
getsamps_high <- names(sort(sample_sums(upf_phyloseq_high), decreasing = TRUE)[1:3])
get_taxa(upf_phyloseq_high, getsamps_high) 


taxa_names(upf_phyloseq_high)
ntaxa(upf_phyloseq_high)
taxa_sums(upf_phyloseq_high)
gettaxa_high <- names(sort(taxa_sums(upf_phyloseq_high), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq_high, gettaxa_high)


######### ANALYZE ##########
upf_phyloseq_filt_high <- subset_taxa(upf_phyloseq_high,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
upf_phyloseq_filt_nolow_high <- filter_taxa(upf_phyloseq_filt_high, function(x) sum(x)>5, prune = TRUE)
upf_phyloseq_nolow_samps_high <- prune_samples(sample_sums(upf_phyloseq_filt_nolow_high)>100, upf_phyloseq_filt_nolow_high)
upf_phyloseq_final_high <- subset_samples(upf_phyloseq_nolow_samps_high, !is.na(upf_status) )

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final_high))), cex=0.1)
upf_phyloseq_rare_high <- rarefy_even_depth(upf_phyloseq_final_high, rngseed = 1, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final_high, file="R_Code/upf_phyloseq_final_high.RData")
save(upf_phyloseq_rare_high, file="R_Code/upf_phyloseq_rare_high.RData")




#### Accessor functions #### UPF LOW
sample_variables(upf_phyloseq_low)
get_variable(upf_phyloseq_low, c("upf_status","asthma"))

sample_names(upf_phyloseq_low)
nsamples(upf_phyloseq_low) 
sample_sums(upf_phyloseq_low)
getsamps_low <- names(sort(sample_sums(upf_phyloseq_low), decreasing = TRUE)[1:3])
get_taxa(upf_phyloseq_low, getsamps_low) 

taxa_names(upf_phyloseq_low)
ntaxa(upf_phyloseq_low)
taxa_sums(upf_phyloseq_low)
gettaxa_low <- names(sort(taxa_sums(upf_phyloseq_low), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq_low, gettaxa_low)


######### ANALYZE ##########
upf_phyloseq_filt_low <- subset_taxa(upf_phyloseq_low,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
upf_phyloseq_filt_nolow_low <- filter_taxa(upf_phyloseq_filt_low, function(x) sum(x)>5, prune = TRUE)
upf_phyloseq_nolow_samps_low <- prune_samples(sample_sums(upf_phyloseq_filt_nolow_low)>100, upf_phyloseq_filt_nolow_low)
upf_phyloseq_final_low <- subset_samples(upf_phyloseq_nolow_samps_low, !is.na(upf_status) )

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final_low))), cex=0.1)
upf_phyloseq_rare_low <- rarefy_even_depth(upf_phyloseq_final_low, rngseed = 2, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final_low, file="R_Code/upf_phyloseq_final_low.RData")
save(upf_phyloseq_rare_low, file="R_Code/upf_phyloseq_rare_low.RData")


#### Combined High and Low UPF Phyloseq Object ####
#### Load data ####
metams <- "MS_Files/final_filtered_ms_metadata.tsv"
meta <- read_delim(metams, delim="\t")

otums <- "QIIME2_Files/export/feature-table.txt"
otu <- read_delim(file = otums, delim="\t", skip=1)

taxms <- "QIIME2_Files/export/taxonomy.tsv"
tax <- read_delim(taxms, delim="\t")

phylotreems <- "QIIME2_Files/export/tree.nwk"
phylotree <- read.tree(phylotreems)

#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
samp_df <- as.data.frame(meta[,-1:-2])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
upf_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)


#### Looking at phyloseq object #####
otu_table(upf_phyloseq_high)
sample_data(upf_phyloseq_high)
tax_table(upf_phyloseq_high)
phy_tree(upf_phyloseq_high)

otu_table(upf_phyloseq_low)
sample_data(upf_phyloseq_low)
tax_table(upf_phyloseq_low)
phy_tree(upf_phyloseq_low)

#### Accessor functions #### - UPF HIGH
sample_variables(upf_phyloseq)
get_variable(upf_phyloseq, c("upf_status","asthma"))


sample_names(upf_phyloseq)
nsamples(upf_phyloseq) #425
sample_sums(upf_phyloseq)
getsamps <- names(sort(sample_sums(upf_phyloseq), decreasing = TRUE)[1:3])
get_taxa(upf_phyloseq, getsamps) 


taxa_names(upf_phyloseq)
ntaxa(upf_phyloseq)
taxa_sums(upf_phyloseq)
gettaxa <- names(sort(taxa_sums(upf_phyloseq), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq, gettaxa)


######### ANALYZE ##########
upf_phyloseq_filt <- subset_taxa(upf_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
upf_phyloseq_filt_nolow <- filter_taxa(upf_phyloseq_filt, function(x) sum(x)>5, prune = TRUE)
upf_phyloseq_nolow_samps <- prune_samples(sample_sums(upf_phyloseq_filt_nolow)>100, upf_phyloseq_filt_nolow)
upf_phyloseq_final <- subset_samples(upf_phyloseq_nolow_samps, !is.na(upf_status) )

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final))), cex=0.1)
upf_phyloseq_rare_high <- rarefy_even_depth(upf_phyloseq_final, rngseed = 1, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final, file="R_Code/upf_phyloseq_final.RData")
save(upf_phyloseq_rare_high, file="R_Code/upf_phyloseq_rare.RData")
