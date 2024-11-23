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
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df_high <- as.data.frame(meta_high[,-1:-2])
# Make sampleids the rownames
rownames(samp_df_high)<- meta_high$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP_high <- sample_data(samp_df_high)
class(SAMP_high)

samp_df_low <- as.data.frame(meta_low[,-1:-2])
# Make sampleids the rownames
rownames(samp_df_low)<- meta_low$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP_low <- sample_data(samp_df_low)
class(SAMP_low)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
upf_phyloseq_high<- phyloseq(OTU, SAMP_high, TAX, phylotree)
upf_phyloseq_low<- phyloseq(OTU, SAMP_low, TAX, phylotree)


#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(upf_phyloseq_high)
sample_data(upf_phyloseq_high)
tax_table(upf_phyloseq_high)
phy_tree(upf_phyloseq_high)

otu_table(upf_phyloseq_low)
sample_data(upf_phyloseq_low)
tax_table(upf_phyloseq_low)
phy_tree(upf_phyloseq_low)

#### Accessor functions #### UPF HIGH
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(upf_phyloseq_high)
# colnames(sample_data(atamaca))
get_variable(upf_phyloseq_high, c("upf_status","asthma")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names:
sample_names(upf_phyloseq_high)
# How many samples do we have?
nsamples(upf_phyloseq_high) #425
# What is the sum of reads in each sample?
sample_sums(upf_phyloseq_high)
# Save the sample names of the 3 samples with the most reads
getsamps_high <- names(sort(sample_sums(upf_phyloseq_high), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(upf_phyloseq_high, getsamps_high) 

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(upf_phyloseq_high)
# How many taxa do we have?
ntaxa(upf_phyloseq_high)
# What is the total read count for each taxa?
taxa_sums(upf_phyloseq_high)
# Let's find the top 3 most abundant asvs (taxa) in our data
gettaxa_high <- names(sort(taxa_sums(upf_phyloseq_high), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq_high, gettaxa_high)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
upf_phyloseq_filt_high <- subset_taxa(upf_phyloseq_high,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
upf_phyloseq_filt_nolow_high <- filter_taxa(upf_phyloseq_filt_high, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
upf_phyloseq_nolow_samps_high <- prune_samples(sample_sums(upf_phyloseq_filt_nolow_high)>100, upf_phyloseq_filt_nolow_high)
# Remove samples where upf_status is na
upf_phyloseq_final_high <- subset_samples(upf_phyloseq_nolow_samps_high, !is.na(upf_status) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final_high))), cex=0.1)
upf_phyloseq_rare_high <- rarefy_even_depth(upf_phyloseq_final_high, rngseed = 1, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final_high, file="R_Code/upf_phyloseq_final_high.RData")
save(upf_phyloseq_rare_high, file="R_Code/upf_phyloseq_rare_high.RData")




#### Accessor functions #### UPF LOW
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(upf_phyloseq_low)
# colnames(sample_data(atamaca))
get_variable(upf_phyloseq_low, c("upf_status","asthma")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names:
sample_names(upf_phyloseq_low)
# How many samples do we have?
nsamples(upf_phyloseq_low) #425
# What is the sum of reads in each sample?
sample_sums(upf_phyloseq_low)
# Save the sample names of the 3 samples with the most reads
getsamps_low <- names(sort(sample_sums(upf_phyloseq_low), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(upf_phyloseq_low, getsamps_low) 

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(upf_phyloseq_low)
# How many taxa do we have?
ntaxa(upf_phyloseq_low)
# What is the total read count for each taxa?
taxa_sums(upf_phyloseq_low)
# Let's find the top 3 most abundant asvs (taxa) in our data
gettaxa_low <- names(sort(taxa_sums(upf_phyloseq_low), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq_low, gettaxa_low)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
upf_phyloseq_filt_low <- subset_taxa(upf_phyloseq_low,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
upf_phyloseq_filt_nolow_low <- filter_taxa(upf_phyloseq_filt_low, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
upf_phyloseq_nolow_samps_low <- prune_samples(sample_sums(upf_phyloseq_filt_nolow_low)>100, upf_phyloseq_filt_nolow_low)
# Remove samples where upf_status is na
upf_phyloseq_final_low <- subset_samples(upf_phyloseq_nolow_samps_low, !is.na(upf_status) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final_low))), cex=0.1)
upf_phyloseq_rare_low <- rarefy_even_depth(upf_phyloseq_final_low, rngseed = 2, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final_low, file="R_Code/upf_phyloseq_final_low.RData")
save(upf_phyloseq_rare_low, file="R_Code/upf_phyloseq_rare_low.RData")

