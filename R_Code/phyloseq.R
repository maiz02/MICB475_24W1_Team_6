library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

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
samp_df <- as.data.frame(meta[,-1:-2])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

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
upf_phyloseq<- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(upf_phyloseq)
sample_data(upf_phyloseq)
tax_table(upf_phyloseq)
phy_tree(upf_phyloseq)


#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(upf_phyloseq)
# colnames(sample_data(atamaca))
get_variable(upf_phyloseq, c("upf_status","asthma")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names:
sample_names(upf_phyloseq)
# How many samples do we have?
nsamples(upf_phyloseq) #425
# What is the sum of reads in each sample?
sample_sums(upf_phyloseq)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(upf_phyloseq), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(upf_phyloseq, getsamps) 

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(upf_phyloseq)
# How many taxa do we have?
ntaxa(upf_phyloseq) #3097
# What is the total read count for each taxa?
taxa_sums(upf_phyloseq)
# Let's find the top 3 most abundant asvs (taxa) in our data
gettaxa <- names(sort(taxa_sums(upf_phyloseq), decreasing = TRUE)[1:3] )
get_sample(upf_phyloseq, gettaxa)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
upf_phyloseq_filt <- subset_taxa(upf_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
upf_phyloseq_filt_nolow <- filter_taxa(upf_phyloseq_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
upf_phyloseq_nolow_samps <- prune_samples(sample_sums(upf_phyloseq_filt_nolow)>100, upf_phyloseq_filt_nolow)
# Remove samples where upf_status is na
upf_phyloseq_final <- subset_samples(upf_phyloseq_nolow_samps, !is.na(upf_status) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(upf_phyloseq_final))), cex=0.1)
upf_phyloseq_rare <- rarefy_even_depth(upf_phyloseq_final, rngseed = 1, sample.size = 1000)

##### Saving #####
save(upf_phyloseq_final, file="R_Code/upf_phyloseq_final.RData")
save(upf_phyloseq_rare, file="R_Code/upf_phyloseq_rare.RData")

