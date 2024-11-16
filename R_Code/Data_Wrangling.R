library(tidyverse)
library(dplyr)
library(readr)

# loading in metadata
msmetafp <- "MS_Files/corrected_ms_metadata.tsv"
msmeta <- read_delim(msmetafp)

msmanifestfp <- "MS_Files/ms_manifest.tsv"
msmanifest <- read_delim(msmanifestfp)

#### SUBSETTING DATA ####
# filtering for the columns we need
msmeta_col <- select(msmeta, "sample-id", "site_x", "allergies", "asthma", "disease")

# combining regions into their countries
x <- msmeta_col$site_x

site_country <- ifelse(x %in% c("San Francisco", "Boston", "New York", "Pittsburgh"), 
                      "USA", 
                      ifelse(x == "Edinburgh", "UK", 
                             ifelse(x == "Buenos Aires", "Argentina", 
                                    ifelse(x == "San Sebastian", "Spain", x))))
msmeta_col$site_x <- site_country

# dropping samples with NA values
msmeta_col <- msmeta_col %>%
  drop_na(asthma, allergies)

# filtering healthy individuals
msmeta_col <- msmeta_col %>%
  filter(disease != "MS")

# combining columns
msmeta_col <- msmeta_col %>%
  unite(country_allergies, site_x, allergies, sep=", ", remove=FALSE) %>%
  unite(country_asthma, site_x, asthma, sep=", ", remove=FALSE) %>%
  unite(allergies_and_asthma, allergies, asthma, sep=", ", remove=FALSE)

# generating filtered ms metadata set with the remaining columns
filtered_msmeta <- msmeta_col %>%
  left_join(msmeta, by = "sample-id")

filtered_msmeta <- filtered_msmeta %>%
  rename("ms_status" = "disease.x", "country" = "site_x.x", "allergies" = "allergies.x", "asthma" = "asthma.x", "region" = "site_x.y")

#adding the new column based on upf status -> high_upf(USA/UK) and low_upf(Spain/Argentina)
filtered_msmeta <- filtered_msmeta %>%
  mutate(upf_status = ifelse(country %in% c("USA", "UK"), "upf_high",
                             ifelse(country %in% c("Argentina", "Spain"), "upf_low", NA)))

#adding the new columns that indicates upf status and allergy/asthma conditions
filtered_msmeta <- filtered_msmeta %>%
  mutate(upf_allergies = paste(upf_status, allergies, sep = ", ")) %>%
  mutate(upf_asthma = paste (upf_status, asthma, sep = ", "))

#duplicating sample-id column for proper alpha rarefaction and reordering columns
filtered_msmeta$`sample-id_2` <- filtered_msmeta$`sample-id`
filtered_msmeta <- filtered_msmeta %>%
  select("sample-id", "sample-id_2", "upf_status","upf_allergies", "upf_asthma", everything())

#adding new allergy columns for pcoa plot
filtered_msmeta <- filtered_msmeta %>%
  mutate(allergies_yn = recode(allergies, `0`="no", `1`="yes")) %>%
  select("sample-id", "sample-id_2", "upf_status","upf_allergies", "upf_asthma", "allergies", "allergies_yn", everything())

#adding new asthma columns for pcoa plot
filtered_msmeta <- filtered_msmeta %>%
  mutate(asthma_yn = recode(asthma, `0`="no", `1`="yes")) %>%
  select("sample-id", "sample-id_2", "upf_status","upf_allergies", "upf_asthma", 
         "allergies", "allergies_yn", "asthma", "asthma_yn", everything())

# exporting the filtered_msmeta to a TSV file
output_filepath <- "MS_Files/final_filtered_ms_metadata.tsv"
write_tsv(filtered_msmeta, output_filepath)

#### RECONCILING MANIFEST BASED ON FILTERED DATA ####
# taking only the sample-id column
filtered_metadata_only_samples <- select(filtered_msmeta, `sample-id`)

# joining the filtered sample ids with the manifest file
filtered_manifest <- left_join(filtered_metadata_only_samples, msmanifest)

# exporting the filtered_manifest to a TSV file
filtered_manifest_filepath <- "filtered_ms_manifest.tsv"
write_tsv(filtered_manifest, filtered_manifest_filepath)


#####Creating phyloseq object #####
#load in metadata
meta <- read_delim(file = "filtered_ms_metadata.tsv", delim = "\t")
#load in features table
otu <- read_delim(file="feature-table.txt", delim = "\t", skip=1)

#load in your taxonomy file
tax <- read_delim(file = "taxonomy.tsv", delim="\t")

#load in tree
phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table into phyloseq object ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata into phyloseq object ####
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)

#### Formatting taxonomy table ####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)


#Create phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#### PRUNE out the bad ASVs ####
# Remove ASVs that have less than 5 counts total with prune (typically mistakes/seq error)
mpt_filt_nolow <- filter_taxa(mpt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads (ex. poor sequencing run) 
mpt_filt_nolow_samps <- prune_samples(sample_sums(mpt_filt_nolow)>100, mpt_filt_nolow)
# !is.na(month) --> keeping anything that ISN'T N/A!
mpt_final <- subset_samples(mpt_filt_nolow_samps, !is.na(month) )


#save 
save(mpt, file="mpt.RData")
save(mpt_final, file="mpt_final.RData")



