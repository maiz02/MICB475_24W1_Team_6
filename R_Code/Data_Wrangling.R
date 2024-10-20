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
                      ifelse(x == "Edinburgh", "United Kingdom", 
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

# exporting the filtered_msmeta to a TSV file
output_filepath <- "filtered_ms_metadata.tsv"
write_tsv(filtered_msmeta, output_filepath)

#### RECONCILING MANIFEST BASED ON FILTERED DATA ####
# taking only the sample-id column
filtered_metadata_only_samples <- select(filtered_msmeta, `sample-id`)

# joining the filtered sample ids with the manifest file
filtered_manifest <- left_join(filtered_metadata_only_samples, msmanifest)

# exporting the filtered_manifest to a TSV file
filtered_manifest_filepath <- "filtered_ms_manifest.tsv"
write_tsv(filtered_manifest, filtered_manifest_filepath)
