library(tidyverse)

# loading in metadata
msmetafp <- "MICB475_24W1_Team_6/MS_Files/corrected_ms_metadata.tsv"
msmeta <- read_delim(msmetafp)

msmanifestfp <- "MICB475_24W1_Team_6/MS_Files/ms_manifest.tsv"
msmanifest <- read_delim(msmanifestfp)

# filtering for the columns we need
msmeta_col <- select(msmeta, "sample-id", "site_x", "allergies", "asthma")

# combining regions into their countries
x <- msmeta_col$site_x

site_country <- ifelse(x %in% c("San Francisco", "Boston", "New York", "Pittsburgh"), 
                      "USA", 
                      ifelse(x == "Edinburgh", "Scotland", 
                             ifelse(x == "Buenos Aires", "Argentina", 
                                    ifelse(x == "San Sebastian", "Spain", x))))
msmeta_col$site_x <- site_country

library(dplyr)
library(readr)

# dropping samples with NA values
msmeta_col <- msmeta_col %>%
  drop_na(asthma, allergies)

#filtering healthy individuals
msmeta_col <- msmeta_col %>%
  filter(disease != "MS")

#generating filtered ms metadata set with the remaining columns
filtered_msmeta <- msmeta_col %>%
  left_join(msmeta, by = "sample-id")

filtered_msmeta <- filtered_msmeta %>%
  rename("ms_status" = "disease.x", "country" = "site_x.x", "allergies" = "allergies.x", "asthma" = "asthma.x", "region" = "site_x.y")

# exporting the filtered_msmeta to a TSV file
output_filepath <- "filtered_ms_metadata.tsv"
write_tsv(filtered_msmeta, output_filepath)
