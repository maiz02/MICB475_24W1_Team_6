library(tidyverse)

# loading in metadata
msmetafp <- "corrected_ms_metadata.tsv"
msmeta <- read_delim(msmetafp)

msmanifestfp <- "ms_manifest.tsv"
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
