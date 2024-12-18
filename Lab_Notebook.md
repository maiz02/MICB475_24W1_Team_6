**October 5** <br>
Yna <br>
- filtered metadata file by selecting only for the columns to be used in our project (sample id, site, allergies, asthma)
- renamed each region into their respective country and updated the filtered metadata file
- uploaded the metadata and manifest files to the GitHub repository

**October 6** <br>
Julia <br>
- filtered out the rows of MS patients to obtain samples of healthy individuals from the metadata filtered by Yna
- filtered out the rows containing NA values
- Merged the remaining columns of metadata with the columns of interest and exported the filtered dataset in TSV format, which was then uploaded to the GitHub repository

**October 7** <br>
Andrea <br>
- reconciled the filtered metadata and manifest files to create a new, filtered manifest file

**October 13** <br>
Yna <br>
- combined columns that we will be using for our data processing in QIIME2
- updated the metadata and manifest files to reflect this
- uploaded the corrected files into the GitHub repository

**October 19** <br>
Yna <br>
- updated the code and the metadata file to change Scotland to UK

Andrea <br>
- completed all data processing in QIIME2
- ran initial beta diversity metrics for asthma and allergy groups

**October 23** <br>
Julia <br>
- updated the code and the metadata file to combine high/low UPF countries

**November 15** <br>
Yna <br>
- updated metadata to add new columns (allergies_yn, asthma_yn) to be used for the PCoA plots
- updated phyloseq object code in data wrangling code to reflect the updated metadata
- generated PCoA plots

**November 17** <br>
Yna <br>
- made new PCoA plots that show UPF status, allergies, and asthma all in one plot

**November 18** <br>
Chaeyoon <br>
- ran PERMANOVA with the beta diversity metrics (weighted unifrac, bray-curtis, jaccard) and tabulated results

**November 17** <br>
Yna <br>
- updated PCoA plot for weighted unifrac
- generated alpha diversity metrics

**November 22** <br>
Julia <br>
- generated phyloseq objects of upf_high / upf_low groups

**November 23** <br>
Andrea <br>
- ran ISA on high and low UPF groups 

**November 24** <br>
Chaeyoon <br>
- ran core microbiome analysis and generated venn diagrams for high and low upf groups
- ran ancom on low upf group

**November 26** <br>
Julia <br>
- merged alpha diversity plots by placing them on the same panel for the manuscript figure

**December 9** <br>
Yna <br>
- generated stats for new alpha diversity analysis

Andrea <br>
- generated updated PCoA plots to show clusters according to UPF status
