# October 7 - 3rd meeting

## Attendees
- Sam
- Andrea Garcia
- Yna Ortiz
- Julia Jung
- Vivian Tan
- Chaeyoon Chang

## Agenda
1) Making a new manifest file? What do we do with it?
2) New column combining country and allergies (yes or no) / country and asthma (yes or no) 
3) Difference between research objective vs experimental aims and rationale
4) Clarify steps: old metadata file --> qiime2 (denoising/demultiplexing) --> make new columns? (new metadata file) --> R (phyloseq)
5) Rarefaction: first through qiime to see the number, then through R
6) Our research question: exploring the regional differences in incidents of asthma and allergies
7) the # of samples for each countries?
   
## Minutes
1) We have done the wrangling (regions to countries, removed NA rows and MS patient rows), what do we do about the new manifest file?  Combining manifest and metadata either at the beginning or at the qiime level
2) Metadata and manifest: metadata has variables (control, response), directs code (what variable are associated with what). manifest has sample ids and reads. anything in your manifest not in the metadata will be ignored. 
3) Denoising and demultiplexing on qiime on the old file.
4) Rarefaction: in the manuscript, pick one graph (Evelyn likes R)
5) Columns: make combined columns (countries & allergies and countries & asthma --> 2 columns). Do numbers (1 for US, 2 for Spain) and for allergies/asthma variable
6) Our aims currently are technique-based not biology-based. Our aim: Data processing through qiime2, what we should have as aim: how does xxx affect yyy microbiome? Our research question: looking at regional differences first (differences in microbiome based on region (regions as variable), and then two separate aims (asthma and allergies), we should also include microbiome and diversity in our aims
   -1st aim: region as variable (maybe remove this as not as relevant to our research question)  
   -2nd aim: country and allergies: yes or no, split to 8 groups, all at once
   -3rd aim: country and asthma: yes or no (same as aim 2)
   -4th aim: relation between asthma and allergies (occurring together, as literatures states)
      -another column combining asthma and allergies (1,1)
   -Comparing the microbiome of people with and without allergies. Disease markers? or exploratory question or drug targets? Diversity of microbes in each country and asthma and allergies. Characterizing microbiomes, different diets/culture/cuisine for what region is associated with less allergies or asthma?  
9) We filtered the metadata (only the healthy,so ~400), if statement or for loop to keep whatever is in the metadata. Sample IDs match. WE have to subset the manifest file first.
10) Even though we have different sample sizes, just go with it (statistics takes it into account)
11) If we still need to wrangle, we can change proposal requirements
12) Our final metadata is the filtered one

 ## Action Items
1) Gantt chart
2) Background research and writing the intro
3) Subsetting the manifest file!!
4) Convert strings to numbers
