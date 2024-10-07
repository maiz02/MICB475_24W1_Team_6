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
2) Metadata and manifest: metadata has variables (control, response), directs code (what variable are associated with what). manifest has sample ids and reads. anything in your manifest not in the metadata will be ignored. Our final metadata is the filtered one
3) Rarefaction: in the manuscript, pick one graph (Evelyn likes R). Rarefaction on R with the filtered file. 
4) Columns: make combined columns (countries & allergies and countries & asthma --> 2 columns). Do numbers (1 for US, 2 for Spain) and for allergies/asthma variable. (In R)
      - Filter manifest in R also.
      - manifest is the rawest form of our data. All subsetting, and then qiime!
5) Our aims currently are technique-based not biology-based. Our aim: Data processing through qiime2, what we should have as aim: how does xxx affect yyy microbiome? Our research question: looking at regional differences first (differences in microbiome based on region (regions as variable), and then two separate aims (asthma and allergies), we should also include microbiome and diversity in our aims
   - (removed) aim: region as variable (remove this as not as relevant to our research question)  
   - 1st aim: country and allergies: yes or no, split to 8 groups, all at once
   - 2nd aim: country and asthma: yes or no (same as aim 2)
   - 3rd aim: relation between asthma and allergies (occurring together as stated by literatures)
      - another column combining asthma and allergies (1,1)
   - Comparing the microbiome of people with and without allergies. Disease markers? or exploratory question or drug targets? Diversity of microbes in each country and asthma and allergies. Characterizing microbiomes, different diets/culture/cuisine for what region is associated with less allergies or asthma?  
6) We filtered the metadata (only the healthy,so ~400), if statement or for loop to keep whatever is in the metadata. Sample IDs match. We have to subset the manifest file first.
7) Even though we have different sample sizes, just go with it (statistics takes it into account)
8) If we still need to wrangle, we can change proposal requirements
9) Next meeting: Oct 15, Tuesday 2pm


 ## Action Items
1) Background research and writing the introduction
2) Subsetting the manifest file and Convert strings to numbers
3) Gantt chart
