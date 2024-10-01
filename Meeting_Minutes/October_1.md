# October 1 - 2nd meeting

## Attendees
* Dr. Sun
* Sam
* Andrea Garcia
* Vivian Tan
* Yna Ortiz
* Julia Jung
* Chaeyoon Chang

## Agenda
1) Discuss possible research questions
   - Possible research questions: different regions and...
     - asthma
     - allergy_specific: pollen/season/hay_fever --> respiratory or environmental allergies
     - diet_specific_needs: lactose_intolerant, vegan + vegetarian
3) Ask if it is possible to clean allergy_specific and diet_special_needs columns
5) Access to UJEMI --> Had issues on September 30, need to check if people have explored this research question previously
6) How small can we go with sample size?
7) How much data processing is expected for the project proposal? --> What steps should we expect to take a long time?
8) Do we have to specify in our rationale why all of our samples come from an MS population?
  
## Minutes
#### Research questions
- asthma
  - 5 regions --> Might be easier to use countries instead of regions
  - Pros: easy to pick out in terms of the data
  - Cons: 173 samples with asthma across the whole data set --> Not a bad number, will have to trim down anyways because the dataset is so big
    - Metadata will have to be filtered as well as manifest initally --> Will make processing faster
- allergy_specific
  - Pros: Can be sold as more of theme with asthma
  - Cons: column is not standardized --> might be harder to clean the data initially
    - To make it simpler, make it binary --> If they have allergies or not 
- diet_specific_needs
  - Not a strong choice
       - Asthma and allergy can be sold more as a theme
   
#### Dataset rationale
- Well annotated data set --> rationale
- Use only the healthy samples

#### Aim 1: Data wrangling 
- Filtering metadata and manifest
- Compare regions and then add in asthma stuff later
- Create 2 new metadata categories: Combine regions and asthma, then regions and allergies
- Do this through R --> Do not use Excel!
  - Ask R to go into categories --> ex. US, yes and US, no (2 columns: for allergies and asthma) 
- Convert allergies column into a binary column
- New metadata and manifest file
- Both allergies and asthma would be good --> Higher chance that there is a significant difference and is more interesting
- Metadata: sample_id, region/asthma, region/allergies --> Keep other metadata columns, just make new ones
  - Also: changing allergies column to yes or no 

#### Ideal minimum sample size
- The dataset is pretty big --> Will probably be in the 100s, or a little less

#### Expectations for project proposal 
- Project is data wrangling heavy --> Can write a section on data wrangling instead of processing
- Can be adjusted because so much overhead is needed

#### Ideas for aims
- Aim 2: Data processing through QIIME2
- Aim 3: Do diversity metrics only on region --> Taxonomic analysis
  - See if whether regionally, there is any differences before bringing asthma + allergies in
  - Ex. US has more diverse microbiomes --> Will inform last aim 
- Aim 4: Bring in allergies and asthma --> Can pick our own analysis (diversity metrics, others we will learn in future R modules)
- Project might have to pivot depending on results
- Possible 5th aim: Functional analysis
- Most relevant analysis: Indicator taxa --> unique taxa to a cohort
 
 ## Action Items
* Start doing data wrangling by next week! --> So we can get as much help as possible 

* Send email to Dr. Sun to get link to UJEMI articles
   - But question is novel!
   - Get link for MS/eczema article
   - Check novelty in broader literature --> Novel in the sense it has not been done on this dataset
   - Can just check if our results matches that of the broader literature 
