# November 4 - 8th meeting

## Attendees
- Andrea Garcia
- Yna Ortiz
- Julia Jung
- Chaeyoon Chang
- Vivian Tan 

## Agenda
1) Rarefaction curve --> generated based on sample-id, how many samples is enough?
2) Generated unweighted UniFrac visualization --> significance, check code
3) Beta diversity --> should it have been filtered based on metadata?
4) Server problems --> set new timeline based on this?
5) Lab notebook code for QIIME
   
## Minutes
Rarefaction curve
- duplicating column didn't let us see all the samples --> changed the sampling depth (8000)
- based on UPF for allergies: 26 samples, sampling depth lessens by ~4300 (more samples)
- 20 samples at 8000 sampling depth which is not bad
- reducing ASVs will reduce what we get from our analysis (unless we get no significance from the 20 samples)
- based on UPF for asthma: even less, but still ok at 8000
- if we keep sampling depth we see on the rarefaction curve (8000) we lose 20-30% of our samples
- at 6000 sampling depth: not bad, losing 2% of features but close to 90% of samples
- for allergies: 6000 sampling depth
- for asthma: 6000 sampling depth
- toggle until we get a significant drop in ASVs - we want to get a specific number for the sampling depth
- justification: retaining as many samples per group without a significant drop in the features
- not sure if teaching team is gonna look through the code to check how good the sampling depth is

Unweighted UniFrac - allergies
- significance is present!
- still need to generate the other ones
- need to plot to see the clustering to see if it is in line with our hypothesis
- done at a sampling depth of 6000
- need to generate PCoA
- if we find that one metric looks better than the other, do we continue using it? - do all for both
- not sure how long it took to generate since it was run overnight but visualizing took ~20 mins

Checking code for beta diversity
- are we supposed to be filtering based on metadata before we run beta diversity metrics? - yes
- generate table --> frequency based
- were we supposed to create a new table based on filtered metadata? - as long as we specify it should be fine

Server problems
- server won't be gone but if things go wrong, Evelyn won't be able to help us troubleshoot
- whole class might be on server so it could wrong very slow
- QIIME can be downloaded onto computer and everything can be done locally but need to download all the files

Lab notebook code
- just copy-paste code from QIIME2

Others
- DESeq: DESeq2 is used for RNA seq analysis rather than 16S --> creates a lot of false positives
- ANCOM recommended to use instead of DESeq (found on QIIME and there is a tutorial on how to use it on the QIIME website)
- running ANCOM should take the same amount of time as running DESeq (might take ~30 hrs for 100 samples </3)
- Evelyn will join our meetings after reading week

## Action Items
- make sure .qza and .qzv files are uploaded on GitHub (teaching team might open them)
- generate PCoA plot/s
- do the rest of the beta diversity analyses
- email Sam if we have questions
- likely no meeting during reading week 
