# November 18 - 9th meeting

## Attendees
- Yna Ortiz
- Julia Jung
- Chaeyoon Chang
- Vivian Tan 

## Agenda
1) Discuss the results from the beta diversity analyses (QIIME2)
2) PCoA plots: discuss results, ask about the plot for all the different countries
4) Discuss the PERMANOVA results
5) Ask about instructions for presentation
   
## Minutes

PCoA:
Asthma shows a difference. [Bray Curtis & Jaccard] --> can see this on PCoA and statistics. 
Dr. Sun likes the 3 channels for a PCoA, but it's hard to gauge with the current PCoA.
- We should do a separate colour one for UPS, Asthma, Allergies for better visualizing to portray the stats.
Weighted Unifrac vs Unweighted Unifrac: showss that abundance is driving major differences.

Weighted Unifrac PCoA:
Play around with the axes, to expand the important parts of PCoA (just note that 3 samples are excluded from the image.  
- Ex. Log the axes 
Add ellipses, to show different groups 
- Stats ellipses will automatically do it.


Fixing visualization:
It's likely that only Weighted Unifrac is going to be in there, so we will try to fix that only. 
Dr. Sun advised us to get out analysis done, and decide which other figures will actually go in the manuscript before fixing some of the visualization. 


PERMANOVA:
Bray curtis and weighted unifac show statistical significance. 
The UPF status is significant for all of them, for every table. 

UPS status + asthma = significant
asthma by itself = not significant.
- So more processed food consumption, asthma has an effect on their microbiome / or the other way around. 
- So one of them is driving a difference in microbiome.
- Status of country can affect asthma.

This is only seen in asthma, not allergies. 
- Allergies is signiicant in some cases (ex. weighted unifrac) by itself. But not when combined with high/low UPF. 

Experimental Aim 3: 
We should just move forward with just asthma. Because asthma shows a more clear correlation, we will move forward with it. [describe it in the paper]
We will likely not get much from an allergies analysis. 


Important: Alpha diversity
Beta diversity is showing us there is significance, but it doesn't tell us "how" --> therefore alpha diversity is important. 
- Make two phyloseq for asthma: high and low UPF.
- then run alpha diversity for core microbiome, DESeq, etc.
- then cross-compare (dr. sun said "lay them together")
Ex. If UPF low or UPF high affects the microbiome.

DESeq2: 
- Can give false positives. 
   - Ex. 0 in group A, 100 in GroupB. DESeq will say that group A will have 10X more when in reality there is none in GroupA
ANCON:
- Better if we have more 0s in groups. 


Countries:
- It's better if we don't consider the countries. Dr. Sun just suggests keeping to the high UPF and low UPF. 

Presentation:
If another group is presenting your project, prepare something for coaching the other group. 
- Maybe if there's something that's unique in your analysis, you can teach them more about that. 
- You can coach them on whatever you'd like --> in the end they will decide on what to present. 
Feedback --> It will give us an idea on how to fix our visual, and narrative. [how clear the story is]


Manuscript:
Main figures will be our "skeletons" --> ex. 5 figures 
1st figure will likely be the weighted unifrac PCoA plot. 





## Action Items
We can put together skeleton for presentation or manuscript starting today if we'd like. 
Fix the Weighed UniFrac PCoA only. 
Next meeting: We can determine the key figures in our manuscript/presentation. 





