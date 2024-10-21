# October 21 - 6th meeting

## Attendees
- Andrea Garcia
- Yna Ortiz
- Julia Jung
- Chaeyoon Chang
- Vivian Tan 

## Agenda
1) Double check chosen sampling depth --> sample-id category was not visible when rarefaction plot was generated
2) Confirm if response variables are good
3) Confirm if we can do both Bray-Curtis and Weighted Unifrac --> we noticed in the paper we got the dataset from that they used both
4) Confirm if the PCoA plot is the best option
5) Double check that the biological relevance is enough --> we brought up the difference in diet between the regions
6) Double check that difference in sample size is okay --> 4 and 5 for Spain, 1 in both country and allergies/asthma
   
## Minutes
Sampling Depth/Alpha rarefraction:
- no metadata for sample ID, to look at each individual sample
- it's possible the input file was geneated differently, and "sample ID" can't be read as a metadata column.. or not categorical. (Sam will check)
- But looking at graph will still be fine, as we're still looking at the plateau.

Sample Size:
- our limiting group was for Spain, will it be an issue for disparity between samples for each group?
- It's something we can discuss in discussion as a limitation, but we should be concern where there's a different between Spain
- It's likely that what it will not be significant, due to the difference in sample size. We can say we saw a trend, but it's possible it's not a significance...

Bray-Curtis and Weighted Unifrac:
- it should be fine if other papers have done both beta diversity analysis.

PCoA plot:
- better for seeing correlations, thus we decided on PCoA plot than box plot.
- Better to put all countries together into a PCoA plot to see how they cluster, to add to our indicator species experimental aim #3 discussion/results.
- Include all the different countries (0,1) into a PCoA plot, one for asthma and one for allergy . But may take alot of time to process, due to the amount of data in one graph.

Biological relevance: for regional differences, can we only focus on dietary factors?
- In discussion, we may still have to focus on different aspects, but it works for a hypothesis.
- For experimental aim #3, we can see if the biomarkers + higher UPF consumption countries Ex. more diversity in mediterrean diet + Western diet
- So for example, in PCoA, lets say we may see a bigger different in those with/without allergies in USA, in compared to UK (with/without allergies). These results could be explained with diet and literature, but we also have to brief discuss that there's other reasons.
- We could also do other variables with weight and age to control it. If good distriubtion (ex. normal distriubtion, bell-shaped), don't need to control for the variable. If it's skewed, we should mention it and fix the data -- thus to only have "diet" as a variable. (ex. lots of 60 year olds, and very few 20 year olds -- then we'd remove the 20 years old.)
- Be mindful of what your data represents. If there's two huge group distriubtion, we should just decide on one population and discuss why it's more important -- likely choosing the bigger sample size in the end.
- The biggest hurdle is controlling for genetic differences, as we'd want a population with very close genetics. But it's still valid to say of why we're doing this comparison using diet regardless of genetics -- we just have to make really strong comparisons. Ex. That people's diets in USA is very different from Argentina.
- Age distribution, should be done for each country (we'd want equal age distriubtion across different countries). 
- -----> if the controls AREN'T done, we could talk about it in the discussion. 


Dataset removal:
- we can talk about how Spain had too few samples, so we had to remove it. In [Dataset section]

If too  difficult:
- we could compare within a country if it's too broad of a comparison. And if we DO see a significance in oure preliminary results cross-country, we could go back to cross-country comparison.
- Focus within a country, or we'd have to control as many variables as we need.
- Likely to not talk about dietary differences, due to so many confounding variables of a country and there isn't a section for a "diet".
- ----->If we make the differentiation clear, it should be okay...


FINAL CONCLUSION:

How to organize dataset:
So we should just look at beta diversiy of between countries rather than having a healthy control group and asthma. In addition we can group Spain/Argentina + USA/UK together for UPFs consumption. Then we don't have to do any controls, as it will control for genetics if we combine countries together.

Keypoints: 
- Emphasize UPF consumption and diet. 
- For indicator species, we should still look at high and low UPF indicator species. I
- We can also do DESeq/core microbiome for aim experimental #3.

  
 ## Action Items
 Focus on combining low and high UPFs countries in experimental aim and data set. 
 Work on revision of proposal, likely given a week. 

note:
 Meeting moved to Wed 3PM. 
