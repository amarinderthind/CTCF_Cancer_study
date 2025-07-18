## CTCF Binding Site Mutations: Linking Topologically Associated Domains Dysregulation to Cutaneous Squamous Cell Carcinoma Progression

Preprint : https://www.researchsquare.com/article/rs-6844715/v1

This repository contains R scripts used in the analysis of CTCF binding site mutation for its role in cancer progression.

**Background**:
The high mutational burden in Cutaneous squamous cell carcinoma (cSCC) has made it difficult to under-stand the significance of variants in the noncoding and regulatory genome. This study presents the first investigation of mutations at CCCTC-binding factor binding sites (CTCFbs) of topologically associating domains (TADs) across defined stages of disease progression.
 
 <p align="center">
<img src="https://github.com/user-attachments/assets/52891726-edc7-41a4-bb34-6b6835643466" width=86% height=700>&nbsp; &nbsp; 
</p>

**Description of the workflow (See figure above):**

1) `ctcf-gene-extarctions-both-inside-outside.r`  ctcf loop genes extraction both inside and outside loop upto 1000bp either side with ensemble-id.  

2) `various-Genomic-regions-Backgroud-rate-comparison.R`  Comparison of mutation density between genomic regions and their background using various scaling factors.  
3)  `Mutatinal-Density-overlapping-regions.R` Compare mutational Density between overlapping regions.
4)  `motif-CTCFbs-specific-plots.R`  Script for ploting CTCFbs (normalized position) position specific mutational density (average across chorot).
5) `Loop-genes-association.r` is used to explore the association between mutations in CTCF binding sites (CTCFbs) and changes in expression for genes located within TAD loops under investigation.  
6) `ctcfmut-vs-gene-expression-permutations.R` This script performs a permutation test to evaluate whether the observed association between CTCF loop mutations and gene expression is statistically significant by comparing against a null distribution generated over 1,000 randomized iterations.  

     **Note: Choose numCores based on available CPUs (e.g., use 75% of total to avoid system freeze; on 32-core machine, use 26).
     Set `chunk_size` (e.g., 48) to balance load — smaller chunks reduce memory per core, larger chunks may run faster if RAM allows.**
 
7)  `cancer-progression-analysis-PNM-vs-PM.R` This script analyzes the association between CTCF binding site mutations and cSCC progression by comparing mutation frequencies across tumor subgroups and identifying differentially mutated loops and associated gene expression changes.
8) `Cor_AUC_analysis-PNM-PM_loops-mutated_vs_other_within_a_group.r` This script evaluates the association between gene expression and CTCF loop mutation status using Spearman correlation and AUC to identify genes whose expression distinguishes mutated from non-mutated loops.
   
#### Note: Due to ethical reasons,Somatic mutation files relevant to this study are available upon request from the corresponding author.

```
## License

This pipeline code is released under an academic non-commercial use license.  
Use is permitted for academic research with proper citation.  
Commercial use or redistribution (even modified versions) requires written permission.

See the [LICENSE](./LICENSE) file for details.
```
