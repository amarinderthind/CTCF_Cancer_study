# Study Title: CTCF Binding Site Mutations: Linking Topologically Associated Domains Dysregulation to Cutaneous Squamous Cell Carcinoma Progression 

**Background**:
Cutaneous squamous cell carcinoma (cSCC) is the most common lethal malignancy with metastatic potential. The high mutational burden in cSCC has made it difficult to under-stand the significance of variants in the noncoding and regulatory genome. This study presents the first investigation of mutations at CCCTC-binding factor binding sites (CTCFbs) of topologically associating domains (TADs) across defined stages of disease progression —primary tumours that have not metastasized, primary tumours that have metastasized, and lymph node metastasis

This repository contains R scripts used in the analysis of CTCF binding site mutation for its role in cancer progression.

1) **Loop-genes-association.r** is used to explore the association between mutations in CTCF binding sites (CTCFbs) and changes in expression for genes located within TAD loops under investigation.
2) **various-Genomic-regions-Backgroud-rate-comparison.R** : Comparison of mutation density between genomic regions and their background using various scaling factors.
3) **ctcfmut-vs-gene-expression-permutations.R** : This script performs a permutation test to evaluate whether the observed association between CTCF loop mutations and gene expression is statistically significant by comparing against a null distribution generated over 1,000 randomized iterations.
4) ** cancer-progression-analysis-PNM-vs-PM.r** : This script analyzes the association between CTCF binding site mutations and cSCC progression by comparing mutation frequencies across tumor subgroups and identifying differentially mutated loops and associated gene expression changes.
