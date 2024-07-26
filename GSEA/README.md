## Differential Expression Analysis
Differential expression analysis is used to identify changes in gene expression levels between different conditions or groups (e.g., treated vs. control). This analysis is crucial in understanding the biological differences and underlying mechanisms involved in various experimental conditions.

### Differential Expression Analysis using DESeq2

This repository provides an example of differential expression analysis using the DESeq2 package in R. The example uses the RNAseq data obtained from https://www.ebi.ac.uk/gxa/home. This data consists of 54 samples from 18 individuals. Each individual has a primary colorectal cancer sample, a metastatic liver sample, and a normal sample of the surrounding colonic epithilium. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts.

## Key Steps

- **Step 1: Install and Load Necessary Libraries**
  - Install the DESeq2 and other required packages.
  - Load the DESeq2, and ggplot2 libraries.

- **Step 2: Load and Prepare the Data**
  - Load the  dataset.
  - Extract the count data and sample information.

- **Step 3: Create DESeq2 Dataset**
  - Create a DESeq2 dataset object from the count data and sample information.

- **Step 4: Pre-filtering**
  - Remove genes with low counts to improve the power of the analysis.

- **Step 5: Differential Expression Analysis**
  - Run the DESeq function to perform differential expression analysis.
  - View the results.

- **Step 6: Visualize the Results**
  - Create a volcano plot to visualize the log2 fold changes and p-values.
  - Create an MA plot to visualize the mean expression and log2 fold changes.

- **Step 7: Extract Significant Results**
  - Identify significantly differentially expressed genes using a cutoff for adjusted p-value.
 
A complete Jupyter notebook can be found here: https://github.com/BhadraNivedita/Differential-gene-Expression-Analyses/blob/main/Differential%20gene%20Expression%20Analysis%20with%20DESEQ2.ipynb
