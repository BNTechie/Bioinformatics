# Self learning resources

## Tools for bioinformatics analysis

1. Gene ID conversion tool: https://www.syngoportal.org/convert
2. Connectivity Map analysis tool: https://clue.io/query
3. Gprofiler: https://biit.cs.ut.ee/gprofiler/gost
4. Gene set enrichment tool: http://bioinformatics.sdstate.edu/go/
5. Pantherdb: https://pantherdb.org/tools/compareToRefList.jsp?&showAll=false

## GWAS studies practical guide
1. https://www.youtube.com/watch?v=nrbgly0Bcv8
2. https://www.r-bloggers.com/2017/10/genome-wide-association-studies-in-r/

### What is GWAS study? 

A Genome-Wide Association Study (GWAS) is a research approach used to identify genetic variations associated with specific diseases or traits. Here's a detailed explanation of what GWAS involves and its significance:

### What is GWAS?

1. **Objective:**
   - The primary goal of a GWAS is to uncover the genetic basis of complex traits or diseases by scanning the genome for single nucleotide polymorphisms (SNPs) that occur more frequently in individuals with a particular condition compared to those without.

2. **Methodology:**
   - **Sample Collection:** GWAS begins with the collection of DNA samples from two groups: individuals with the disease or trait of interest (cases) and individuals without it (controls).
   - **Genotyping:** The DNA samples are genotyped to identify SNPs across the genome. Modern GWAS typically use high-throughput genotyping arrays that can examine hundreds of thousands to millions of SNPs simultaneously.
   - **Statistical Analysis:** Each SNP is statistically analyzed to determine if there is a significant association between the SNP and the disease or trait. This involves comparing the frequency of each SNP in cases versus controls.
   - **Correction for Multiple Testing:** Given the large number of SNPs tested, corrections for multiple comparisons are necessary to reduce the likelihood of false positives. Common methods include the Bonferroni correction or the False Discovery Rate (FDR) approach.
   - **Replication:** Findings from the initial analysis are often validated in independent cohorts to confirm the associations.

3. **Output:**
   - The results of a GWAS are typically presented as a Manhattan plot, where each dot represents a SNP and its association with the trait. Peaks in the plot indicate regions of the genome that are significantly associated with the trait.

### Significance of GWAS

1. **Understanding Genetic Architecture:**
   - GWAS has helped identify numerous genetic loci associated with a wide range of diseases and traits, providing insights into their genetic architecture and biological pathways.

2. **Disease Mechanisms:**
   - By pinpointing genetic variants linked to diseases, GWAS can reveal new biological mechanisms and pathways involved in disease development, which can inform the development of new therapeutic targets.

3. **Personalized Medicine:**
   - GWAS findings contribute to the field of personalized medicine by identifying genetic markers that can predict disease risk, treatment response, or adverse drug reactions, allowing for more tailored healthcare strategies.

4. **Polygenic Risk Scores:**
   - GWAS data can be used to create polygenic risk scores, which aggregate the effects of multiple genetic variants to estimate an individual's genetic predisposition to a particular disease.

### Challenges and Limitations

1. **Complex Traits:**
   - Many complex traits and diseases are influenced by numerous genetic variants, each contributing a small effect, as well as environmental factors. This makes it challenging to identify all relevant variants.

2. **Population Stratification:**
   - Genetic differences between populations can lead to spurious associations if not properly controlled for, making it crucial to include diverse populations in GWAS to ensure findings are broadly applicable.

3. **Missing Heritability:**
   - Despite identifying many genetic associations, a large portion of the heritability of complex traits remains unexplained. This "missing heritability" suggests that other factors, such as rare variants, gene-gene interactions, and gene-environment interactions, also play significant roles.

### Conclusion

GWAS is a powerful tool for uncovering the genetic underpinnings of diseases and traits. By scanning the genome for associations between genetic variants and specific conditions, GWAS has significantly advanced our understanding of human genetics and contributed to the development of personalized medicine. However, it also faces challenges that require ongoing research and methodological improvements to fully realize its potential.

About PLINK: https://www.cog-genomics.org/plink/2.0/input


