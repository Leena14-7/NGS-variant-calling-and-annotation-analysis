# NGS-variant-calling-and-annotation-analysis

# Description
Variant calling analysis of whole genome sequencing (WGS) data from Bioproject (tongue cancer samples and cell lines) using,
   
   . Tools - GATK , BWA, Samtools, SRA Toolkit, bcftools, VEP and IGV   
   . Reference genome-: hg38(chromosomes 5 and 7)

# Workflow
1.Raw Data is download from public data set by using SRA Toolkit
2.Quality check by using Fastqc
3.Trimming by using Fastp
4.Aligmnet with BWA-MEM 
5.Processed BAM files using Samtools
6.Variant calling with GATK
  (HaplotypeCaller/Mutect2)
7.Filtered variants to retain high-quality results  
8.Annotated variants using VEP  
9.Visualized selected variants in IGV

# Results
- Chromosome 5 showed a higher number of variants compared to chromosome 7.  
- Most variants were located in intergenic and intronic regions.  
- These variants are predicted to have minimal impact on protein function.  
- The majority were classified as having modifier impact.

# Author 
  Leena Patil 
  NGS-variant-calling-and-annotation-analysis
