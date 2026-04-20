#!/bin/bash

echo "===== NGS VARIANT CALLING PIPELINE START ====="

# Project Info
echo "Project: PRJEB62494 (Tongue cancer)"
echo "Reference: GRCh38 (chr5 & chr7)"

# data selection 
https://www.ncbi.nlm.nih.gov/sra/?term=PRJEB62494

# ---------------------------------------------------------
# System Update
# ---------------------------------------------------------
echo "Updating system..."
sudo apt-get update && sudo apt-get upgrade -y

# ---------------------------------------------------------
# Create Directories
# ---------------------------------------------------------
echo "Creating directories..."
mkdir -p \
   raw_data \
   trimmed \
   qc_reports \
   reference_genome \
   alignment_reads \
   variants 
          
# -------------------------------------------------------------------
# Install Tools
# -------------------------------------------------------------------
echo "Installing tools..."
sudo apt-get install -y \
   sra-toolkit \
   fastqc \
   fastp \
   bwa \
   samtools \
   bcftools \
   tabix
         
# ------------------------------------------------------------------------
# Step 1: Download Data
# ------------------------------------------------------------------------
echo "Downloading SRA data..."         

# Configure SRA toolkit to save in the current directory 
vdb-config --interactive

fastq-dump --split-files --gzip -X 1000000 ERR11468775
fastq-dump --split-files --gzip -X 1000000 ERR11468776
fastq-dump --split-files --gzip -X 1000000 ERR11468777

echo "Download complete"


# -------------------------------------------------------------------------
# Step 2: Quality Check
# -------------------------------------------------------------------------
echo "Running FastQC (raw)..."

fastqc raw_data/*.fastq.gz -o qc_reports/


echo "Quality check complete"


# -------------------------------------------------------------------------
# Step 3: Trimming
# --------------------------------------------------------------------------
echo "Running fastp..."

for FILE in raw_data/*_1.fastq.gz
do
  SAMPLE=$(basename ${FILE} _1.fastq.gz)

  fastp \
     -i raw_data/${SAMPLE}_1.fastq.gz \
     -I raw_data/${SAMPLE}_2.fastq.gz \
     -o trimmed/${SAMPLE}_trim_1.fastq.gz \
     -O trimmed/${SAMPLE}_trim_2.fastq.gz \
     --detect_adapter_for_pe \
     --thread 4 \
     --html trimmed/${SAMPLE}.html
done

echo "Trimming complete"


# ------------------------------------------------------------------------
# Step 4: QC After Trim
# ------------------------------------------------------------------------
echo "Running FastQC (trimmed)..."

fastqc trimmed/*_trim_*.fastq.gz o-trimmed/


# ------------------------------------------------------------------------
# Step 5: Reference Genome
# ------------------------------------------------------------------------
echo "Downloading reference genome..."

wget -p reference_genome/https://hgdownload.gi.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz
wget -p reference_genome/https://hgdownload.gi.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz

gunzip reference_genome/*.fa.gz

cat reference_genome/chr5.fa reference_genome/chr7.fa > reference_genome/ref_chr5_chr7.fa 

echo "Indexing reference..."
bwa index -a bwtsw reference_genome/ref_chr5_chr7.fa
samtools faidx reference_genome/ref_chr5_chr7.fa


# -------------------------------------------------------------------------
# Step 6: Alignment
# -------------------------------------------------------------------------
echo "Running BWA alignment..."

for FILE in trimmed/*_trim_1.fastq.gz
do
  SAMPLE=$(basename "$FILE" _trim_1.fastq.gz)

  bwa mem -t 4 reference_genome/ref_chr5_chr7.fa \
    trimmed/${SAMPLE}_trim_1.fastq.gz \
    trimmed/${SAMPLE}_trim_2.fastq.gz \
    > alignment_reads/${SAMPLE}.sam
done

echo "Alignment complete"

# ------------------------------------------------------------------------
# Step 7: SAM → BAM
# -----------------------------------------------------------------------
echo "Converting SAM to BAM..."

for FILE in alignment_reads/*.sam
do
  SAMPLE=$(basenmae "$FILE".sam)

  samtools view -b "$FILE".sam > alignment_reads/${SAMPLE}.bam
  
  samtools sort alignment_reads/${SAMPLE}.bam -o alignment_reads/${SAMPLE}_sorted.bam

done

echo "BAM processing complete"

# ----------------------------------------------------------------------
# Step 8: GATK Preprocessing
# ---------------------------------------------------------------------
echo "GATK Preprocessing..."

gatk CreateSequenceDictionary \
   -R reference_genome/ref_chr5_chr7.fa \
   -O reference_genome/ref_chr5_chr7.dict

# Download the sites for BQSR having known variants 
wget -P reference_genome/ https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

# docker version to run gatk latest for version control
docker pull broadinstitute/gatk:latest
docker run -it -v $PWD:/data broadinstitute/gatk:latest
 cd\variant_calling_project\alignment_reads

# Loop for samples
for BAM in *_sorted.bam; do
  SAMPLE=${BAM%_sorted.bam}

  # A : Add Read Groups
  gatk AddOrReplaceReadGroups \
    -I "$BAM" \
    -O "${SAMPLE}_RG.bam" \
    -RGID "${SAMPLE}" \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM "${SAMPLE}"

  # B : MarkDuplicates
  gatk MarkDuplicates \
    -I "${SAMPLE}_RG.bam" \
    -O "${SAMPLE}_Markdup.bam" \
    -M "${SAMPLE}_metrics.txt"

  # C : Index BAM
  samtools index "${SAMPLE}_Markdup.bam"

  # D : BaseRecalibrator
  gatk BaseRecalibrator \
    -I "${SAMPLE}_Markdup.bam" \
    -R /data/reference_genome/ref_chr5_chr7.fa \
    --known-sites /data/reference_genome/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O "${SAMPLE}_recal.table"

  # E : ApplyBQSR
  gatk ApplyBQSR \
    -R /data/reference_genome/ref_chr5_chr7.fa \
    -I "${SAMPLE}_Markdup.bam" \
    --bqsr-recal-file "${SAMPLE}_recal.table" \
    -O "${SAMPLE}_recal.bam"

done

# F : Index bam file 
   samtools index *_recal.bam 

echo "GATK Preprocessing Complete"


# ----------------------------------------------------------------------------
# Step 9: Variant Calling
# ----------------------------------------------------------------------------
echo "Running variant calling..."

# somatic variant calling
for BAM in alignment_reads/*_recal.bam; do
  BASE=$(basename "$BAM".bam)

  gatk Mutect2 \
    -R reference_genome/ref_chr5_chr7.fa \
    -I "$BAM" \
    -O variants2/${BASE}_somatic.vcf.gz
done 

# Germline variant calling 
for BAM in alignment_reads/*_recal.bam; do
  BASE=$(basename "$BAM".bam)

  gatk Haplotypecaller \
     -R reference_genome/ref_chr5_chr7.fa \
     -I alignment_reads/sample1_recal.bam \
     -O variants/${BASE}sample1_germline.g.vcf.gz \
     -ERC GVCF 
done 

# Repeat for sample2 and sample3

echo "Variant calling complete"

# -----------------------------------------------------------------------------
# Step 10: Variant Filtering
# ----------------------------------------------------------------------------
echo "Running variant Filtering..."

for VCF in variants/*_somatic.vcf.gz; do
  BASE=$(basename $VCF _somatic.vcf.gz)

  ./gatk FilterMutectCalls \
    -R reference_genome/ref_chr5_chr7.fa \
    -V $VCF \
    -O variants/${BASE}_filtered.vcf.gz

   bcftools view \
    -f PASS \
    variants/${BASE}_filtered.vcf.gz \
    > variants/${BASE}_final.vcf.gz
done
 
# bigzip  
for VCF in variants/*_final.vcf; do
  BASE=$(basename $VCF .vcf)

  bgzip -c $VCF > variants/${BASE}.vcf.gz
done

#Indexing
for VCF in variants/*.vcf.gz; do
  tabix -f -p vcf $VCF
done

echo "Running variant Filtering complete"

# -----------------------------------------------------------------------
# Step 10: Split by Chromosome
# ----------------------------------------------------------------------
echo "Splitting variants by chromosome..."

for VCF in variants3/*_final.vcf.gz; do
  BASE=$(basename $VCF _final.vcf.gz)

  bcftools view -r chr5 $VCF -o variants/${BASE}_chr5.vcf.gz
  bcftools view -r chr7 $VCF -o variants/${BASE}_chr7.vcf.gz
done

# bgzip of chr5 and chr 7
for VCF in variants/*_chr5.vcf; do
  BASE=$(basename $VCF .vcf)

  bgzip -c $VCF > variants/${BASE}.vcf.gz
done

for VCF in variants/*_chr7.vcf; do
  BASE=$(basename $VCF .vcf)

  bgzip -c $VCF > variants2/${BASE}.vcf.gz
done

# Indexing of chr5 and chr 7
for VCF in variants/*_chr5.vcf.gz; do
  tabix -p vcf $VCF
done

for VCF in variants2/*_chr7.vcf.gz; do
  tabix -p vcf $VCF
done

echo "Splitting complete"

# -------------------------------------------------------------------
# Step 11: Merge
# -------------------------------------------------------------------
echo "Merging chromosome-wise VCFs..."

bcftools merge variants/*_chr5.vcf.gz -o variants/merge_chr5.vcf
bcftools merge variants/*_chr7.vcf.gz -o variants/merge_chr7.vcf

echo "Merging complete"

# -------------------------------------------------------------------
# END
# --------------------------------------------------------------------
echo "===== PIPELINE COMPLETED ====="




