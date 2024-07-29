# HTSeq-counts

To work with HTSeq counts, you generally need to process raw sequencing data (FASTQ files) to generate alignment files (BAM files) and then use HTSeq to count the reads aligned to each gene. This involves several steps including quality control, alignment, and counting. Below is a step-by-step guide to process a public dataset using common tools on a Mac.

### Step-by-Step Guide

### Prerequisites

1. **Install Homebrew**: If you haven't already, install [Homebrew](https://brew.sh/), the package manager for macOS.
   
   ```sh
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

2. **Install Required Tools**: You will need to install `fastqc`, `hisat2`, `samtools`, and `htseq`.

   ```sh
   brew install fastqc
   brew install hisat2
   brew install samtools
   brew install htseq
   ```

### Step 1: Download Example Data

For this example, we will use a small public dataset from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra). Let's download a FASTQ file from an example study.

```sh
# Create a directory for the data
mkdir -p ~/htseq_example
cd ~/htseq_example

# Download example FASTQ file using SRA Toolkit
brew install sratoolkit
fastq-dump --split-files SRR585570
```

### Step 2: Quality Control with FastQC

Run FastQC to check the quality of the raw sequencing data.

```sh
fastqc SRR585570_1.fastq SRR585570_2.fastq
```

### Step 3: Align Reads to the Reference Genome with HISAT2

Download the reference genome and build the HISAT2 index.

```sh
# Download reference genome (e.g., GRCh38)
curl -O ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build HISAT2 index
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_index
```

Align the reads to the reference genome.

```sh
hisat2 -x GRCh38_index -1 SRR585570_1.fastq -2 SRR585570_2.fastq -S SRR585570.sam
```

### Step 4: Convert SAM to BAM and Sort

Convert the SAM file to a BAM file and sort it.

```sh
samtools view -S -b SRR585570.sam > SRR585570.bam
samtools sort SRR585570.bam -o SRR585570_sorted.bam
samtools index SRR585570_sorted.bam
```

### Step 5: Count Reads with HTSeq

Download the GTF file for the reference genome (annotation file).

```sh
curl -O ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gunzip Homo_sapiens.GRCh38.104.gtf.gz
```

Count the reads aligned to each gene.

```sh
htseq-count -f bam -r pos -s no -t exon -i gene_id SRR585570_sorted.bam Homo_sapiens.GRCh38.104.gtf > counts.txt
```

### Example Script

Here's a complete script to automate the above steps:

```sh
# Create a directory for the data
mkdir -p ~/htseq_example
cd ~/htseq_example

# Download example FASTQ file using SRA Toolkit
brew install sratoolkit
fastq-dump --split-files SRR585570

# Run FastQC
brew install fastqc
fastqc SRR585570_1.fastq SRR585570_2.fastq

# Download reference genome and build HISAT2 index
curl -O ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
brew install hisat2
hisat2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_index

# Align reads to the reference genome
hisat2 -x GRCh38_index -1 SRR585570_1.fastq -2 SRR585570_2.fastq -S SRR585570.sam

# Convert SAM to BAM and sort
brew install samtools
samtools view -S -b SRR585570.sam > SRR585570.bam
samtools sort SRR585570.bam -o SRR585570_sorted.bam
samtools index SRR585570_sorted.bam

# Download the GTF file for the reference genome
curl -O ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gunzip Homo_sapiens.GRCh38.104.gtf.gz

# Count reads with HTSeq
brew install htseq
htseq-count -f bam -r pos -s no -t exon -i gene_id SRR585570_sorted.bam Homo_sapiens.GRCh38.104.gtf > counts.txt
```

### Running the Script

1. Open a terminal on your Mac.
2. Copy and paste the script into the terminal.
3. Replace `SRR585570` with the appropriate SRA accession number if using a different dataset.
4. Run the script to download the data, perform quality control, align reads, and count the reads using HTSeq.

This script provides a comprehensive workflow for processing RNA-seq data and generating HTSeq counts using public datasets on a Mac.
