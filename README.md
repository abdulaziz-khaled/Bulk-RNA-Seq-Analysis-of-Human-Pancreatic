# 🧬 Bulk RNA-seq Analysis Reveals Gene Expression changes in Human Pancreatic β-cells Under Metabolic Stress.
This project performs bulk RNA-Seq data analysis on human pancreatic β-cells under metabolic stress using dataset GSE159984. The pipeline includes data retrieval from GEO using SRA Toolkit, quality control with FastQC and MultiQC, read trimming with Trimmomatic, read alignment to the human reference genome (hg38) using STAR/HISAT2, gene quantification with featureCounts, normalization and differential expression analysis with DESeq2, and visualization (PCA, heatmap, and volcano plot). The analysis compares 3 control samples (paired-end) against 3 stressed samples (single-end). All tools used are open-source and run on Ubuntu via WSL (Windows Subsystem for Linux).

# 🧩 Table of Contents
- [Background](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#background)
- [Data Source](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#-data-source-and-sample-information)
- [WorkFlow Overview](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#%EF%B8%8F-workflow-overview)
- [Tools and Software Used](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#%EF%B8%8F-tools-and-software-used)
- [Results](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#-results)
- [Coding](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic?tab=readme-ov-file#-analysis-execution-code-snippets)

# 💡Background 
Next Generation Sequencing (NGS) technologies have transformed biological research by enabling the genome-wide investigation of gene expression. Among these, RNA sequencing (RNA-seq) has become the gold standard for transcriptomic profiling, providing insights into the regulation of cellular processes under different physiological and pathological conditions.
Human pancreatic β-cells play a central role in maintaining glucose homeostasis by secreting insulin. Dysfunction of β-cells is a hallmark of type 2 diabetes mellitus (T2DM), a metabolic disorder with increasing global prevalence. Metabolic stress, such as prolonged exposure to high levels of fatty acids (e.g., palmitate) or glucose, can impair β-cell survival and function, contributing to disease progression.
The dataset analyzed in this project, GSE159984, investigates the transcriptomic response of human β-cells under metabolic stress compared to untreated controls. By applying Bulk RNA-seq analysis, it becomes possible to identify differentially expressed genes (DEGs) and uncover biological pathways that may underlie stress-induced β-cell dysfunction.
The overarching aim of this analysis is: To identify differentially expressed genes (DEGs) and enriched pathways in human pancreatic β-cells exposed to metabolic stress compared to control conditions.

## 📚 Data Source and Sample Information

| Field | Details |
| :--- | :--- |
| **Dataset Source** | **GSE159984 (GEO)** |
| **Sample Groups** | **3 Control** vs **3 Stress**  |
| **Control Group Accessions** | SRR12885688, SRR12885679, SRR12885709 |
| **Control Read Type** | Paired-end (**6 FASTQ files**) |
| **Stress Group Accessions** | SRR12885579, SRR12885580, SRR12885581 |
| **Stress Read Type** | Single-end (**3 FASTQ files**) |

# ⚙️ Workflow Overview
This section outlines the complete bioinformatics pipeline used for the RNA-seq analysis of the GSE159984 dataset, covering the process from raw data retrieval to functional biological interpretation.

**(1) 📥 Data Retrieval and Preparation**

The workflow begins by acquiring and preparing the data for primary analysis:

• **Dataset Source**: GSE159984 (GEO)

• **Action**: Downloading the raw FASTQ files for the Control (paired-end) and Stress (single-end) samples.

**(2) 🛡️ Quality Control (QC) & Trimming**

Data integrity is crucial. We enforce strict quality standards through these two steps:

1. **Initial Quality Control (QC)**
2. **Sequence Quality Trimming (Read Trimming)**

• **Goal**: Removing low-quality bases and adapter sequences.

• **Verification**: Running QC again on the trimmed reads to confirm high quality.

**(3) 🎯 Indexing & Read Alignment (Mapping)**

Clean reads are then mapped to the reference genome:

• **Indexing**: Creating a searchable index for the reference human genome (hg38) (e.g., using STAR or HISAT2).

• **Read Alignment**: Mapping (aligning) the high-quality reads against the indexed genome to determine their genomic origin.

**(4) 🔢 Quantification (Gene Expression Counts)**

Converting aligned reads into quantifiable expression levels:

• **Goal**: Counting the number of reads that map to each specific gene to obtain Raw Counts.

• **Output**: A count matrix detailing the expression level for every gene across all samples.

**(5) 📈 Differential Expression Analysis (DEGs)**

The core statistical analysis to find significant changes:

1. **Normalization**: Adjusting raw counts for sequencing depth biases.
2. **Dimension Reduction**: Performing PCA (Principal Component Analysis) to visualize sample clustering and overall variability.
3. **DEGs Analysis**: Identifying Differentially Expressed Genes (DEGs) between the Stress and Control groups using established packages (e.g., DESeq2).

**(6) 📑 Functional Enhancement & Biological Interpretation**

The final step to link statistical results back to biological meaning:

**Goal**: Understanding the biological pathways and functions associated with the identified DEGs.

## 🛠️ Tools and Software Used

The table below details the specific tools, their versions, and R packages utilized in each unique stage of the bioinformatics pipeline.

| Workflow Stage | Tool / R Package | Version | Primary Functions |
| :--- | :--- | :--- | :--- |
| **Data Retrieval** | **SRA-toolkit** | **3.2.1** | Downloading raw **FASTQ** files from the SRA database. |
| **Quality Control & Trimming** | **FastQC** | **v0.12.1** | Assessing initial read quality. |
| | **MultiQC** | **1.31** | Summarizing reports from multiple tools. |
| | **Trimmomatic** | **0.40** | Removing low-quality reads and adapters. |
| **Alignment (Mapping)** | **HISAT2** | **2.2.1** | Aligning quality-filtered reads to the human reference genome (**hg38**). |
| **Quantification** | **FeatureCounts** | **2.0.1** | Counting the number of reads mapped to each gene to generate **Raw Counts**. |
| **Differential Expression** | **DESeq2** | **1.42.0** | Performing normalization and statistical analysis to identify **DEGs**. |
| **Functional Enrichment** | **clusterProfiler** | **4.10.1** | Performing **GO** and **KEGG** enrichment analysis. |
| | **org.Hs.eg.db** | **3.18.0** | Providing human annotation data for gene ID conversion. |
| | **enrichplot** | **1.22.0** | Visualizing the results of enrichment analysis. |
| **Data Visualization** | **ggplot2** | **3.5.2** | Creating high-quality custom plots (Volcano, Box plots). |
| | **Pheatmap** | **1.0.13** | Generating visually appealing **Heatmaps**. |


# 📃 Results
**(1) Data Retrieval** 

[Screenshot](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/1a8afe130a5905bd35699d3f3e24b5832b0a2557/Screenshot%20from%202025-09-24%2019-00-53.png)

**(2) QC (Before Read Trimming)**
1. [SRR12885579_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885579_fastqc.html)
2. [SRR12885580_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885580_fastqc.html)
3. [SRR12885581_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885581_fastqc.html)
4. [SRR12885679_1_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885679_1_fastqc.html)
5. [SRR12885679_2_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885679_2_fastqc.html)
6. [SRR12885688_1_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885688_1_fastqc.html)
7. [SRR12885688_2_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885688_2_fastqc.html)
8. [SRR12885709_1_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885709_1_fastqc.html)
9. [SRR12885709_2_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/SRR12885709_2_fastqc.html)
10. [multiqc_report.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC/multiqc_report.html)

**(3) Read Trimming**

[Screenshot](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/41f2f955894c6471fd3aca812b47374e1c08e3da/Screenshot%20from%202025-09-25%2021-03-56.png)

**(4) QC (After Read Trimming)** 
1. [SRR12885579.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885579.trim_fastqc.html#M9)
2. [SRR12885580.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885580.trim_fastqc.html)
3. [SRR12885581.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885581.trim_fastqc.html)
4. [SRR12885679_1.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885679_1.trim_fastqc.html)
5. [SRR12885679_2.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885679_2.trim_fastqc.html)
6. [SRR12885688_1.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885688_1.trim_fastqc.html)
7. [SRR12885688_2.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885688_2.trim_fastqc.html)
8. [SRR12885709_1.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885709_1.trim_fastqc.html)
9. [SRR12885709_2.trim_fastqc.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/SRR12885709_2.trim_fastqc.html)
10. [multiqc_report.html](https://abdulaziz-khaled.github.io/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/QC_trimmed/multiqc_report.html)

**(5) Indexing and Alignment**
- [Genome Index](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/Screenshot%20from%202025-10-02%2007-41-36.png)
- [Alignment Single](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/Screenshot%20from%202025-10-02%2007-45-07.png)
- [Alignment Paired](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/Screenshot%20from%202025-10-02%2007-45-19.png)

**(6) Quantification (Gene Counts)**

[Count](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/counts)

**(7) Normalization and Downstream Analysis**
- [Normalization](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/results/normalized_counts.csv)
- [PCA](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/results/PCA_plot.pdf)
- [VSD](https://github.com/abdulaziz-khaled/Bulk-RNA-Seq-Analysis-of-Human-Pancreatic/blob/main/results/vsd_object.rds)
- [All Results](
- [Top 50 DEGs](

## 💻 Analysis Execution: Code Snippets

### 1. 📥 Data Retrieval (SRA-toolkit)

This step outlines setting up the project structure and downloading the raw sequencing data using the `SRA-toolkit` commands.

```bash
# 1. Setup Project Directory
# Create the main project folder and the raw data folder
mkdir -p Bulk_RNAseq_Project/raw_data
cd Bulk_RNAseq_Project/raw_data

# 2. Download SRA Files (prefetch)
# Download Paired-end (Control) samples
prefetch SRR12885688 SRR12885679 SRR12885709

# Download Single-end (Stress) samples
prefetch SRR12885579 SRR12885580 SRR12885581

# 3. Convert SRA to FASTQ (fasterq-dump)
# Convert SRA files to FASTQ format. 
# Note: --split-files is mandatory for paired-end reads.

# Convert Control files (Paired-end)
fasterq-dump --split-files SRR12885688
fasterq-dump --split-files SRR12885679
fasterq-dump --split-files SRR12885709

# Convert Stress files (Single-end)
fasterq-dump SRR12885579
fasterq-dump SRR12885580
fasterq-dump SRR12885581
```

### 2. 🛡️ Quality Control (QC) - Raw Data

This step performs an initial quality assessment on the raw FASTQ files using `FastQC` and generates a summarized report with `MultiQC`.

```bash
# Ensure we are in the raw_data directory where FASTQ files are located
cd ~/Bulk_RNAseq_Project/raw_data 

# Create the output directory for QC reports
mkdir -p ../QC

# Run FastQC on all FASTQ files
# -t 4: Use 4 CPU threads for faster processing
# -o ../QC: Output reports to the QC directory
fastqc -t 4 -o ../QC *.fastq

# Change directory to the output folder
cd ../QC

# Aggregate all individual FastQC reports into a single, summarized HTML report
multiqc . -o.
```
### 3. ✂️ Read Trimming

This step trims adapters and low-quality bases from the raw FASTQ files using `Trimmomatic`.

```bash
# Install Trimmomatic (via Bioconda)
conda install -c bioconda trimmomatic -y

# Create output directory for trimmed reads
mkdir -p ~/Bulk_RNAseq_Project/trimmed

# Paired-end sample trimming
trimmomatic PE -threads 4 \
~/Bulk_RNAseq_Project/raw_data/SRR12112345_1.fastq \
~/Bulk_RNAseq_Project/raw_data/SRR12112345_2.fastq \
~/Bulk_RNAseq_Project/trimmed/SRR12112345_1.trim.fastq \
~/Bulk_RNAseq_Project/trimmed/SRR12112345_1.unpaired.fastq \
~/Bulk_RNAseq_Project/trimmed/SRR12112345_2.trim.fastq \
~/Bulk_RNAseq_Project/trimmed/SRR12112345_2.unpaired.fastq \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Single-end sample trimming
trimmomatic SE -threads 4 \
~/Bulk_RNAseq_Project/raw_data/SRR12112349.fastq \
~/Bulk_RNAseq_Project/trimmed/SRR12112349.trim.fastq \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
### 4. 🛡️ Quality Control (QC) - After Trimming

This step performs quality assessment on the trimmed FASTQ files using `FastQC` and aggregates reports with `MultiQC`.

```bash
# Navigate to the trimmed reads directory
cd ~/Bulk_RNAseq_Project/trimmed

# Create output directory for QC of trimmed reads
mkdir -p ../qc_trimmed

# Run FastQC on all trimmed FASTQ files using 4 threads
fastqc *.trim.fastq -o ../qc_trimmed/ -t 4

# Aggregate all FastQC reports into a single HTML report
multiqc ../qc_trimmed/ -o ../qc_trimmed/
```
### 5. 🧬 Indexing and Alignment

This step prepares the reference genome for alignment and performs read alignment using `HISAT2`. Both paired-end and single-end reads are aligned, and BAM files are sorted and indexed. Finally, QC statistics are generated with `samtools flagstat`.

```bash
# Create directories for genome index and alignment outputs
mkdir -p ~/Bulk_RNAseq_Project/genome_index/hisat2_index
mkdir -p ~/Bulk_RNAseq_Project/alignment

# Download reference genome and annotation
cd ~/Bulk_RNAseq_Project/genome_index
wget -c https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -c https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Unzip downloaded files
cd ~/Bulk_RNAseq_Project/genome_index && gunzip -v *.gz && ls -lh

# Build HISAT2 index
hisat2-build -p 4 ~/Bulk_RNAseq_Project/genome_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/Bulk_RNAseq_Project/genome_index/hisat2_index

# -------------------------------
# Paired-end reads alignment
for f in ~/Bulk_RNAseq_Project/trimmed/*_1.trim.fastq; do
  base=$(basename "$f" _1.trim.fastq)
  echo "=== $base (paired) ==="
  hisat2 -p 4 -x ~/Bulk_RNAseq_Project/genome_index/hisat2_index \
    -1 ~/Bulk_RNAseq_Project/trimmed/${base}_1.trim.fastq \
    -2 ~/Bulk_RNAseq_Project/trimmed/${base}_2.trim.fastq 2> ~/Bulk_RNAseq_Project/alignment/${base}.hisat2.log \
    | samtools view -@ 4 -bS - \
    | samtools sort -@ 2 -o ~/Bulk_RNAseq_Project/alignment/${base}.sorted.bam
  samtools index ~/Bulk_RNAseq_Project/alignment/${base}.sorted.bam
done

# -------------------------------
# Single-end reads alignment
for f in ~/Bulk_RNAseq_Project/trimmed/*trim.fastq; do
  [[ "$f" == *"_1.trim.fastq" || "$f" == *"_2.trim.fastq" ]] && continue
  base=$(basename "$f" .trim.fastq)
  echo "=== $base (single) ==="
  hisat2 -p 4 -x ~/Bulk_RNAseq_Project/genome_index/hisat2_index -U "$f" 2> ~/Bulk_RNAseq_Project/alignment/${base}.hisat2.log \
    | samtools view -@ 4 -bS - \
    | samtools sort -@ 2 -o ~/Bulk_RNAseq_Project/alignment/${base}.sorted.bam
  samtools index ~/Bulk_RNAseq_Project/alignment/${base}.sorted.bam
done

# -------------------------------
# Alignment QC
for f in ~/Bulk_RNAseq_Project/alignment/*.sorted.bam; do
  base=$(basename "$f" .sorted.bam)
  echo "=== $base ==="
  samtools index "$f"
  samtools flagstat "$f" > ~/Bulk_RNAseq_Project/alignment/${base}.flagstat.txt
  cat ~/Bulk_RNAseq_Project/alignment/${base}.flagstat.txt
done
```

### 6. 📊 Quantification (Raw Counts)

This step generates raw gene counts from the aligned BAM files using `featureCounts`, combines single-end and paired-end counts into a `count matrix`, and creates a `sample info table`.

```bash
# Create output directory for counts
mkdir -p ~/Bulk_RNAseq_Project/counts

# -------------------------------
# Single-end samples
featureCounts -T 4 \
  -a ~/Bulk_RNAseq_Project/genome_index/Homo_sapiens.GRCh38.109.gtf \
  -o ~/Bulk_RNAseq_Project/counts/counts_single.txt \
  ~/Bulk_RNAseq_Project/alignment/single/*.sorted.bam

# Paired-end samples
featureCounts -T 4 -p \
  -a ~/Bulk_RNAseq_Project/genome_index/Homo_sapiens.GRCh38.109.gtf \
  -o ~/Bulk_RNAseq_Project/counts/counts_paired.txt \
  ~/Bulk_RNAseq_Project/alignment/paired/*.sorted.bam

# -------------------------------
# Merge paired and single counts into a single CSV (count matrix)
cd ~/Bulk_RNAseq_Project/counts
Rscript -e "
paired <- read.delim('counts_paired.txt', header=TRUE, row.names=1, check.names=FALSE, comment.char='#', fill=TRUE)
single <- read.table('counts_single.txt', header=TRUE, row.names=1, check.names=FALSE, comment.char='', fill=TRUE)
paired_counts <- paired[,6:ncol(paired)]
single_counts <- single[,6:ncol(single)]
all_counts <- merge(paired_counts, single_counts, by=0, all=TRUE)
rownames(all_counts) <- all_counts\$Row.names
all_counts <- all_counts[,-1]
all_counts[is.na(all_counts)] <- 0
write.csv(all_counts,'counts_matrix.csv')
"

# -------------------------------
# Create sample info file
ls ~/Bulk_RNAseq_Project/alignment/paired/*.sorted.bam ~/Bulk_RNAseq_Project/alignment/single/*.sorted.bam \
  | xargs -n 1 basename \
  | sed 's/.sorted.bam//' \
  | awk 'BEGIN{print "sample,condition"} {print $1",Unknown"}' \
  > ~/Bulk_RNAseq_Project/counts/sample_info.csv
  ```
  
### 7. 📈 Normalization, PCA, DEGs, and Functional Enrichment

This step performs bulk RNA-seq downstream analysis using `DESeq2` and generates normalized counts, PCA plots, differential expression results, top 50 DEGs, heatmaps, volcano plots, and GO/KEGG enrichment plots.

```r
# -----------------------------
# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# -----------------------------
# Create output folders
dir.create("~/Bulk_RNAseq_Project/results", showWarnings = FALSE, recursive = TRUE)
dir.create("~/Bulk_RNAseq_Project/functional", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Step 1: Read counts and sample info
counts <- read.csv('~/Bulk_RNAseq_Project/counts/counts_matrix.csv', row.names=1)
coldata <- read.csv('~/Bulk_RNAseq_Project/counts/sample_info.csv', row.names=1, header=TRUE, sep=',')
stopifnot(all(rownames(coldata) == colnames(counts)))

# -----------------------------
# Step 2: DESeq2 object, normalization, VST, PCA
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts,'~/Bulk_RNAseq_Project/results/normalized_counts.csv')

vsd <- vst(dds, blind=FALSE)
saveRDS(vsd,'~/Bulk_RNAseq_Project/results/vsd_object.rds')

pca <- plotPCA(vsd, intgroup='condition', returnData=TRUE)
pca_plot <- ggplot(pca, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  ggtitle('PCA Plot (VST Transformed Data)') +
  theme_minimal() +
  theme(text=element_text(size=14))
ggsave('~/Bulk_RNAseq_Project/results/PCA_plot.pdf', plot=pca_plot, width=6, height=5)
ggsave('~/Bulk_RNAseq_Project/results/PCA_plot.png', plot=pca_plot, width=6, height=5, dpi=300)

# -----------------------------
# Step 3: DESeq2 differential expression
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$padj[is.na(res_df$padj)] <- 1
write.csv(res_df, '~/Bulk_RNAseq_Project/results/all_results_df.csv')

# -----------------------------
# Step 4: Top50 DEGs and Heatmap
res_sorted <- res_df[order(res_df$padj), ]
top50_DEGs_data <- head(res_sorted, 50)
write.csv(top50_DEGs_data, '~/Bulk_RNAseq_Project/results/top50_DEGs.csv')

top50_genes_names <- rownames(top50_DEGs_data)
if(all(top50_genes_names %in% rownames(assay(vsd)))) {
  pheatmap(assay(vsd)[top50_genes_names, ], scale='row', cluster_cols=TRUE,
           fontsize_row=6, fontsize_col=8,
           filename='~/Bulk_RNAseq_Project/results/heatmap_top50_DEGs.pdf')
  pheatmap(assay(vsd)[top50_genes_names, ], scale='row', cluster_cols=TRUE,
           fontsize_row=6, fontsize_col=8,
           filename='~/Bulk_RNAseq_Project/results/heatmap_top50_DEGs.png')
}

# -----------------------------
# Step 5: Volcano Plot
res_df$color <- ifelse(res_df$padj<0.05 & res_df$log2FoldChange>1, 'red',
                       ifelse(res_df$padj<0.05 & res_df$log2FoldChange< -1, 'blue', 'grey'))
volcano_plot <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=color)) +
  geom_point(size=2) +
  scale_color_identity() +
  theme_minimal() +
  ggtitle('Volcano Plot') +
  xlab('Log2 Fold Change') +
  ylab('-log10(Adjusted P-value)') +
  theme(text=element_text(size=14))
ggsave('~/Bulk_RNAseq_Project/results/volcano_plot.pdf', plot=volcano_plot, width=6, height=5)
ggsave('~/Bulk_RNAseq_Project/results/volcano_plot.png', plot=volcano_plot, width=6, height=5, dpi=300)

# -----------------------------
# Step 6: GO & KEGG enrichment plots
go_results <- read.csv('~/Bulk_RNAseq_Project/functional/GO_BP_top50_DEGs.csv')
kegg_results <- read.csv('~/Bulk_RNAseq_Project/functional/KEGG_top50_DEGs.csv')
genes_clean <- sub('\\..*', '', rownames(top50_DEGs_data))
gene_list <- bitr(genes_clean, fromType='ENSEMBL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID

if(nrow(go_results) > 0) {
  go_object <- new('enrichResult', result=go_results, gene=gene_list, organism='Homo sapiens',
                   keytype='ENTREZID', ontology='BP')
  go_object <- pairwise_termsim(go_object)
  pdf('~/Bulk_RNAseq_Project/functional/GO_BP_Dotplot.pdf', width=8, height=6); dotplot(go_object, showCategory=15); dev.off()
  pdf('~/Bulk_RNAseq_Project/functional/GO_BP_Emapplot.pdf', width=8, height=6); emapplot(go_object, showCategory=15); dev.off()
  pdf('~/Bulk_RNAseq_Project/functional/GO_BP_Cnetplot.pdf', width=8, height=6); cnetplot(go_object, showCategory=10, circular=TRUE, colorEdge=TRUE); dev.off()
} else {
  print('GO WARNING: No significant terms to plot.')
}

if(nrow(kegg_results) > 0) {
  kegg_object <- new('enrichResult', result=kegg_results, gene=gene_list, organism='hsa', keytype='ENTREZID')
  kegg_object <- pairwise_termsim(kegg_object)
  pdf('~/Bulk_RNAseq_Project/functional/KEGG_Dotplot.pdf', width=8, height=6); dotplot(kegg_object, showCategory=15); dev.off()
  pdf('~/Bulk_RNAseq_Project/functional/KEGG_Emapplot.pdf', width=8, height=6); emapplot(kegg_object, showCategory=15); dev.off()
  pdf('~/Bulk_RNAseq_Project/functional/KEGG_Cnetplot.pdf', width=8, height=6); cnetplot(kegg_object, showCategory=10, circular=TRUE, colorEdge=TRUE); dev.off()
} else {
  print('KEGG WARNING: No significant pathways to plot.')
}

print('All results and plots have been generated successfully.')
```