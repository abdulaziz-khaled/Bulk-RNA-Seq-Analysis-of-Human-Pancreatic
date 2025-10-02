# 🧬 Bulk RNA-seq Analysis Reveals Gene Expression changes in Human Pancreatic β-cells Under Metabolic Stress.
This project performs bulk RNA-Seq data analysis on human pancreatic β-cells under metabolic stress using dataset GSE159984. The pipeline includes data retrieval from GEO using SRA Toolkit, quality control with FastQC and MultiQC, read trimming with Trimmomatic, read alignment to the human reference genome (hg38) using STAR/HISAT2, gene quantification with featureCounts, normalization and differential expression analysis with DESeq2, and visualization (PCA, heatmap, and volcano plot). The analysis compares 3 control samples (paired-end) against 3 stressed samples (single-end). All tools used are open-source and run on Ubuntu via WSL (Windows Subsystem for Linux).

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