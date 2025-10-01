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
