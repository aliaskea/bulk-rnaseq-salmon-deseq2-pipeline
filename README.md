# PTMA CD8 RNA-seq Analysis (Salmon + DESeq2)
Reproducible bulk RNA-seq workflow for differential expression analysis of PTMA overexpression vs vector control in mouse OT1 CD8+ T cells.

PTMA (Prothymosin Alpha) has been implicated in regulating T cell persistence and mitochondrial integrity.
This project evaluates whether PTMA overexpression induces global transcriptomic changes in OT1 CD8+ T cells using bulk RNA-seq.

We test:

H0: No gene expression difference between vector and PTMA OE  
H1: PTMA OE significantly alters transcriptional programs

## Dataset
- GEO accession: `GSE307982`
- Organism: `Mus musculus`
- Design: 6 samples total (3 `vector`, 3 `ptma_oe`), paired-end RNA-seq

## Workflow Overview
1. FastQC on raw reads
2. Build Salmon transcriptome index
3. Quantify transcript abundance with Salmon
4. Build transcript-to-gene mapping from GTF
5. Differential expression with DESeq2
6. Aggregate QC with MultiQC

## Project Structure
- `config/` - sample sheet, contrasts, and analysis paths
- `workflow/` - Snakemake entrypoint and rules
- `scripts/` - R scripts for tx2gene and DESeq2
- `data/` - raw inputs
- `refs/` - reference FASTA/GTF
- `results/` - pipeline outputs

## Setup
```bash
conda env create -f environment.yml
conda activate ptma-cd8-rnaseq
```

## Configure Inputs
Edit:
- `config/samples.tsv`
- `config/contrasts.tsv`
- `config/config.yaml`

Expected key inputs:
- FASTQ files in `data/fastq/`
- Transcript FASTA in `refs/Mus_musculus.GRCm39.cdna.all.fa.gz`
- GTF in `refs/Mus_musculus.GRCm39.110.gtf.gz`

## Run
```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

## Main Outputs
- Salmon quantification: `results/salmon/<sample>/quant.sf`
- FastQC and MultiQC: `results/qc/`
- tx2gene map: `results/reference/tx2gene.tsv`
- DESeq2 results: `results/deseq2/deseq2_results.tsv`

## Key Results
- Differential expression analysis identified **2 significantly upregulated genes** (FDR<0.05) in PTMA-overexpressing CD8+ T cells.
- No genes were significantly downregulated at the same threshold.
- PCA did not show clear global separation between conditions, suggesting limited transcriptome-wide remodeling.
- The two significant genes (B4galt5 and Gcc2) are involved in protein glycosylation and Golgi-associated trafficking, indicating a potentially targeted cellular effect rather than broad transcriptional reprogramming.

## Note
This repository is an independent implementation inspired by public Salmon+DESeq2 Snakemake workflows (https://github.com/niekwit/rna-seq-salmon-deseq2/tree/main), adapted for this specific dataset and analysis goals.
