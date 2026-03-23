# QIIME2 Microbiome Analysis Pipeline

A pipeline for 16S rRNA microbiome data analysis, from raw FASTQ files to statistical analysis and visualization. This pipeline is optimized for datasets processed at Karolinska Institutet as part of analyses commissioned by EpiEndo.

---

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [QIIME2 Installation](#qiime2-installation)
  - [R and Required Packages](#r-and-required-packages)
- [Pipeline Structure](#pipeline-structure)
- [Running the Pipeline](#running-the-pipeline)
  - [1. QIIME2 Workflow (Bash)](#1-qiime2-workflow-bash)
  - [2. Visualization and Statistics (R)](#2-visualization-and-statistics-r)
- [Folder Structure](#folder-structure)
- [Expected Runtime](#expected-runtime)
- [Output Files](#output-files)
- [Citation](#citation)

---

## Overview

This pipeline performs:

1. **Data Import and Quality Control**: Import paired-end FASTQ files, merge reads, and perform quality filtering
2. **Feature Table Construction**: Denoise sequences using Deblur
3. **Phylogenetic Analysis**: Construct phylogenetic trees using MAFFT and FastTree
4. **Diversity Analysis**: Calculate alpha and beta diversity metrics
5. **Taxonomic Classification**: Assign taxonomy using pre-trained classifiers
6. **Statistical Testing**: PERMANOVA, alpha diversity tests, and LinDA differential abundance analysis
7. **Visualization**: Generate publication-ready plots in R

---

## System Requirements

- **Operating System**: Linux or macOS
- **Memory**: Minimum 16 GB RAM (32 GB recommended for large datasets)
- **Storage**: At least 50 GB free disk space
- **Processor**: Multi-core processor recommended (4+ cores)
- **Software**:
  - QIIME2 2024.10
  - R ≥ 4.0.0
  - Bash shell

---

## Installation

### QIIME2 Installation

Install QIIME2 2024.10 using conda:

```bash
# Download the QIIME2 conda environment file
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

# Create the conda environment
conda env create -n qiime2-amplicon-2024.10 --file qiime2-amplicon-2024.10-py310-linux-conda.yml

# Activate the environment
conda activate qiime2-amplicon-2024.10

# Verify installation
qiime --version
```

For macOS or other installation methods, see the [official QIIME2 installation guide](https://docs.qiime2.org/2024.10/install/).

### R and Required Packages

Install R (≥ 4.0.0) from [CRAN](https://cran.r-project.org/), then install required packages:

```r
# Install CRAN packages
install.packages(c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "tidyverse",
  "ggtext",
  "ggdist",
  "ggsignif",
  "ggrepel",
  "ggforce",
  "ggside",
  "ggpubr",
  "cowplot",
  "vegan",
  "reshape2"
))

# Install color palette packages
install.packages(c("pals"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("qiime2R"))

# Verify installations
library(qiime2R)
library(ggplot2)
library(dplyr)
```

**Note**: If you encounter issues installing `qiime2R`, you may need to install it from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install qiime2R from GitHub
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("zhouhj1994/LinDA")
```

---

## Pipeline Structure

```
q2-pipeline/
├── workflow.sh              # Main QIIME2 analysis script
├── R/
├── adonis_pairwise.Rmd      # Pairwise PERMANOVA analysis (beta-diversity) with subject-level stratification
├── alpha_div_test.Rmd       # Alpha diversity analysis using mixed-effects models (richness and evenness)
├── linDA.Rmd                # Differential abundance analysis across taxonomic levels using LinDA
├── phenotypic_analysis.Rmd  # Analysis of quantitative culture data (logCFU) using mixed-effects models
├── phylo_tree.Rmd           # Construction and annotation of ASVs phylogenetic tree to be visualized in iTOL
├── plot_taxa_diff.Rmd       # Visualization of selected taxa (CLR abundance, boxplots, longitudinal trends)
└── utility.Rmd              # Helper functions for data processing, transformation, and plotting
├── README.md               # This file
└── classifier/            # Taxonomic classifier directory
   └── classifier.qza     # Pre-trained classifier
```
---

## Running the Pipeline

### 1. QIIME2 Workflow (Bash)

The main workflow script (`workflow.sh`) performs all QIIME2 analyses.

#### Required Inputs:

1. **FASTQ directory** (`-f`): Directory containing paired-end FASTQ files
   - Files should be named: `{filename}_1.fq.gz` and `{filename}_2.fq.gz`

2. **Sample mapping file** (`-s`): CSV file mapping FASTQ filenames to sample IDs
   ```csv
   filename,sample_id
   Bacteria_BX896-001M0001,S001
   Bacteria_BX896-001M0002,S002
   ```

3. **Metadata file** (`-m`): TSV file with sample metadata (QIIME2 format)
   ```tsv
   sample-id	subject-id	visit	group
   S001	subject_1	v2	control
   S002	subject_1	v6	treatment
   ```

4. **Output directory** (`-o`): Directory for output files

5. **Classifier file** (`-c`): Directory containing the `.qza` classifier file

6. **Grouping variable** (`-g`): Metadata column name for statistical tests (e.g., "visit", "time_point")

#### Optional Parameters:

- **Sampling depth** (`-d`): Rarefaction depth (default: 74470)
- **Threads** (`-t`): Number of CPU threads (default: 4)

#### Example Command:

```bash
# Activate QIIME2 environment
conda activate qiime2-amplicon-2024.10

# Run the pipeline
bash workflow.sh \
  -f ./data/fastq \
  -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv \
  -o ./output \
  -c ./classifier \
  -g time_point \
  -d 50000 \
  -t 8
```

#### Help:

```bash
bash workflow.sh -h
```


### Output Structure (after running workflow.sh):

```
output/
├── artifacts/
│   ├── imported_reads_demux.qza
│   ├── demux-merged.qza
│   ├── table.qza
│   ├── rep-seqs.qza
│   ├── rooted-tree.qza
│   ├── metadata.tsv
│   └── diversities/
│       ├── weighted_unifrac_pcoa_results.qza
│       ├── faith_pd_vector.qza
│       ├── evenness_vector.qza
│       └── distance_matrix/
│           └── weighted_unifrac_distance_matrix.tsv
├── taxonomy_classification/
│   ├── taxonomy.qza
│   ├── taxonomy.qzv
│   └── taxa-bar-plots.qzv
├── statistical_tests/
│   ├── permanova_*.qzv
│   ├── alpha_significance_*.qzv
│   ├── ancombc_*.qza
│   └── summary_report.txt
└── R_export/
    ├── alpha_diversity.pdf
    ├── pcoa.pdf
    ├── read_stats_boxplot-heatmap.pdf
    ├── ancombc_ASVs_all_ci95.csv
    ├── sig_Genus_lf1_v7-v8.pdf
    ├── etc.
```

---

## Expected Runtime
 The EpiEndo dataset required approximately 3 days on a Linux system (Intel Xeon Gold 6430, 64GB RAM, 1TB storage), with the majority of time spent in the Deblur denoising step.

---

## Output Files

### QIIME2 Artifacts (.qza)

- **table.qza**: Feature table (ASV counts)
- **rep-seqs.qza**: Representative sequences
- **rooted-tree.qza**: Phylogenetic tree
- **taxonomy.qza**: Taxonomic assignments

### QIIME2 Visualizations (.qzv)

View these files at [https://view.qiime2.org/](https://view.qiime2.org/)

- **demux-merged.qzv**: Read quality and count statistics
- **taxa-bar-plots.qzv**: Interactive taxonomic composition plots
- **permanova_*.qzv**: Beta diversity statistical test results
- **alpha_significance_*.qzv**: Alpha diversity test results

### R Outputs (PDF and CSV)

- **alpha_diversity.pdf**: Alpha diversity metrics with statistical significance
- **pcoa.pdf**: Principal Coordinates Analysis visualization
- **read_stats_boxplot-heatmap.pdf**: Read count distribution
- **ancombc_*_all_ci95.csv**: Differential abundance results with confidence intervals
- **normalized_otu_counts.csv**: Bias-corrected abundance values
- **sig_*.pdf**: Plots of significantly differentially abundant taxa

### Getting Help:

- **QIIME2 Documentation**: [https://docs.qiime2.org/](https://docs.qiime2.org/)
- **QIIME2 Forum**: [https://forum.qiime2.org/](https://forum.qiime2.org/)
- **LinDA**: [https://github.com/zhouhj1994/LinDA](https://github.com/zhouhj1994/LinDA)

---

## Citation

If you use this pipeline, please cite:

**QIIME2:**
- Bolyen E, et al. (2019) Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology* 37: 852-857. doi: 10.1038/s41587-019-0209-9

**Deblur:**
- Amir A, et al. (2017) Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems* 2: e00191-16. doi: 10.1128/mSystems.00191-16

**LinDA**
- Zhou H., et al. (2022) LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biol 23, 95 . doi.org/10.1186/s13059-022-02655-5

**Visualization Packages:**
- Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

---

## License

This pipeline is provided for research use. Please ensure you comply with the licenses of all included software packages.

---
