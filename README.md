# QIIME2 Microbiome Analysis Pipeline

A multi-branch 16S rRNA amplicon analysis pipeline, from raw paired-end FASTQ files through denoising, phylogenetics, diversity analysis, and taxonomic classification — with optional 16S rRNA gene copy-number correction (RasperGade16S) and frequency-based contaminant removal (decontam). Downstream statistical analyses and visualizations are performed in R using separately maintained `.Rmd` scripts.

> **Note**: Default parameters are optimized for datasets processed at Karolinska Institutet as part of analyses commissioned by EpiEndo. Adjust rarefaction depth, quality thresholds, and grouping columns to suit your own dataset.

---

## Table of Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [QIIME2](#qiime2)
  - [RasperGade16S Dependencies](#raspergade16s-dependencies-branch-b-only)
  - [R and Required Packages](#r-and-required-packages)
- [Pipeline Architecture](#pipeline-architecture)
  - [Shared Pre-processing (Steps 1–9)](#shared-pre-processing-steps-19)
  - [Branch A — Original (Steps 10–19)](#branch-a--original-steps-1019)
  - [Branch B — RasperGade16S (Steps 20–27)](#branch-b--raspergade16s-steps-2027)
  - [Branch C — decontam (Steps 28–33)](#branch-c--decontam-steps-2833)
  - [Step 34 — Cross-branch QC](#step-34--cross-branch-qc)
- [Inputs](#inputs)
- [Running the Pipeline](#running-the-pipeline)
  - [All Parameters](#all-parameters)
  - [Example Commands](#example-commands)
  - [Checkpoint and Resume System](#checkpoint-and-resume-system)
- [Output Structure](#output-structure)
- [Downstream R Analysis](#downstream-r-analysis)
- [Expected Runtime](#expected-runtime)
- [Citations](#citations)

---

## Overview

The pipeline performs the following analyses, organized into shared steps and three optional branches:

**Shared (all branches):**
1. Build QIIME2 manifest and import paired-end reads
2. Merge paired-end reads (vsearch)
3. Quality filter (Q-score ≥ 30)
4. Deblur denoising (trim length 400 bp)
5. Phylogenetic tree construction (MAFFT + FastTree)

**Branch A — Original:** Standard rarefaction-based diversity analysis, taxonomic classification, PERMANOVA, alpha diversity significance, and Emperor PCoA plots.

**Branch B — RasperGade16S:** 16S rRNA gene copy-number (GCN) correction using [RasperGade16S](https://github.com/vmikk/metagMisc), followed by rarefaction on the corrected integer table and phyloseq object construction.

**Branch C — decontam:** Frequency-based contaminant identification and removal using [decontam](https://benjjneb.github.io/decontam/), requiring a metadata column with DNA concentration values, followed by rarefaction and phyloseq object construction.

**Step 34 (always runs):** Cross-branch QC table comparing sample counts, ASV counts, and library depth distributions across all branches that were executed, plus a unified `phyloseq_all_branches.rds` list for downstream R work.

---

## System Requirements

- **Operating System**: Linux or macOS
- **Memory**: 16 GB RAM minimum; 32 GB recommended for large datasets
- **Storage**: ≥ 50 GB free disk space
- **Processor**: 4+ cores recommended
- **Software**:
  - QIIME2 2024.10 (amplicon distribution)
  - R ≥ 4.0.0
  - Bash shell (POSIX-compliant; uses `set -euo pipefail`)
  - HMMER + Easel (`hmmalign`, `esl-reformat`) — required for Branch B only

---

## Installation

### QIIME2

Install QIIME2 2024.10 using conda:

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

conda env create -n qiime2-amplicon-2024.10 \
    --file qiime2-amplicon-2024.10-py310-linux-conda.yml

conda activate qiime2-amplicon-2024.10

qiime --version
```

For macOS or alternative installation methods, see the [official QIIME2 documentation](https://docs.qiime2.org/2024.10/install/).

### RasperGade16S Dependencies (Branch B only)

Branch B requires HMMER and Easel command-line tools (`hmmalign`, `esl-reformat`) in addition to R packages. Install them into the active conda environment:

```bash
mamba install -c conda-forge -c bioconda easel hmmer

# Verify
which hmmalign && which esl-reformat
```

The workflow validates these dependencies automatically before launching the RasperGade16S R step and will exit with a clear error message if they are missing.

### R and Required Packages

Install R (≥ 4.0.0) from [CRAN](https://cran.r-project.org/), then install the required packages:

```r
# Core tidyverse and plotting
install.packages(c(
  "tidyverse", "dplyr", "ggplot2", "tidyr",
  "ggtext", "ggdist", "ggsignif", "ggrepel",
  "ggforce", "ggside", "ggpubr", "cowplot",
  "vegan", "reshape2", "pals", "scales", "ape"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "biomformat", "decontam"))

# GitHub packages
install.packages("devtools")
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("zhouhj1994/LinDA")

# Branch B: RasperGade16S and sequence parsing
devtools::install_github("wu-lab-uva/RasperGade16S")
```

> **RasperGade16S**: Required package → [wu-lab-uva/RasperGade16S](https://github.com/wu-lab-uva/RasperGade16S)



> If `qiime2R` fails to install from Bioconductor, install it from GitHub as shown above.

---

## Pipeline Architecture

### Shared Pre-processing (Steps 1–9)

All branches run these steps first. Outputs are written to `output/artifacts/`.

| Step | Description |
|------|-------------|
| 1 | Build QIIME2 manifest from sample map and import paired-end reads |
| 2 | Merge paired-end reads with `qiime vsearch merge-pairs` |
| 3 | Summarize demultiplexed data (`.qzv` for quality inspection) |
| 4 | Quality filter with `qiime quality-filter q-score` (min Q = 30) |
| 5 | Tabulate and visualize filter statistics |
| 6 | Deblur denoising (`denoise-16S`, trim length = 400 bp) |
| 7 | Feature-table summary |
| 8 | Deblur statistics visualization |
| 9 | Phylogenetic tree construction: MAFFT alignment → FastTree → rooted tree |

### Branch A — Original (Steps 10–19)

Standard QIIME2 diversity pipeline using the raw Deblur feature table. Steps 17–19 are reserved for downstream R analyses handled separately in R Markdown scripts.

| Step | Description |
|------|-------------|
| 10 | `core-metrics-phylogenetic` — rarefaction at `-d DEPTH` (default: 74,470) with per-sample dropout report |
| 11 | Taxonomic classification with pre-trained sklearn classifier |
| 12 | Taxonomy visualization |
| 13 | Taxa barplot |
| 14 | PERMANOVA for four distance matrices: unweighted UniFrac, weighted UniFrac, Jaccard, Bray–Curtis |
| 15 | Alpha diversity significance tests: Faith's PD, observed features, Shannon, Pielou's evenness |
| 16 | Emperor PCoA plots for all four distance matrices |
| 17–19 | *(Reserved — handled in downstream R scripts)* |

### Branch B — RasperGade16S (Steps 20–27)

Corrects ASV counts for predicted 16S rRNA gene copy number before rarefaction. Taxonomy from Branch A is reused if already computed.

| Step | Description |
|------|-------------|
| 20 | Taxonomic classification (reused from Branch A if `taxonomy.qza` exists) |
| 21 | Export Deblur table, taxonomy, and representative-sequence FASTA for R |
| 22 | **RasperGade16S GCN prediction and correction (R)** — predicts per-ASV copy number from sequences using HMMER/Easel, divides counts by GCN, rounds to integers preserving library depth; outputs both fractional and integer corrected tables |
| 23 | Re-import integer-corrected table into QIIME2; tabulate summary |
| 24 | **Rarefaction depth decision (R)** — compares original vs corrected library sizes, flags samples newly at risk, and saves a diagnostic PDF; *always runs, no checkpoint* |
| 25 | `core-metrics-phylogenetic` on corrected table at rarefaction depth (review step 24 output first; re-run with `-d NEW_DEPTH` if needed) |
| 26 | Export final Branch B artifacts (BIOM, TSV, tree) for R |
| 27 | Build `phyloseq_raspergade.rds` (R) |

> **Before running step 25**, inspect `scripts/step24_rarefaction_decision.pdf`. If samples are newly below the rarefaction depth after GCN correction, re-run Branch B from step 25 with an adjusted `-d` value.

### Branch C — decontam (Steps 28–33)

Identifies and removes contaminant ASVs using the frequency method in decontam. Requires a metadata column with DNA concentration values (`-q`).

| Step | Description |
|------|-------------|
| 28 | Export Deblur table to BIOM/TSV for R |
| 29 | **decontam frequency correction (R)** — calls `isContaminant(..., method="frequency", threshold=0.1)`; outputs `contaminant_scores.tsv` and the filtered table |
| 30 | Re-import decontam-filtered table into QIIME2 |
| 31 | `core-metrics-phylogenetic` on filtered table at rarefaction depth |
| 32 | Export final Branch C artifacts for R |
| 33 | Build `phyloseq_decontam.rds` (R); reuses tree and taxonomy from Branch B |

### Step 34 — Cross-branch QC

Always runs at the end, regardless of which branches were active. Reads all available branch phyloseq `.rds` files, compiles a per-branch QC table (sample count, ASV count, min/median/max depth), and saves:

- `r_objects/cross_branch_qc.tsv` — QC table
- `r_objects/phyloseq_all_branches.rds` — named list of all branch phyloseq objects, ready for downstream comparison

---

## Inputs

### 1. FASTQ Directory (`-f`)

Directory containing gzip-compressed paired-end reads. Files must follow the naming convention:

```
{filename_prefix}_1.fq.gz   # forward read
{filename_prefix}_2.fq.gz   # reverse read
```

### 2. Sample Map (`-s`)

A CSV or TSV file (auto-detected) mapping filename prefixes to QIIME2 sample IDs. Lines beginning with `#`, `sample`, or `sample-id` are treated as headers and skipped.

```csv
filename,sample_id
Bacteria_BX896-001M0001,S001
Bacteria_BX896-001M0002,S002
```

### 3. Metadata File (`-m`)

A QIIME2-compatible TSV file. The first column must be named `sample-id`. All columns referenced by `-g` or `-q` must be present.

```tsv
sample-id	subject-id	visit	group	dna_conc
S001	subject_1	v2	control	12.4
S002	subject_1	v6	treatment	8.7
```

### 4. Output Directory (`-o`)

Created automatically if it does not exist.

### 5. Classifier Directory (`-c`)

Directory containing a pre-trained `.qza` classifier file. The script selects the first `.qza` file found (sorted). Only one classifier per directory is expected.

---

## Running the Pipeline

Activate the QIIME2 conda environment before running the script:

```bash
conda activate qiime2-amplicon-2024.10
```

### All Parameters

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-f` | path | *(required)* | Directory containing FASTQ files |
| `-s` | file | *(required)* | Sample map (CSV or TSV) |
| `-m` | file | *(required)* | QIIME2 metadata file |
| `-o` | dir | *(required)* | Output directory |
| `-c` | dir | *(required)* | Classifier directory |
| `-d` | int | `74470` | Rarefaction depth |
| `-t` | int | `4` | CPU threads |
| `-g` | str | `timepoint` | Metadata column for PERMANOVA grouping |
| `-b` | str | `all` | Branches to run: `original`, `raspergade`, `decontam`, or `all` |
| `-q` | str | *(unset)* | Metadata column with DNA concentrations (required for Branch C) |
| `-S` | int | `1` | Resume from this step number (earlier steps are validated) |
| `-F` | flag | `false` | Force re-run even if output artifacts already exist |
| `-h` | flag | — | Print help and exit |

### Example Commands

**Full run, all branches (decontam requires `-q`):**
```bash
bash workflow.sh \
  -f ./data/fastq \
  -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv \
  -o ./output \
  -c ./classifier \
  -b all \
  -q dna_conc \
  -d 50000 \
  -t 8
```

**Branch A only:**
```bash
bash workflow.sh \
  -f ./data/fastq -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv -o ./output -c ./classifier \
  -b original
```

**Branches A and B only (no decontam):**
```bash
bash workflow.sh \
  -f ./data/fastq -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv -o ./output -c ./classifier \
  -b raspergade
```

**Resume after a failed Deblur step (steps 1–5 artifacts must exist):**
```bash
bash workflow.sh \
  -f ./data/fastq -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv -o ./output -c ./classifier \
  -S 6
```

**Jump directly to RasperGade16S correction (steps 1–21 must already exist):**
```bash
bash workflow.sh \
  -f ./data/fastq -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv -o ./output -c ./classifier \
  -b raspergade -S 22
```

**Redo rarefaction and everything downstream (force re-run from step 10):**
```bash
bash workflow.sh \
  -f ./data/fastq -s ./data/sample_mapping.csv \
  -m ./data/metadata.tsv -o ./output -c ./classifier \
  -F -S 10
```

### Checkpoint and Resume System

Every step is gated by a checkpoint that checks for its expected output artifact. The workflow can be resumed at any step using `-S N`:

- Steps **before** `-S N` are **validated** (the pipeline exits with an error if an expected artifact is missing) but not re-run.
- Steps **at or after** `-S N` are run normally, unless the output already exists and `-F` is not set.
- `-F` (force) bypasses the existence check, causing every step at or after `-S N` to re-run even if its output artifact is present.

This makes it safe to restart after a crash or to re-run a single stage without repeating expensive upstream steps.

---

## Output Structure

```
output/
├── artifacts/                        # Shared QIIME2 artifacts (steps 1–9, Branch A)
│   ├── imported_reads_demux.qza
│   ├── demux-merged.qza / .qzv
│   ├── demux-merged.filtered.qza
│   ├── table.qza / .qzv              # Raw Deblur feature table
│   ├── rep-seqs.qza
│   ├── rooted-tree.qza
│   └── diversities/                  # Branch A core-metrics output
│       ├── weighted_unifrac_distance_matrix.qza
│       ├── weighted_unifrac_pcoa_results.qza
│       ├── weighted_unifrac_emperor.qzv
│       ├── faith_pd_vector.qza
│       ├── evenness_vector.qza
│       └── ...
│
├── taxonomy_classification/
│   ├── taxonomy.qza / .qzv
│   └── taxa-bar-plots.qzv
│
├── statistical_tests/                # Branch A statistical outputs
│   ├── permanova_unweighted_unifrac_distance_matrix.qzv
│   ├── permanova_weighted_unifrac_distance_matrix.qzv
│   ├── permanova_jaccard_distance_matrix.qzv
│   ├── permanova_bray_curtis_distance_matrix.qzv
│   ├── alpha_significance_faith_pd.qzv
│   ├── alpha_significance_observed_features.qzv
│   ├── alpha_significance_shannon.qzv
│   └── alpha_significance_evenness.qzv
│
├── branch_raspergade/                # Branch B QIIME2 artifacts
│   ├── table-corrected.qza / -summary.qzv
│   └── diversities/
│
├── branch_decontam/                  # Branch C QIIME2 artifacts
│   ├── table-decontam.qza / -summary.qzv
│   └── diversities/
│
├── r_objects/                        # Exported tables and phyloseq objects
│   ├── original/
│   │   ├── feature-table.biom / .tsv
│   │   └── phyloseq_original.rds
│   ├── raspergade/
│   │   ├── feature-table.biom / .tsv    # Raw Deblur table (pre-correction)
│   │   ├── corrected-table.tsv          # Integer GCN-corrected table
│   │   ├── corrected-table-fractional.tsv
│   │   ├── raspergade16s-gcn.tsv        # Per-ASV predicted GCN + taxonomy
│   │   ├── asvs-missing-raspergade-gcn.tsv  # ASVs without GCN prediction (if any)
│   │   ├── taxonomy/taxonomy.tsv
│   │   ├── rep-seqs/dna-sequences.fasta
│   │   ├── final/feature-table.biom / .tsv
│   │   └── phyloseq_raspergade.rds
│   ├── decontam/
│   │   ├── contaminant_scores.tsv       # ASV-level decontam p-values
│   │   ├── decontam-table.tsv / .rds    # Filtered table
│   │   ├── final/feature-table.biom / .tsv
│   │   └── phyloseq_decontam.rds
│   ├── cross_branch_qc.tsv             # Step 34: per-branch QC summary
│   └── phyloseq_all_branches.rds       # Step 34: named list of all phyloseq objects
│
└── scripts/                           # Auto-generated R scripts (with logs)
    ├── step22_raspergade.R / .log
    ├── step24_rarefaction_decision.R / .log
    ├── step24_rarefaction_decision.pdf  # Library size comparison plot
    ├── step27_phyloseq_raspergade.R / .log
    ├── step29_decontam.R / .log
    ├── step33_phyloseq_decontam.R / .log
    └── step34_cross_branch_qc.R / .log
```

### Viewing QIIME2 Visualizations

All `.qzv` files can be viewed interactively at [https://view.qiime2.org/](https://view.qiime2.org/) — no local QIIME2 installation required for viewing.

---

## Downstream R Analysis

After `workflow.sh` completes, the `phyloseq_all_branches.rds` object in `r_objects/` is the primary input to the R Markdown analysis scripts located in `R/`:

| Script | Purpose |
|--------|---------|
| `adonis_pairwise.Rmd` | Pairwise PERMANOVA (beta diversity) with subject-level stratification |
| `alpha_div_test.Rmd` | Alpha diversity analysis using linear mixed-effects models (richness and evenness) |
| `linDA.Rmd` | Differential abundance analysis across taxonomic levels using LinDA |
| `phenotypic_analysis.Rmd` | Analysis of quantitative culture data (log CFU) with mixed-effects models |
| `phylo_tree.Rmd` | Construction and annotation of ASV phylogenetic tree for visualization in iTOL |
| `plot_taxa_diff.Rmd` | Visualization of selected taxa (CLR-transformed abundance, boxplots, longitudinal trends) |
| `utility.Rmd` | Helper functions for data processing, transformation, and plotting |

---

## Expected Runtime

The EpiEndo dataset (Karolinska Institutet) required approximately **3 days** on a Linux system (Intel Xeon Gold 6430, 64 GB RAM, 1 TB storage). The majority of runtime is consumed by the Deblur denoising step (step 6). The RasperGade16S GCN prediction step (step 22) can also be time-consuming depending on the number of ASVs and available CPU cores.

---

## Citations

If you use this pipeline, please cite:

**QIIME2:**
Bolyen E, et al. (2019) Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology* 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9

**Deblur:**
Amir A, et al. (2017) Deblur rapidly resolves single-nucleotide community sequence patterns. *mSystems* 2: e00191-16. https://doi.org/10.1128/mSystems.00191-16

**RasperGade16S:**
Consult the [RasperGade16S GitHub repository](https://github.com/vmikk/metagMisc) for the appropriate citation.

**decontam:**
Davis NM, et al. (2018) Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6: 226. https://doi.org/10.1186/s40168-018-0605-2

**phyloseq:**
McMurdie PJ and Holmes S (2013) phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE* 8: e61217. https://doi.org/10.1371/journal.pone.0061217

**LinDA:**
Zhou H, et al. (2022) LinDA: linear models for differential abundance analysis of microbiome compositional data. *Genome Biology* 23: 95. https://doi.org/10.1186/s13059-022-02655-5

**ggplot2:**
Wickham H (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

---

## Getting Help

- **QIIME2 Documentation**: https://docs.qiime2.org/
- **QIIME2 Forum**: https://forum.qiime2.org/
- **decontam**: https://benjjneb.github.io/decontam/
- **LinDA**: https://github.com/zhouhj1994/LinDA

---

## License

This pipeline is provided for research use. Please ensure you comply with the licenses of all included software packages.
