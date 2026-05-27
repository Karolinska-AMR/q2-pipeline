#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# render_comparison.R
# Run this script to set up the directory structure and render the Rmd
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. Install missing packages ────────────────────────────────────────────────
pkgs <- c("tidyverse", "vegan", "ggplot2", "patchwork", "cowplot",
          "UpSetR", "scales", "ggrepel", "pheatmap", "RColorBrewer",
          "knitr", "kableExtra", "gt", "gtsummary", "ape", "rmarkdown")

for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
        message("Installing: ", p)
        install.packages(p, quiet = TRUE)
    }
}

if (!requireNamespace("qiime2R", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("jbisanz/qiime2R")
}

# ── 2. Directory structure ─────────────────────────────────────────────────────
# Set this to your project base directory
PROJECT_DIR <- "~/projects/microbiome"

# Source files from workflow output
WORKFLOW_DIR <- file.path(PROJECT_DIR, "workflow_results")

# Analysis output directory
ANALYSIS_DIR <- file.path(WORKFLOW_DIR, "branch_comparison")
dir.create(ANALYSIS_DIR, showWarnings = FALSE, recursive = TRUE)

setwd(ANALYSIS_DIR)

# Create required input directories
dir.create("input/original",              showWarnings = FALSE, recursive = TRUE)
dir.create("input/raspergade",            showWarnings = FALSE, recursive = TRUE)
dir.create("input/stats/original",        showWarnings = FALSE, recursive = TRUE)
dir.create("input/stats/raspergade",      showWarnings = FALSE, recursive = TRUE)

# ── 3. Copy input files ────────────────────────────────────────────────────────
message("Copying input files...")

copy_if_exists <- function(from, to) {
    if (file.exists(from)) {
        file.copy(from, to, overwrite = TRUE)
        message("  ✓ ", basename(to))
    } else {
        message("  ✗ Missing: ", from)
    }
}

# Original branch
copy_if_exists(
    file.path(WORKFLOW_DIR, "artifacts/table.qza"),
    "input/original/table.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "taxonomy_classification/taxonomy.qza"),
    "input/original/taxonomy.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "artifacts/diversities/faith_pd_vector.qza"),
    "input/original/faith_pd_vector.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "artifacts/diversities/evenness_vector.qza"),
    "input/original/evenness_vector.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "artifacts/diversities/weighted_unifrac_distance_matrix.qza"),
    "input/original/weighted_unifrac_distance_matrix.qza"
)
copy_if_exists(
  file.path(WORKFLOW_DIR, "artifacts/diversities/weighted_unifrac_pcoa_results.qza"),
  "input/original/weighted_unifrac_pcoa_results.qza"
)
# RasperGade branch
copy_if_exists(
    file.path(WORKFLOW_DIR, "branch_raspergade/table-corrected.qza"),
    "input/raspergade/table-corrected.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "branch_raspergade/diversities/faith_pd_vector.qza"),
    "input/raspergade/faith_pd_vector.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "branch_raspergade/diversities/evenness_vector.qza"),
    "input/raspergade/evenness_vector.qza"
)
copy_if_exists(
    file.path(WORKFLOW_DIR, "branch_raspergade/diversities/weighted_unifrac_distance_matrix.qza"),
    "input/raspergade/weighted_unifrac_distance_matrix.qza"
)
copy_if_exists(
  file.path(WORKFLOW_DIR, "branch_raspergade/diversities/weighted_unifrac_pcoa_results.qza"),
  "input/raspergade/weighted_unifrac_pcoa_results.qza"
)
# Metadata
copy_if_exists(
    file.path(PROJECT_DIR, "input_data/meta_data.tsv"),
    "input/meta_data.tsv"
)

# ── 4. Copy CSV statistics ─────────────────────────────────────────────────────
message("\nCopying statistical outputs...")

# Adjust these paths to where your R downstream scripts saved their CSVs
STATS_ORIG <- file.path(WORKFLOW_DIR, "export_orig")
STATS_RASP <- file.path(WORKFLOW_DIR, "branch_raspergade/export")

stat_files <- c(
    "faith_pd_lmm_anova.csv",
    "faith_pd_lmm_pairwise_emmeans.csv",
    "pielou_evenness_lmm_anova.csv",
    "pielou_evenness_lmm_pairwise_emmeans.csv",
    "weighted_unifrac_pairwise_adonis.csv",
    "linDA_ASVs_Sig_LFC1.csv",
    "linDA_Genus_Sig_LFC1.csv",
    "linDA_Family_Sig_LFC1.csv",
    "linDA_Order_Sig_LFC1.csv",
    "linDA_Class_Sig_LFC1.csv",
    "linDA_Phylum_Sig_LFC1.csv"
)

for (f in stat_files) {
    copy_if_exists(file.path(STATS_ORIG, f), file.path("input/stats/original", f))
    copy_if_exists(file.path(STATS_RASP, f), file.path("input/stats/raspergade", f))
}

message("\n✓ Report saved: ", file.path(ANALYSIS_DIR, "EP395_branch_comparison.html"))
message("✓ Figures:      ", file.path(ANALYSIS_DIR, "figures/"))
message("✓ Tables:       ", file.path(ANALYSIS_DIR, "tables/"))
