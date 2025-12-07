# ANCOM-BC Analysis Script
# Differential abundance analysis of microbiome data using ANCOM-BC

# Load required libraries ----
library(dplyr)
library(qiime2R)
library(ANCOMBC)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggtext)
library(pals)


args <- commandArgs(trailingOnly = TRUE)
results_dir <<- args[1]
metadata_tsv <<- args[2]

print(results_dir)
print(metadata_tsv)

#Set working directory
setwd(results_dir)
export_dir <- "./R_export"
# Create export directory if it doesn't exist
if (!dir.exists(export_dir)) {
  dir.create(export_dir, recursive = TRUE)
}

# Global function to load and process metadata
VISIT_LEVELS <<- c("Baseline 1", "Baseline 2", "1–week EP395", "2–week EP395")

load_metadata <- function() {
  cat("Loading metadata...\n")
  
  metadata <- read.table(metadata_tsv, header = TRUE,
                         stringsAsFactors = FALSE,sep = '\t') %>% 
    as.data.frame() %>% mutate(timepoint = gsub("-", "–", timepoint, fixed = TRUE)) %>%
    mutate(timepoint = factor(timepoint, levels = VISIT_LEVELS))
  
  cat("Metadata loaded successfully\n")
  cat(sprintf("  - Samples: %d\n", nrow(metadata)))
  
  return(metadata)
}

# SECTION 1: ASV-level ANCOM-BC Analysis ----

## Load and prepare data ----
load_and_prepare_data <- function() {
  # Load metadata
  metadata <- load_metadata()
  row.names(metadata)<- metadata$`sample.id`
  print(head(metadata))
  # Load feature table
  feature_table <- read_qza("./artifacts/table.qza")$data
  
  # Load and process taxonomy
  taxonomy_annot <- read_qza('./taxonomy_classification/taxonomy.qza')$data
  taxonomy_annot <- taxonomy_annot %>% 
    mutate(taxa = str_remove_all(Taxon, "D_\\d__")) %>%
    separate_wider_delim(taxa, "; ", 
                        names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                        too_few = "align_start") %>%
    mutate_at(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             list(~ substr(., 4, nchar(.))))
  
  row.names(taxonomy_annot) <- taxonomy_annot$Feature.ID
  feature_table <- as.matrix(feature_table)
  
  return(list(metadata = metadata, feature_table = feature_table, taxonomy_annot = taxonomy_annot))
}

#clean the output columns of ANCOMBC
clean_column_names <- function(colnames_vector) {
  
  # Replace long timepoint strings with readable forms
  colnames_vector <- str_replace_all(colnames_vector, "timepointBaseline 2", "baseline_2")
  colnames_vector <- str_replace_all(colnames_vector, "timepoint1–week EP395", "week1_ep395")
  colnames_vector <- str_replace_all(colnames_vector, "timepoint2–week EP395", "week2_ep395")
  
  # Replace spaces and special characters with underscores
  colnames_vector <- str_replace_all(colnames_vector, "[ –]", "_")  # en-dash or space to underscore
  
  return(colnames_vector)
}

## Run ANCOM-BC analysis ----
run_ancombc_asv <- function(feature_table, metadata) {
  cat("Running ANCOM-BC analysis at ASV level...\n")
  
  results <- ancombc(
    data = feature_table, 
    meta_data = metadata,
    formula = "timepoint",
    p_adj_method = "bonferroni", 
    prv_cut = 0.05, 
    lib_cut = 0,
    group = "timepoint", 
    struc_zero = TRUE, 
    neg_lb = TRUE, 
    tol = 1e-5,
    max_iter = 100, 
    conserve = TRUE, 
    alpha = 0.05, 
    global = TRUE,
    n_cl = 5, 
    verbose = TRUE
  )
  
  return(results)
}

## Process and export normalized counts ----
export_normalized_counts <- function(results, suffix = "") {
  # Extract sampling fractions
  samp_frac <- results$samp_frac
  samp_frac[is.na(samp_frac)] <- 0
  
  # Calculate bias-corrected abundances
  log_obs_abn <- log(results$feature_table + 1)
  log_corr_abn <- t(t(log_obs_abn) - samp_frac)
  corr_abn <- exp(log_corr_abn)
  
  # Export files
  
  write.csv(corr_abn, file.path(export_dir,paste0('normalized_otu_counts', suffix, '.csv')))
  write.csv(log_corr_abn, file.path(export_dir,paste0('log_normalized_otu_counts', suffix, '.csv')))
  
  cat("Normalized counts exported successfully.\n")
  return(log_corr_abn)
}

## Process ANCOM-BC results and calculate confidence intervals ----
process_ancombc_results <- function(results, taxonomy_annot, rank_name = "ASV") {
  res <- results$res
  # Extract individual result tables
  tab_lfc <- res$lfc %>% rename_with(~ paste0("lfc_", .x), -taxon)
  tab_q <- res$q %>% rename_with(~ paste0("q_", .x), -taxon)
  tab_diff <- res$diff_abn %>% rename_with(~ paste0("diff_", .x), -taxon)
  tab_w <- res$W %>% rename_with(~ paste0("w_", .x), -taxon)
  tab_se <- res$se %>% rename_with(~ paste0("se_", .x), -taxon)
  
  # Calculate 95% confidence intervals
  z <- 1.96
  
  stat_df <- tab_diff %>% 
    merge(tab_q, by = 'taxon') %>% 
    merge(tab_lfc, by = 'taxon') %>% 
    merge(tab_w, by = 'taxon') %>% 
    merge(tab_se, by = 'taxon') %>%
    mutate(across(starts_with("w_"), 
                  .names = "lci_{.col}", 
                  ~ . - z * get(paste0("se_", sub("w_", "", cur_column()))))) %>%
    mutate(across(starts_with("w_"),  
                  .names = "lcu_{.col}",
                  ~ . + z * get(paste0("se_", sub("w_", "", cur_column())))))
 
  # Add taxonomy information for ASV-level analysis
  if (rank_name == "ASV") {
    stat_df <- stat_df %>% 
      merge(taxonomy_annot, by.x = 'taxon', by.y = 'Feature.ID') %>% 
      rename('feature_id' = 'taxon', 'taxon' = 'Taxon')

        # Define column order for ASV analysis
    columns <- c("feature_id", 'taxon', "w_(Intercept)", "se_(Intercept)", "lci_w_(Intercept)", 
                "lcu_w_(Intercept)", "lfc_(Intercept)", "q_(Intercept)", "diff_(Intercept)", 
                "w_timepointBaseline 2", "se_timepointBaseline 2", "lci_w_timepointBaseline 2", "lcu_w_timepointBaseline 2", "lfc_timepointBaseline 2", 
                "q_timepointBaseline 2", "diff_timepointBaseline 2", "w_timepoint1–week EP395", "se_timepoint1–week EP395", "lci_w_timepoint1–week EP395", 
                "lcu_w_timepoint1–week EP395", "lfc_timepoint1–week EP395", "q_timepoint1–week EP395", "diff_timepoint1–week EP395", "w_timepoint2–week EP395", 
                "se_timepoint2–week EP395", "lci_w_timepoint2–week EP395", "lcu_w_timepoint2–week EP395", "lfc_timepoint2–week EP395", "q_timepoint2–week EP395", "diff_timepoint2–week EP395")
    
    stat_df <- stat_df[columns]
    filename <- file.path(export_dir,'ancombc_ASVs_all_ci95.csv')
  } else {
    # For aggregated analysis
    stat_df <- stat_df %>% rename(!!rank_name := 'taxon')
    
    columns <- c(rank_name, "w_(Intercept)", "se_(Intercept)", "lci_w_(Intercept)", 
                "lcu_w_(Intercept)", "lfc_(Intercept)", "q_(Intercept)", "diff_(Intercept)", 
                "w_timepointBaseline 2", "se_timepointBaseline 2", "lci_w_timepointBaseline 2", "lcu_w_timepointBaseline 2", "lfc_timepointBaseline 2", 
                "q_timepointBaseline 2", "diff_timepointBaseline 2", "w_timepoint1–week EP395", "se_timepoint1–week EP395", "lci_w_timepoint1–week EP395", 
                "lcu_w_timepoint1–week EP395", "lfc_timepoint1–week EP395", "q_timepoint1–week EP395", "diff_timepoint1–week EP395", "w_timepoint2–week EP395", 
                "se_timepoint2–week EP395", "lci_w_timepoint2–week EP395", "lcu_w_timepoint2–week EP395", "lfc_timepoint2–week EP395", "q_timepoint2–week EP395", "diff_timepoint2–week EP395")
    
    stat_df <- stat_df[columns]
    filename <- file.path(export_dir,paste0('ancombc_', rank_name, '_all_ci95.csv'))
  }
  # clean columns
  colnames(stat_df) <- clean_column_names(colnames(stat_df))
  
  write_csv(stat_df, filename)
  cat("Statistical results exported to:", filename, "\n")
  
  return(stat_df)
}


# SECTION 2: Aggregated Analysis (Genus-level) ----

## Aggregate features by taxonomic rank ----
aggregate_features <- function(feature_table, taxonomy_annot, rank = 'Genus') {
  cat("Aggregating features by", rank, "level...\n")
  
  taxa_meta <- taxonomy_annot %>% select(c(Feature.ID, all_of(rank)))
  
  feature_table_agg <- as.data.frame(feature_table) %>% 
    rownames_to_column(var = "Feature.ID") %>%
    merge(taxa_meta, by = "Feature.ID", all.x = TRUE) %>% 
    select(-'Feature.ID') %>%
    mutate(across(!starts_with(rank), as.numeric)) %>%
    group_by(across(all_of(rank))) %>% 
    summarize(across(where(is.numeric), sum)) %>% 
    ungroup() %>% 
    as.data.frame()
  
  rank_val <- feature_table_agg[, rank]
  row.names(feature_table_agg) <- replace(rank_val, is.na(rank_val), 'unknown')
  feature_table_agg <- feature_table_agg[, -1]
  feature_table_agg <- as.matrix(feature_table_agg)
  
  return(feature_table_agg)
}

## Run ANCOM-BC on aggregated data ----
run_ancombc_aggregated <- function(feature_table_agg, metadata) {
  cat("Running ANCOM-BC analysis on aggregated data...\n")
  
  results <- ancombc(
    data = feature_table_agg, 
    meta_data = metadata,
    formula = "timepoint",
    p_adj_method = "holm", 
    prv_cut = 0.10, 
    lib_cut = 1000,
    group = "timepoint", 
    struc_zero = TRUE, 
    neg_lb = TRUE, 
    tol = 1e-5,
    max_iter = 100, 
    conserve = TRUE, 
    alpha = 0.05, 
    global = TRUE,
    n_cl = 5, 
    verbose = TRUE
  )
  
  return(results)
}

## Create bar plot for significant aggregated taxa ----
plot_significant_aggregated <- function(stat_df, rank = 'Genus', lfc_threshold = 1) {
  df_filter <- stat_df %>% 
    filter((diff_visitv7 == TRUE & abs(lfc_visitv7) > lfc_threshold) | 
           (diff_visitv8 == TRUE & abs(lfc_visitv8) > lfc_threshold)) %>% 
    select(all_of(rank), 'lfc_visitv8', 'lfc_visitv7', "se_visitv8", "se_visitv7") %>% 
    rename('v8' = 'lfc_visitv8', 'v7' = 'lfc_visitv7', 
           'se8' = "se_visitv8", 'se7' = "se_visitv7", 
           'taxa_lvl' := !!rank)
  
  if (nrow(df_filter) == 0) {
    cat("No significant", rank, "found with LFC threshold of", lfc_threshold, "\n")
    return(NULL)
  }
  
  # Prepare melted data
  melt_df <- bind_cols(
    df_filter %>% select(taxa_lvl, v7, v8) %>% 
      melt() %>% rename('visit' = 'variable') %>% select(taxa_lvl, visit, value),
    df_filter %>% select(taxa_lvl, se7, se8) %>% 
      melt() %>% rename('visit_se' = 'value') %>% select(visit_se)
  )
  
  # Order taxa by effect size
  melt_df$taxa_lvl <- factor(melt_df$taxa_lvl, 
                            levels = unique(melt_df[order(melt_df$value), "taxa_lvl"]))
  
  # Create plot
  p <- ggplot(data = melt_df, 
              aes(x = taxa_lvl, y = value, fill = visit, color = visit, group = visit)) + 
    geom_bar(stat = "identity", width = 0.7, color = "#5E686D", alpha = 0.8,
             position = position_dodge()) +
    geom_errorbar(aes(ymin = value - visit_se, ymax = value + visit_se), 
                  width = 0.2, position = position_dodge(0.7), color = "#1D1616") + 
    labs(x = NULL, y = "Log fold change", 
         title = paste0("Differentially Abundant ", rank)) + 
    scale_fill_manual(name = "", values = c("#A19AD3", "#A31D1D")) +
    scale_color_manual(name = "", values = c("#A19AD3", "#A31D1D")) +
    coord_flip() +
    theme_bw() + 
    theme(
      plot.title = element_markdown(hjust = 0.5),
      axis.text = element_text(size = 14),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position = 'right'
    )
  
  ggsave(paste0("./export/sig_", rank, "_lf1_v7-v8.pdf"), 
         height = 6, width = 8, device = 'pdf')
  cat("Significant", rank, "plot saved.\n")
  
  return(p)
}

# MAIN EXECUTION WORKFLOW ----

main_analysis <- function() {
  # Create export directory if it doesn't exist

  
  cat("=== Starting ANCOM-BC Analysis Workflow ===\n\n")
  
  # PART 1: ASV-level analysis
  cat("PART 1: ASV-level Analysis\n")
  cat("==========================\n")
  
  # Load data
  data_list <- load_and_prepare_data()
  metadata <- data_list$metadata
  feature_table <- data_list$feature_table
  taxonomy_annot <- data_list$taxonomy_annot
  
  # Run ANCOM-BC
  results_asv <- run_ancombc_asv(feature_table, metadata)
  
  # Export normalized counts
  log_corr_abn <- export_normalized_counts(results_asv)
  
  # Process results and export statistics
  stat_df_asv <- process_ancombc_results(results_asv, taxonomy_annot, "ASV")
  

  cat("\nPART 2: Other taxonomic-level Aggregated Analysis\n")
  cat("========================================\n")
  
  # PART 2: other taxonomy-level analysis
  ranks <- ranks <- c("Phylum","Class","Order","Family","Genus")
  for (rank in ranks) {
    
    print(rank)
    # Aggregate features by genus
    feature_table_agg <- aggregate_features(feature_table, taxonomy_annot, rank)
    
    # Run ANCOM-BC on aggregated data
    results_genus <- run_ancombc_aggregated(feature_table_agg, metadata)
    
    # Export normalized counts for aggregated data
    export_normalized_counts(results_genus, paste0("_", rank))
    
    # Process results
    stat_df <- process_ancombc_results(results_genus, taxonomy_annot, rank)
  }

  
  cat("\n=== Analysis Complete ===\n")
  cat("Check the '", export_dir,"' directory for output files and plots.\n")
  
 
}


 main_analysis()