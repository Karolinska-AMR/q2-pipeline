# ==============================================================================
# EpiEndo Microbiome Analysis Pipeline
# Author: Mohammad Razavi
# Date: 2025-01-14
# ==============================================================================

# Load required libraries
library(dplyr)
library(qiime2R)
library(ggplot2)
library(ggtext)
library(tidyr)
library(tibble)
library(purrr)
library(ggsignif)
library(stringr)
library(cowplot)
library(ggside)
library(reshape2)
library(patchwork)
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

# ==============================================================================
# GLOBAL CONFIGURATIONS AND HELPER FUNCTIONS
# ==============================================================================

# Define visit labels mapping
VISIT_MAPPING <- c(
  "v2" = "Baseline 1",
  "v6" = "Baseline 2",
  "v7" = "1–week EP395",
  "v8" = "2–week EP395"
)

VISIT_LEVELS <<- c("Baseline 1", "Baseline 2", "1–week EP395", "2–week EP395")

# Define color schemes
VISIT_COLORS <<- c("#5CB338", "#FFB200", "#A19AD3", "#A31D1D")
VISIT_SHAPES <<- c(21, 24, 22, 23)

# Global function to load and process metadata
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


# ==============================================================================
# SECTION 1: ALPHA DIVERSITY ANALYSIS
# ==============================================================================

analyze_alpha_diversity <- function() {
  cat("\n=== ALPHA DIVERSITY ANALYSIS ===\n")
  
  # Clear environment and load data
  metadata <- load_metadata()
  
  cat("Loading alpha diversity metrics...\n")
  richness <- read_qza("./artifacts/diversities/faith_pd_vector.qza")$data
  richness <- richness %>% 
    mutate(divtype = "Richness (Faith's PD)") %>% 
    add_rownames(var = "sample.id") %>% 
    rename('value' = 'faith_pd')
  
  evenness <- read_qza("./artifacts/diversities/evenness_vector.qza")$data
  evenness <- evenness %>% 
    mutate(divtype = "Evenness") %>% 
    add_rownames(var = "sample.id") %>% 
    rename('value' = 'pielou_evenness')
  
  # Combine diversity metrics
  alpha_div_df <- rbind(richness, evenness) %>%
    merge(metadata[c("sample.id", "subject.id", "timepoint")], by = 'sample.id')
  
  
  # Save output
  fout<- file.path(export_dir, "alphadiversities.csv")
  write.csv(alpha_div_df, fout, row.names = FALSE)
  cat("Alpha diversity data saved to:",fout,"\n")
  
  # Perform statistical analysis
  cat("\nPerforming ANOVA and Tukey HSD tests...\n")
  annotation_stats <- compute_alpha_statistics(alpha_div_df)
  
  # Generate plots
  cat("\nGenerating alpha diversity plots...\n")
  plot_alpha_diversity(alpha_div_df, annotation_stats)
  
  cat("=== ALPHA DIVERSITY ANALYSIS COMPLETE ===\n")
  
  return(alpha_div_df)
}

# Function to compute alpha diversity statistics
compute_alpha_statistics <- function(data) {
  metrics <- c("Richness (Faith's PD)", "Evenness")
  p_threshold <- 0.05
  
  process_metric <- function(metric, data) {
    cat(sprintf("  Processing metric: %s\n", metric))
    
    df_filter <- data %>% filter(divtype == metric)
    max_val <- max(df_filter$value)
    y_step <- 0.05 * max_val
    
    # Perform ANOVA and Tukey HSD
    anova_result <- aov(value ~ timepoint, data = df_filter)
    tukey_result <- TukeyHSD(anova_result, p.adjust.method = "fdr")

        # Process Tukey results
    stat <- as.data.frame(tukey_result$timepoint) %>%
      rownames_to_column(var = "group") %>%
      mutate(
        split_groups = strsplit(group, "-"),
        end = sapply(split_groups, `[`, 1),
        start = sapply(split_groups, `[`, 2)
      ) %>%
      filter(start == 'Baseline 1') %>%
      mutate(
        divtype = metric,
        y = max_val + (row_number() * y_step),
        timepoint = start,
        label = paste0("p=", round(`p adj`, 3))
      )
    
    if (nrow(stat) == 0) {
      cat(sprintf("    No significant results for %s\n", metric))
      return(NULL)
    }
    
    cat(sprintf("    Found %d significant comparisons\n", nrow(stat)))
    return(stat)
  }
  
  annotation_stats <- map_dfr(metrics, process_metric, data = data, .id = NULL)
  return(annotation_stats)
}

# Function to plot alpha diversity
plot_alpha_diversity <- function(data, stat) {

  plot_alpha_div <- function(data, stat, metric_name, hide_axis_x = FALSE) {
    custom_theme <- theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 12),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        axis.text.x = element_text(size = 11, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 13),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        panel.grid.major.x = element_line(linewidth = 0.3, linetype = "solid"),
        panel.grid.major.y = element_line(linetype = "dashed", colour = "#bdbdbd"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "#bdbdbd"),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
      )
    print(head(data))
    p <- ggplot(data, aes(x = timepoint, y = value, fill = timepoint)) +
      geom_dotplot(
        shape = 21,
        binaxis = "y",
        stackdir = "center",
        dotsize = 0.6,
        alpha = 0.8,
        position = position_nudge(x = 0)
      ) +
      geom_boxplot(
        width = 0.3,
        outlier.color = NA,
        alpha = 0.65
      ) +
      geom_signif(
        data = stat,
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        textsize = 5,
        vjust = -0.2,
        manual = TRUE
      ) +
      scale_fill_manual(name = "", values = VISIT_COLORS) +
      scale_color_manual(name = "", values = VISIT_COLORS) +
      custom_theme
    
    x_label <- ifelse(hide_axis_x, "", "Time Point")
    p <- p + labs(y = metric_name, x = x_label)
    
    return(p)
  }
  
  # Create plots for each metric
  p1 <- plot_alpha_div(
    data %>% filter(str_detect(divtype, "^Richness")),
    stat %>% filter(str_detect(divtype, "^Richness")),
    "Richness (Faith's PD)", TRUE
  )
  
  p2 <- plot_alpha_div(
    data %>% filter(str_detect(divtype, "^Evenness")),
    stat %>% filter(str_detect(divtype, "^Evenness")),
    "Evenness"
  )
  
  # Combine plots
  combined_plot <- plot_grid(p1, p2, labels = "AUTO", ncol = 1)
  
  # Save plot
  fout <- file.path(export_dir, "alpha_diversity.pdf")
  ggsave(fout, 
         plot = combined_plot, 
         height = 12, 
         width = 6, 
         device = 'pdf')
  
  cat("Alpha diversity plot saved to: ",fout,"\n")
}

# ==============================================================================
# SECTION 2: BETA DIVERSITY (PCoA) ANALYSIS
# ==============================================================================

analyze_beta_diversity <- function() {
  cat("\n=== BETA DIVERSITY (PCoA) ANALYSIS ===\n")
  
  # Load data
  metadata <- load_metadata()
  
  cat("Loading PCoA results...\n")
  pcoaResults <- read_qza("./artifacts/diversities/weighted_unifrac_pcoa_results.qza")
  pca_tbl <- pcoaResults$data$Vectors
  
  # Calculate variance proportions
  pc_variances <- pca_tbl %>%
    summarize(across(everything(), var))
  total_variance <- sum(pc_variances, na.rm = TRUE)
  proportion_variance <- round(100.0 * pc_variances / total_variance, 2)
  
  cat(sprintf("PC1 explains %.2f%% of variance\n", proportion_variance$PC1))
  cat(sprintf("PC2 explains %.2f%% of variance\n", proportion_variance$PC2))
  
  print(head(metadata))
  # Merge with metadata
  pca_tbl <- merge(
    x = pca_tbl %>% select('SampleID', 'PC1', 'PC2'),
    y = metadata,
    by.x = 'SampleID',by.y = 'sample.id'
  ) %>%
    mutate(timepoint = factor(timepoint, levels = VISIT_LEVELS))
  
  names(pca_tbl)
  str(pca_tbl)
  summary(pca_tbl$PC1)
  summary(pca_tbl$PC2)
  head(pca_tbl)
  
  # Create PCoA plot
  cat("Generating PCoA scatter plot...\n")
  p <- ggplot(data = pca_tbl, aes(shape = timepoint, x = PC2, y = PC1)) +
    geom_point(aes(fill = timepoint), size = 5, alpha = 0.7, stroke = 0.5, color = 'black') +
    scale_shape_manual(name = "", values = VISIT_SHAPES) +
    scale_fill_manual(name = "", values = VISIT_COLORS) +
    labs(
      y = paste0("PCoA1 [", proportion_variance$PC1, "%]"),
      x = paste0("PCoA2 [", proportion_variance$PC2, "%]")
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 12),
      legend.text = element_markdown(size = 14),
      legend.title = element_markdown(size = 10),
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.title.x = element_text(margin = margin(t = -50, unit = "pt")),
      axis.title.y = element_text(margin = margin(r = -10, unit = "pt"))
    ) +
  geom_ysideboxplot(aes(x = timepoint, y = PC1, fill = timepoint),  orientation = "x",# Local aes: use PC1 on shared y
                    inherit.aes = FALSE, varwidth = TRUE, alpha = 0.7, show.legend = FALSE) +
    geom_xsideboxplot(aes(y = timepoint, x = PC2, fill = timepoint), orientation = "y", # Local aes: use PC2 on shared x
                      inherit.aes = FALSE, varwidth = TRUE, alpha = 0.7, show.legend = FALSE) +
    theme(ggside.axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust = 0)) +
    guides(shape = guide_legend(ncol = 4))
  
  # Save plot
  fout <- file.path(export_dir, "pcoa_scatter.pdf")
  ggsave(fout, 
         plot = p, 
         height = 10, 
         width = 10, 
         device = 'pdf')
  
  cat("PCoA plot saved to: ",fout,"\n")
  cat("=== BETA DIVERSITY ANALYSIS COMPLETE ===\n")
  
  return(list(d=pca_tbl,fig=p))
}

# ==============================================================================
# SECTION 3: DIFFERENTIAL ABUNDANCE ANALYSIS
# ==============================================================================

analyze_differential_abundance <- function() {
  cat("\n=== DIFFERENTIAL ABUNDANCE ANALYSIS ===\n")
  
  # Load data
  metadata <- load_metadata()
  
  cat("Loading taxonomy and differential abundance data...\n")
  taxonomy_annot <- read_qza('./taxonomy_classification/taxonomy.qza')$data
  taxonomy_annot <- taxonomy_annot %>% 
    mutate(taxa = str_remove_all(Taxon, "D_\\d__")) %>%
    separate_wider_delim(
      taxa, "; ",
      names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      too_few = "align_start"
    ) %>%
    mutate_at(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
              list(~ substr(., 4, nchar(.))))
  
  diff_otus <- read.csv("./R_export/ancombc_ASVs_all_ci95.csv")
  log_nrm_counts <- read.csv("./R_export/log_normalized_otu_counts.csv")
  
  cat(sprintf("Loaded %d differential OTUs\n", nrow(diff_otus)))
  
  # Get top differential OTUs
  diff_w1w2 <- rbind(
    get_top_otus(diff_otus, 'week1_ep395', n = 200, lfc_threshold = 1),
    get_top_otus(diff_otus, 'week2_ep395', n = 200, lfc_threshold = 1)
  )
  
  cat(sprintf("Selected %d unique differential OTUs\n", length(unique(diff_w1w2$feature_id))))
  
  # Generate heatmap
  cat("Generating differential abundance heatmap...\n")
  plot_differential_heatmap(diff_w1w2, log_nrm_counts, metadata, taxonomy_annot)
  
  # Generate taxonomic rank plots
  cat("Generating taxonomic rank plots...\n")
  plot_taxonomic_ranks()
  
  cat("=== DIFFERENTIAL ABUNDANCE ANALYSIS COMPLETE ===\n")
}

# Function to get top differential OTUs
get_top_otus <- function(df, group, n = 10, lfc_threshold = 1, id = 'feature_id') {
  rank_col <- paste0("lfc_", group)
  padj_col <- paste0("q_", group)
  
  top_otus <- df %>% 
    filter(.data[[padj_col]] < 0.05) %>%
    filter(abs(.data[[rank_col]]) >= lfc_threshold) %>%
    mutate(group_label = group, lfc = .data[[rank_col]]) %>%
    select(all_of(id), lfc, group_label) %>%
    arrange(desc(abs(lfc))) %>% 
    slice_head(n = n)
  
  if (nrow(top_otus) < n) {
    cat(sprintf("  Warning: Only %d OTUs met criteria for group %s (requested n = %d)\n", 
                nrow(top_otus), group, n))
  } else {
    cat(sprintf("  Found %d top OTUs for group %s\n", nrow(top_otus), group))
  }
  
  return(top_otus)
}

# Function to plot differential abundance heatmap
plot_differential_heatmap <- function(diff_v7v8, log_nrm_counts, metadata, taxonomy_annot) {
  # Select and aggregate counts
  sel_counts <- log_nrm_counts %>% 
    filter(X %in% diff_v7v8$feature_id) %>% 
    t()
  colnames(sel_counts) <- sel_counts[1, ]
  sel_counts <- sel_counts[-1, , drop = FALSE]
  
  agg_counts <- sel_counts %>%
    merge(metadata %>% select("sample.id", "timepoint"), by.y = "sample.id", by.x = "row.names") %>%
    mutate(across(!starts_with("Row.names") & !starts_with('timepoint'), as.numeric)) %>%
    group_by(timepoint) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    t() %>%
    as.data.frame()
  


  colnames(agg_counts) <- agg_counts[1, ]
  agg_counts <- agg_counts[-1, , drop = FALSE]
  
  print(head(agg_counts))
  
  # Prepare data for plotting
  agg_counts_mlt <- agg_counts %>% 
    rownames_to_column(var = "fid") %>%
    pivot_longer(
      cols = -fid,
      names_to = "timepoint",
      values_to = "count"
    ) %>% 
    merge(taxonomy_annot, by.y = 'Feature.ID', by.x = "fid") %>% 
    mutate(count = as.numeric(count))
  
  # Create labels
  labels <- agg_counts_mlt %>% 
    distinct() %>%
    mutate(
      otu_name = paste0('ASV_', str_sub(fid, -4)),
      species = coalesce(Species, otu_name),
      taxon_label = paste0("***", Genus, "***:<span style='font-size:10px'> ", species, "</span>")
    ) %>%
    select(fid, taxon_label)
  
  fid_to_genus_species <- set_names(labels$taxon_label, labels$fid)
  
  agg_counts_mlt <- agg_counts_mlt %>%
    mutate(
      timepoint = factor(timepoint, levels = VISIT_LEVELS)
    ) %>%
    mutate(
      fid = factor(
        fid,
        levels = unique(fid[order(Genus, decreasing = TRUE)])
      )
    )
  
  # Create heatmap
  p <- ggplot(agg_counts_mlt, aes(timepoint, fid, fill = count)) +
    geom_tile(color = "#bdbdbd", lwd = 0.2, linetype = 1) +
    scale_y_discrete(labels = fid_to_genus_species) +
    scale_fill_distiller(name = "Log(Normalized Abundance)", palette = "RdPu", direction = 1) +
    labs(y = "**ASVs**, Genus: Species", x = "") +
    theme(
      legend.position = "bottom",
      legend.title = element_text(vjust = 0.8),
      axis.text.y = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(size = rel(1.2), angle = 0, hjust = 0.5, vjust = 0.5)
    )
  fout = file.path(export_dir,'diff_ASVs_week1-week2_lfc1.pdf')
  ggsave(fout, 
         plot = p, 
         height = 15, 
         width = 9, 
         device = 'pdf')
  
  cat("Heatmap saved to: ",fout,"\n")
}

# Function to load taxonomic statistics
load_tx_stat <- function(rank) {
  
  stat_df <- read.csv(file.path(export_dir,paste0('ancombc_', rank, '_all_ci95.csv')))
  
  df_filter <- stat_df %>%
    filter((diff_week1_ep395 == TRUE & abs(lfc_week1_ep395) > 1) | 
             (diff_week2_ep395 == TRUE & abs(lfc_week2_ep395) > 1)) %>%
    select(all_of(rank), 'lfc_week2_ep395', 'lfc_week1_ep395', "se_week2_ep395", "se_week1_ep395") %>%
    rename('w2' = 'lfc_week2_ep395', 'w1' = 'lfc_week1_ep395', 
           'se_w2' = "se_week2_ep395", 'se_w1' = "se_week1_ep395") %>%
    rename('taxa_lvl' := !!rank)
  
  melt_df <- bind_cols(
    df_filter %>% select(taxa_lvl, w1, w2) %>% melt() %>%
      rename('timepoint' = 'variable') %>% select(taxa_lvl, timepoint, value),
    df_filter %>% select(taxa_lvl, se_w1, se_w2) %>% melt() %>%
      rename('timepoint_se' = 'value') %>% select(timepoint_se)
  ) %>%
    mutate(rank = rank)
  
  melt_df$taxa_lvl <- factor(
    melt_df$taxa_lvl, 
    levels = unique(melt_df[order(melt_df$value), "taxa_lvl"])
  )
  
  cat(sprintf("  Loaded %d differential %s\n", 
              length(unique(melt_df$taxa_lvl)), tolower(rank)))
  
  return(melt_df)
}

# Function to plot taxonomic bars
plot_taxonomic_ranks <- function() {
  plot_bars <- function(data, title) {

    visit_labels <- c("w1"="1–week EP395", "w2"="2–week EP395")
    color_labels <- c("w1"="#A19AD3", "w2"="#A31D1D")
    p <- ggplot(data = data, 
                aes(x = taxa_lvl, y = value, fill = timepoint, color = timepoint, group = timepoint)) +
      geom_bar(stat = "identity", width = 0.7, color = "#5E686D", alpha = 0.8,
               position = position_dodge()) +
      geom_errorbar(aes(ymin = value - timepoint_se, ymax = value + timepoint_se), 
                    width = 0.2,
                    position = position_dodge(0.7), color = "#1D1616") +
      labs(x = NULL, y = "LFC") +
      scale_fill_manual(name = "", values = color_labels,labels=visit_labels) +
      scale_color_manual(name = "", values = color_labels,labels=visit_labels) +
      ylim(c(-4.2, 1.1)) +
      coord_flip() +
      theme_bw() +
      theme(
        plot.title = element_blank(),
        axis.text = element_text(size = rel(1)),
        axis.text.y = element_text(face = 'bold.italic', size = rel(0.7)),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = 'bottom'
      )
    
    return(p)
  }
  
  # Load and plot each taxonomic rank
  p1 <- plot_bars(load_tx_stat('Genus'), "Genera")
  p2 <- plot_bars(load_tx_stat('Family'), "Families")
  p3 <- plot_bars(load_tx_stat('Class'), "Classes")
  p4 <- plot_bars(load_tx_stat('Phylum'), "Phyla")
  
 
  # Combine plots
  col1 <- p1/p3 + plot_layout(heights = c(1,0.3))
  col2 <- p2/p4 + plot_layout(heights = c(1,0.3))
  
  combined <- (col1 | col2) & theme(legend.position = "bottom")
  combined + plot_layout(guides = "collect") +  plot_annotation(tag_levels = 'A')
  
               
  fout = file.path(export_dir,'diff_otus_w1w2_lfc1_taxonomicRank.pdf')
  ggsave(fout, 
         height = 8, 
         width = 9, 
         device = 'pdf')
  
  cat("Taxonomic rank plots saved to: ", fout,"\n")
}

## Create visualization for significant ASVs ----
# plot_significant_asvs <- function(results, taxonomy_annot, log_corr_abn, metadata, lfc_threshold = 1) {
#   res <- results$res
#   
#   # Extract results tables
#   tab_lfc <- res$lfc %>% rename_with(~ paste0("lfc_", .x), -taxon) %>% 
#     merge(taxonomy_annot, by.x = 'taxon', by.y = 'Feature.ID')
#   tab_q <- res$q %>% rename_with(~ paste0("q_", .x), -taxon)
#   tab_diff <- res$diff_abn %>% rename_with(~ paste0("diff_", .x), -taxon)
#   
#   # Filter significant ASVs
#   sig_taxa <- tab_lfc %>% 
#     merge(tab_q, by = 'taxon') %>% 
#     merge(tab_diff, by = 'taxon') %>%
#     filter((diff_visitv8 == TRUE & abs(lfc_visitv8) > lfc_threshold) | 
#              (diff_visitv7 == TRUE & abs(lfc_visitv7) > lfc_threshold))
#   
#   if (nrow(sig_taxa) == 0) {
#     cat("No significant taxa found with LFC threshold of", lfc_threshold, "\n")
#     return(NULL)
#   }
#   
#   # Prepare data for plotting
#   sigcount_data <- log_corr_abn[sig_taxa$taxon, ] %>% 
#     melt() %>% 
#     rename('taxon' = 'Var1', 'SampleID' = 'Var2', 'ncount' = 'value') %>%
#     merge(taxonomy_annot[, c('Feature.ID', 'Family', 'Genus')], 
#           by.x = 'taxon', by.y = 'Feature.ID') %>%
#     merge(metadata, by = 'SampleID') %>% 
#     arrange(Family)
#   
#   # Filter families with sufficient taxa
#   family_counts <- sig_taxa %>% 
#     group_by(Family) %>% 
#     summarise(fcount = n()) %>% 
#     filter(fcount > 4)
#   
#   sigcount_data <- sigcount_data %>% 
#     filter(Family %in% family_counts$Family)
#   
#   if (nrow(sigcount_data) == 0) {
#     cat("No families with >4 significant taxa found.\n")
#     return(NULL)
#   }
#   
#   # Create color palette
#   n_colors <- length(unique(sigcount_data$taxon))
#   colors_palette <- colorRampPalette(c(glasbey(), kelly(), alphabet()))(n_colors)
#   color_clustering <- dist(t(col2rgb(colors_palette))) %>% hclust()
#   reshuffled_colors <- colors_palette[color_clustering$order]
#   
#   # Create plot
#   p <- ggboxplot(sigcount_data, x = 'visit', y = 'ncount', 
#                  fill = 'taxon', color = 'taxon', width = 0.5, alpha = 0.5) +
#     facet_grid(Family ~ ., switch = 'y') +
#     scale_fill_manual(name = "", values = reshuffled_colors) +
#     scale_color_manual(name = "", values = reshuffled_colors) +
#     labs(y = 'log(Normalized Abundance)', x = 'Visit') +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "none",
#       axis.text.x = element_text(size = 12, hjust = 0.5, vjust = 0.5),
#       axis.text.y = element_text(size = 12),
#       axis.title.y = element_text(size = 14),
#       axis.title.x = element_text(size = 13),
#       strip.background = element_blank(),
#       strip.placement = "outside",
#       strip.text = element_text(size = 10),
#       panel.grid.major.x = element_line(linewidth = 0.1, linetype = "solid"),
#       panel.grid.major.y = element_line(linewidth = 0.1, linetype = "dashed", colour = "gray80"),
#       panel.grid.minor.y = element_line(linewidth = 0.1, linetype = "dotted", colour = "gray80")
#     )
#   
#   ggsave("./export/sig_ASVs_gfamily_lf1_v7-v8_c5.pdf", height = 13, width = 10, device = 'pdf')
#   cat("Significant ASVs plot saved.\n")
#   
#   return(p)
# }

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

run_complete_analysis <- function() {
  cat("\n")
  cat("===============================================================================\n")
  cat("                    EpiEndo Microbiome Analysis Pipeline                       \n")
  cat("===============================================================================\n")
  
  # Run all analyses
  alpha_results <- analyze_alpha_diversity()
  beta_results <- analyze_beta_diversity()
  analyze_differential_abundance()
  
  cat("\n")
  cat("===============================================================================\n")
  cat("                         ANALYSIS COMPLETE                                      \n")
  cat("===============================================================================\n")
  cat("\nAll figures and supplementary files have been saved to the report directory.\n")
  
  return(list(
    alpha = alpha_results,
    beta = beta_results
  ))
}

# Execute the complete analysis pipeline

#alpha_results <- analyze_alpha_diversity()
#beta_results <- analyze_beta_diversity()
#analyze_differential_abundance()



results <- run_complete_analysis()