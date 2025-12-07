#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------
# Defaults and globals
# ------------------------------
DEPTH=74470
THREADS=4
GROUP_COL="timepoint"

# Paths (empty until parsed)
FASTQ_DIR=""
SAMPLE_MAP=""
METADATA=""
OUTPUT_DIR=""
CLASSIFIER_DIR=""

# Helper for colored output
info() { printf "\e[1;32m%s\e[0m\n" "$*"; }
warn() { printf "\e[1;33m%s\e[0m\n" "$*" >&2; }
err() { printf "\e[1;31m%s\e[0m\n" "$*" >&2; exit 1; }
step() { printf "\e[1;34mStep %s: %s\e[0m\n" "$1" "$2"; }

usage() {
    cat <<EOF
Usage: $0 -f FASTQ_DIR -s SAMPLE_MAP -m METADATA -o OUTPUT_DIR -c CLASSIFIER_DIR [options]

Required:
  -f PATH    directory with FASTQ files
  -s FILE    sample mapping file (CSV or TSV). Header expected; first column sample filename prefix, second column sample-id/name
  -m FILE    metadata file (QIIME2 compatible)
  -o DIR     output directory (will be created)
  -c DIR     classifier directory (will pick first .qza)

Options:
  -d INT     sampling depth (default: ${DEPTH})
  -t INT     threads (default: ${THREADS})
  -g STR     metadata column for PERMANOVA (default: ${GROUP_COL})
  -h         show help
EOF
    exit 1
}

# ------------------------------
# Argument parsing
# ------------------------------
while getopts ":f:s:m:o:c:d:t:g:h" opt; do
    case "$opt" in
        f) FASTQ_DIR="$OPTARG" ;;
        s) SAMPLE_MAP="$OPTARG" ;;
        m) METADATA="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CLASSIFIER_DIR="$OPTARG" ;;
        d) DEPTH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        g) GROUP_COL="$OPTARG" ;;
        h) usage ;;
        :) err "Missing argument for -$OPTARG" ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

# Validate required args
[[ -z "$FASTQ_DIR" || -z "$SAMPLE_MAP" || -z "$METADATA" || -z "$OUTPUT_DIR" || -z "$CLASSIFIER_DIR" ]] && usage

mkdir -p "$OUTPUT_DIR"
# Normalize to absolute paths where useful
FASTQ_DIR=$(realpath "$FASTQ_DIR")
SAMPLE_MAP=$(realpath "$SAMPLE_MAP")
METADATA=$(realpath "$METADATA")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
CLASSIFIER_DIR=$(realpath "$CLASSIFIER_DIR")

# Make output dirs
mkdir -p "$OUTPUT_DIR"/{artifacts,taxonomy_classification,statistical_tests}

# Find classifier
CLASSIFIER_FILE=$(find "$CLASSIFIER_DIR" -maxdepth 1 -type f -name "*.qza" | sort | head -n1 || true)
[[ -z "$CLASSIFIER_FILE" ]] && err "No classifier .qza found in $CLASSIFIER_DIR"

# ------------------------------
# Header & status
# ------------------------------
cat <<-EOF
IMPORTANT: Default parameters are optimized for microbiome datasets processed at Karolinska Institutet.
===== QIIME2 Microbiome Analysis Workflow =====
Input FASTQ directory: $FASTQ_DIR
Sample map file: $SAMPLE_MAP
Metadata file: $METADATA
Output directory: $OUTPUT_DIR
Classifier file: $CLASSIFIER_FILE
Sampling depth: $DEPTH
Threads: $THREADS
Metadata column for PERMANOVA: $GROUP_COL
===============================================
EOF

# if false; then ## IGNORE - jump to R analysis, for testing purposes
# ------------------------------
# Utility wrappers
# ------------------------------
run_qiime() {
    # Pass args directly to qiime; log and fail-on-error due to set -e
    step "QIIME" "$*"
    qiime "$@"
}


# ------------------------------
# Create manifest from sample map (supports CSV or TSV)
# Expected sample map: header + rows where first column is filename prefix (used to find *_1.fq.gz and *_2.fq.gz)
# second column is sample-id to write to manifest
# ------------------------------
MANIFEST="$OUTPUT_DIR/manifest_v2.csv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"

# detect delimiter from first non-empty line (after optional header)
first_line=$(sed -n '1p' "$SAMPLE_MAP" || true)
if [[ "$first_line" == *","* && "$first_line" != *$'\t'* ]]; then DELIM=","; else DELIM=$'\t'; fi

# Read sample map skipping the header (tail -n +2) safely
while IFS="$DELIM" read -r col1 col2 _rest || [[ -n "$col1" ]]; do
    # skip header if it looks like header (contains letters not digits and not typical filename characters)
    if [[ "$col1" =~ ^(#|sample|Sample|sample-id) ]]; then
        continue
    fi

    sample_prefix=$(echo "$col1" | xargs)
    sample_id=$(echo "${col2:-$col1}" | xargs)

    [[ -z "$sample_prefix" || -z "$sample_id" ]] && { warn "Skipping malformed line: $col1"; continue; }
   
    f1=$(find "$FASTQ_DIR" -type f -name "${sample_prefix}_1.fq.gz" | head -n1 || true)
    f2=$(find "$FASTQ_DIR" -type f -name "${sample_prefix}_2.fq.gz" | head -n1 || true)

    if [[ -n "$f1" && -n "$f2" ]]; then
        printf "%s\t%s\t%s\n" "$sample_id" "$f1" "$f2" >> "$MANIFEST"
    else
        warn "Warning: Missing files for filename '${sample_prefix}'. Skipping."
    fi

done < <(tail -n +2 "$SAMPLE_MAP" | sed '/^[[:space:]]*$/d')

info "Manifest written to $MANIFEST"

# ------------------------------
# Run QIIME pipeline steps (straightforward, no eval)
# ------------------------------

# 1. Import
step 1 "Importing sequence data"
run_qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST" \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path "$OUTPUT_DIR/artifacts/imported_reads_demux.qza"

# 2. Merge pairs
step 2 "Merging paired-end reads"
run_qiime vsearch merge-pairs \
    --i-demultiplexed-seqs "$OUTPUT_DIR/artifacts/imported_reads_demux.qza" \
    --o-merged-sequences "$OUTPUT_DIR/artifacts/demux-merged.qza" \
    --o-unmerged-sequences "$OUTPUT_DIR/artifacts/dmux-unmerged.qza" \
    --p-threads "$THREADS"

# 3. Summarize
step 3 "Summarizing demultiplexed data"
run_qiime demux summarize \
    --i-data "$OUTPUT_DIR/artifacts/demux-merged.qza" \
    --o-visualization "$OUTPUT_DIR/artifacts/demux-merged.qzv"

# 4. Quality filter
step 4 "Quality filtering"
run_qiime quality-filter q-score \
    --i-demux "$OUTPUT_DIR/artifacts/demux-merged.qza" \
    --o-filtered-sequences "$OUTPUT_DIR/artifacts/demux-merged.filtered.qza" \
    --o-filter-stats "$OUTPUT_DIR/artifacts/demux-merged.filter_stats.qza" \
    --p-min-quality 30

# 5. Tabulate filter stats
step 5 "Tabulating filter stats"
run_qiime metadata tabulate \
    --m-input-file "$OUTPUT_DIR/artifacts/demux-merged.filter_stats.qza" \
    --o-visualization "$OUTPUT_DIR/artifacts/demux-merged.filter_stats.qzv"

# 6. Deblur denoising
step 6 "Deblur denoising"
run_qiime deblur denoise-16S \
    --i-demultiplexed-seqs "$OUTPUT_DIR/artifacts/demux-merged.filtered.qza" \
    --p-trim-length 400 \
    --p-sample-stats \
    --o-representative-sequences "$OUTPUT_DIR/artifacts/rep-seqs.qza" \
    --o-table "$OUTPUT_DIR/artifacts/table.qza" \
    --o-stats "$OUTPUT_DIR/artifacts/deblur-stats.qza"

# 7. Feature table summary
step 7 "Summarizing feature table"
run_qiime feature-table summarize \
    --i-table "$OUTPUT_DIR/artifacts/table.qza" \
    --o-visualization "$OUTPUT_DIR/artifacts/table.qzv" \
    --m-sample-metadata-file "$METADATA"

# 8. Deblur stats visualize
step 8 "Visualizing deblur stats"
run_qiime deblur visualize-stats \
    --i-deblur-stats "$OUTPUT_DIR/artifacts/deblur-stats.qza" \
    --o-visualization "$OUTPUT_DIR/artifacts/deblur-stats.qzv"

# 9. Phylogenetic tree
step 9 "Constructing phylogenetic tree"
run_qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$OUTPUT_DIR/artifacts/rep-seqs.qza" \
    --o-alignment "$OUTPUT_DIR/artifacts/aligned-rep-seqs.qza" \
    --o-masked-alignment "$OUTPUT_DIR/artifacts/masked-aligned-rep-seqs.qza" \
    --o-tree "$OUTPUT_DIR/artifacts/unrooted-tree.qza" \
    --o-rooted-tree "$OUTPUT_DIR/artifacts/rooted-tree.qza" \
    --p-n-threads "$THREADS"

# 10. Diversity core-metrics
step 10 "Running diversity analyses"
run_qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$OUTPUT_DIR/artifacts/rooted-tree.qza" \
    --i-table "$OUTPUT_DIR/artifacts/table.qza" \
    --p-sampling-depth "$DEPTH" \
    --m-metadata-file "$METADATA" \
    --output-dir "$OUTPUT_DIR/artifacts/diversities"

# 11. Assign taxonomy
step 11 "Assigning taxonomy"
run_qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER_FILE" \
    --i-reads "$OUTPUT_DIR/artifacts/rep-seqs.qza" \
    --o-classification "$OUTPUT_DIR/taxonomy_classification/taxonomy.qza" \
    --p-n-jobs "$THREADS"

# 12. Visualize taxonomy
step 12 "Visualizing taxonomy"
run_qiime metadata tabulate \
    --m-input-file "$OUTPUT_DIR/taxonomy_classification/taxonomy.qza" \
    --o-visualization "$OUTPUT_DIR/taxonomy_classification/taxonomy.qzv"

# 13. Taxa barplot
step 13 "Creating taxonomy barplot"
run_qiime taxa barplot \
    --i-table "$OUTPUT_DIR/artifacts/table.qza" \
    --i-taxonomy "$OUTPUT_DIR/taxonomy_classification/taxonomy.qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$OUTPUT_DIR/taxonomy_classification/taxa-bar-plots.qzv"

# 14. PERMANOVA tests for distance matrices
step 14 "Running PERMANOVA tests"
# list of expected distance matrix file bases
DIST_NAMES=(unweighted_unifrac_distance_matrix weighted_unifrac_distance_matrix jaccard_distance_matrix bray_curtis_distance_matrix)
for dbase in "${DIST_NAMES[@]}"; do
    dm_qza="$OUTPUT_DIR/artifacts/diversities/${dbase}.qza"
    out_viz="$OUTPUT_DIR/statistical_tests/permanova_${dbase}.qzv"
    if [[ -f "$dm_qza" ]]; then
        run_qiime diversity beta-group-significance \
            --i-distance-matrix "$dm_qza" \
            --m-metadata-file "$METADATA" \
            --m-metadata-column "$GROUP_COL" \
            --o-visualization "$out_viz" \
            --p-pairwise
    else
        warn "Distance matrix not found: $dm_qza (skipping)"
    fi
done


# 15. Alpha diversity tests
step 15 "Running alpha diversity significance tests"
ALPHA_METRICS=(faith_pd observed_features shannon evenness)
for metric in "${ALPHA_METRICS[@]}"; do
    vec_qza="$OUTPUT_DIR/artifacts/diversities/${metric}_vector.qza"
    out_qzv="$OUTPUT_DIR/statistical_tests/alpha_significance_${metric}.qzv"
    if [[ -f "$vec_qza" ]]; then
        run_qiime diversity alpha-group-significance \
            --i-alpha-diversity "$vec_qza" \
            --m-metadata-file "$METADATA" \
            --o-visualization "$out_qzv"
    else
        warn "Alpha diversity vector not found: $vec_qza (skipping)"
    fi
done

# 16. Emperor PCoA plots
step 16 "Generating Emperor PCoA plots with grouping"
# Look for pcoa results corresponding to distance names (strip suffix _distance_matrix)
for dbase in "${DIST_NAMES[@]}"; do
    pcoa_base=${dbase%_distance_matrix}
    pcoa_qza="$OUTPUT_DIR/artifacts/diversities/${pcoa_base}_pcoa_results.qza"
    out_qzv="$OUTPUT_DIR/artifacts/diversities/${pcoa_base}_emperor.qzv"
    if [[ -f "$pcoa_qza" ]]; then
        run_qiime emperor plot \
            --i-pcoa "$pcoa_qza" \
            --m-metadata-file "$METADATA" \
            --o-visualization "$out_qzv"
    else
        warn "PCoA results not found for ${pcoa_base} (expected $pcoa_qza)"
    fi
done

# 17. Summary report
step 17 "Generating summary report of statistical tests"
summary="$OUTPUT_DIR/statistical_tests/summary_report.txt"
{
    echo "QIIME2 Statistical Tests Summary"
    echo "================================"
    echo
    echo "Analysis performed on: $(date --iso-8601=seconds)"
    echo "Metadata column used for grouping: $GROUP_COL"
    echo "Sampling depth: $DEPTH"
    echo
    echo "PERMANOVA Tests:"
    echo "---------------"
    printf -- "- %s\n" "Unweighted UniFrac" "Weighted UniFrac" "Jaccard" "Bray-Curtis"
    echo
    echo "Alpha Diversity Tests:"
    echo "--------------------"
    printf -- "- %s\n" "faith_pd" "observed_features" "shannon" "evenness"
    echo
    echo "To view results, open the .qzv files in QIIME2 View: https://view.qiime2.org/"
} > "$summary"

info "===== QIIME2 analysis with statistical tests complete! ====="
info "Output artifacts are in: $OUTPUT_DIR"
info "Statistical test results are in: $OUTPUT_DIR/statistical_tests"
info "To view .qzv files, use 'qiime tools view file.qzv'"

#fi ## End of 'if false' block

# ------------------------------
# Downstream analysis with R
# ------------------------------


# 18. Run R downstream analysis script
step 18 "Running R downstream analysis: ancombc.R"
Rscript "$SCRIPT_DIR/ancombc.R" "$OUTPUT_DIR" "$METADATA"

# 19. Run R downstream analysis script
step 19 "Running R downstream analysis: plot_manager.R"
Rscript "$SCRIPT_DIR/plot_manager.R" "$OUTPUT_DIR" "$METADATA"
exit 0
