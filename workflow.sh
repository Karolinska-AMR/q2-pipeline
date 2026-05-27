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
BRANCHES="all"       # original | raspergade | decontam | all
START_STEP=1         # jump to this step; earlier steps are checkpointed
FORCE=false          # -F: re-run even if output artifact already exists
QUANT_COL=""         # metadata column with DNA concentrations (decontam branch only)

FASTQ_DIR=""
SAMPLE_MAP=""
METADATA=""
OUTPUT_DIR=""
CLASSIFIER_DIR=""

# ------------------------------
# Logging helpers (same style as original workflow)
# ------------------------------
info()    { printf "\e[1;32m%s\e[0m\n"                  "$*"; }
warn()    { printf "\e[1;33mWARN: %s\e[0m\n"            "$*" >&2; }
err()     { printf "\e[1;31mERROR: %s\e[0m\n"           "$*" >&2; exit 1; }
step()    { printf "\e[1;34mStep %s: %s\e[0m\n"         "$1" "$2"; }
bhead()   { printf "\e[1;35m\n>>> Branch: %s <<<\e[0m\n" "$*"; }
skipped() { printf "\e[0;90m  [checkpoint] %-55s already exists\e[0m\n" "$*"; }
qrun()    { step "QIIME" "$*"; qiime "$@"; }

# ------------------------------
# Usage / step map
# ------------------------------
usage() {
    cat <<EOF
Usage: $0 -f FASTQ_DIR -s SAMPLE_MAP -m METADATA -o OUTPUT_DIR -c CLASSIFIER_DIR [options]

Required:
  -f PATH    directory with FASTQ files
  -s FILE    sample map (CSV or TSV); col1=filename prefix, col2=sample-id
  -m FILE    QIIME 2 compatible metadata file
  -o DIR     output directory (will be created)
  -c DIR     classifier directory (first .qza found is used)

Options:
  -d INT     sampling depth for rarefaction        (default: ${DEPTH})
  -t INT     threads                               (default: ${THREADS})
  -g STR     metadata column for PERMANOVA         (default: ${GROUP_COL})
  -b STR     branch(es) to run                     (default: ${BRANCHES})
               original    full original pipeline (steps 1–19)
               raspergade  RasperGade16S correction (steps 1–9 shared, 20–27)
               decontam    decontam frequency       (steps 1–9 shared, 28–33; needs -q)
               all         run all applicable branches
  -q STR     metadata column with DNA concentrations (required for decontam branch)
  -S INT     start from step N; earlier steps are skipped but their output
             artifacts are validated; use to resume after failure or jump to
             a specific stage (default: ${ -})
  -F         force re-run even if output artifacts already exist
  -h         show this help

Step map (for -S):
  ── Shared pre-processing (all branches) ────────────────────────────────────
   1  Build manifest + import reads
   2  Merge paired-end reads (vsearch)
   3  Summarize demultiplexed reads
   4  Quality filter (q-score ≥ 30)
   5  Tabulate filter stats
   6  Deblur denoising           ← restart after failed denoising
   7  Feature-table summary
   8  Deblur stats visualize
   9  Phylogenetic tree (MAFFT + FastTree)
  ── Branch A: original (steps 10–19) ────────────────────────────────────────
  10  Core-metrics-phylogenetic  ← RAREFACTION at depth ${DEPTH}
  11  Taxonomy classification
  12  Taxonomy visualize
  13  Taxa barplot
  14  PERMANOVA (all distance matrices)
  15  Alpha diversity significance
  16  Emperor PCoA plots
  ***
  17  Summary report ← SKIPED- WILL DO IN R
  18  R: ancombc.R ← SKIPED- WILL DO IN R
  19  R: plot_manager.R ← SKIPED- WILL DO IN R
  ***
  ── Branch B: raspergade (steps 20–27) ──────────────────────────────────────
  20  Taxonomy classification    ← shared artifact; skipped if step 11 ran
  21  Export Deblur table, taxonomy, and rep-seqs FASTA for R
  22  RasperGade16S copy-number correction (R; uses rep-seqs FASTA)
  23  Reimport corrected table into QIIME 2
  24  Rarefaction depth decision (R) ← prints recommendation before rarefying
  25  Core-metrics on corrected table ← RAREFACTION (depth may differ from A)
  26  Export corrected artifacts for downstream R
  27  Build phyloseq object — raspergade branch (R)
  ── Branch C: decontam (steps 28–33) ────────────────────────────────────────
  28  Export Deblur table to BIOM/TSV for R
  29  decontam frequency correction (R; requires -q QUANT_COL)
  30  Reimport decontam-filtered table into QIIME 2
  31  Core-metrics on decontam table ← RAREFACTION
  32  Export decontam artifacts for downstream R
  33  Build phyloseq object — decontam branch (R)
  ── Final ────────────────────────────────────────────────────────────────────
  34  Cross-branch QC table and phyloseq list (R)

Examples:
  # Full run, all branches (decontam skipped if -q not given):
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -b all

  # Original branch only, full pipeline:
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -b original

  # Resume from deblur after crash (steps 1–5 artifacts must exist):
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -S 6

  # Jump directly to RasperGade16S R correction (steps 1–21 must already exist):
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -b raspergade -S 22

  # Jump to rarefaction review for raspergade branch:
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -b raspergade -S 24

  # Force redo tree + everything downstream:
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -F -S 9

  # decontam branch only (requires DNA concentration column in metadata):
  $0 -f data/ -s map.tsv -m meta.tsv -o out/ -c clf/ -b decontam -q quant_reading
EOF
    exit 1
}

# ------------------------------
# Argument parsing
# ------------------------------
while getopts ":f:s:m:o:c:d:t:g:b:q:S:Fh" opt; do
    case "$opt" in
        f) FASTQ_DIR="$OPTARG"    ;;
        s) SAMPLE_MAP="$OPTARG"   ;;
        m) METADATA="$OPTARG"     ;;
        o) OUTPUT_DIR="$OPTARG"   ;;
        c) CLASSIFIER_DIR="$OPTARG" ;;
        d) DEPTH="$OPTARG"        ;;
        t) THREADS="$OPTARG"      ;;
        g) GROUP_COL="$OPTARG"    ;;
        b) BRANCHES="$OPTARG"     ;;
        q) QUANT_COL="$OPTARG"    ;;
        S) START_STEP="$OPTARG"   ;;
        F) FORCE=true             ;;
        h) usage                  ;;
        :) err "Missing argument for -$OPTARG" ;;
        \?) err "Unknown option: -$OPTARG" ;;
    esac
done

# Validate required args
[[ -z "$FASTQ_DIR" || -z "$SAMPLE_MAP" || -z "$METADATA" || -z "$OUTPUT_DIR" || -z "$CLASSIFIER_DIR" ]] && usage

# Validate branch value
case "$BRANCHES" in
    original|raspergade|decontam|all) ;;
    *) err "Invalid -b value '${BRANCHES}'. Choose: original | raspergade | decontam | all" ;;
esac

# Decontam requires quant column
if [[ "$BRANCHES" == "decontam" && -z "$QUANT_COL" ]]; then
    err "Branch 'decontam' requires -q QUANT_COL (metadata column with DNA concentration values)."
fi
if [[ "$BRANCHES" == "all" && -z "$QUANT_COL" ]]; then
    warn "No -q QUANT_COL provided; decontam branch will be skipped even with -b all."
fi

mkdir -p "$OUTPUT_DIR"
FASTQ_DIR=$(realpath "$FASTQ_DIR")
SAMPLE_MAP=$(realpath "$SAMPLE_MAP")
METADATA=$(realpath "$METADATA")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
CLASSIFIER_DIR=$(realpath "$CLASSIFIER_DIR")

# Path shortcuts (original dirs preserved exactly)
ART="$OUTPUT_DIR/artifacts"
TAX="$OUTPUT_DIR/taxonomy_classification"
STAT="$OUTPUT_DIR/statistical_tests"
RASP="$OUTPUT_DIR/branch_raspergade"
DECO="$OUTPUT_DIR/branch_decontam"
ROBJ="$OUTPUT_DIR/r_objects"
SCPTS="$OUTPUT_DIR/scripts"

# Branch B representative-sequence FASTA exported from QIIME 2 rep-seqs.qza.
# RasperGade16S predicts GCN from sequences, so this file is required for step 22.
RASP_REP_SEQS_FASTA="$ROBJ/raspergade/rep-seqs/dna-sequences.fasta"

mkdir -p "$ART" "$TAX" "$STAT" "$RASP" "$DECO" \
         "$ROBJ/original" "$ROBJ/raspergade" "$ROBJ/decontam" "$SCPTS"

# Find classifier
CLASSIFIER_FILE=$(find "$CLASSIFIER_DIR" -maxdepth 1 -type f -name "*.qza" | sort | head -n1 || true)
[[ -z "$CLASSIFIER_FILE" ]] && err "No classifier .qza found in $CLASSIFIER_DIR"

# ------------------------------
# Header
# ------------------------------
cat <<EOF
IMPORTANT: Default parameters are optimized for microbiome datasets processed at Karolinska Institutet.
===== QIIME2 Microbiome Revision Workflow =====
Input FASTQ directory:   $FASTQ_DIR
Sample map:              $SAMPLE_MAP
Metadata:                $METADATA
Output directory:        $OUTPUT_DIR
Classifier:              $CLASSIFIER_FILE
Branches:                $BRANCHES
Sampling depth:          $DEPTH
Threads:                 $THREADS
PERMANOVA column:        $GROUP_COL
Start step:              $START_STEP
Force rerun:             $FORCE
DNA conc. column:        ${QUANT_COL:-(not set; decontam branch inactive)}
===============================================
EOF

# ============================================================
# CHECKPOINT SYSTEM
# ============================================================
# checkpoint STEP_NUM ARTIFACT
#   Advances global step counter to STEP_NUM.
#   Returns 0 = proceed with step / 1 = skip this step.
#   Before START_STEP: validates artifact exists, then skips.
#   After  START_STEP: skips only if artifact exists AND FORCE=false.
# ============================================================

_STEP=0

checkpoint() {
    local n="$1"
    local artifact="${2:-}"
    _STEP=$n

    if [[ $_STEP -lt $START_STEP ]]; then
        # We are before the requested start: validate prerequisite then skip
        if [[ -n "$artifact" && ! -e "$artifact" ]]; then
            err "Step $_STEP: required artifact missing.\n  Path: $artifact\n  Re-run from step $_STEP (use -S $_STEP or earlier)."
        fi
        [[ -n "$artifact" ]] && skipped "$artifact"
        return 1   # skip
    fi

    # At or after start: skip only when artifact exists and not forcing
    if [[ "$FORCE" == "false" && -n "$artifact" && -e "$artifact" ]]; then
        skipped "$artifact"
        return 1   # skip
    fi

    return 0   # proceed
}

# branch_active NAME — returns 0 if this branch should run, 1 otherwise
branch_active() {
    local b="$1"
    if [[ "$BRANCHES" == "$b" ]]; then return 0; fi
    if [[ "$BRANCHES" == "all" ]]; then
        # decontam needs quant col even under 'all'
        if [[ "$b" == "decontam" && -z "$QUANT_COL" ]]; then return 1; fi
        return 0
    fi
    return 1
}

# RasperGade16S calls external HMMER/Easel binaries from R.
# The QIIME2 environment may contain hmmalign but still lack Easel miniapps
# such as esl-reformat, so check them before launching the expensive R step.
check_raspergade_cli_deps() {
    local missing=0
    for exe in hmmalign esl-reformat; do
        if ! command -v "$exe" >/dev/null 2>&1; then
            warn "Missing required executable for RasperGade16S: $exe"
            missing=1
        else
            info "Found $exe: $(command -v "$exe")"
        fi
    done
    if [[ "$missing" -ne 0 ]]; then
        err "RasperGade16S dependency check failed. Install Easel/HMMER tools in the active environment, for example: mamba install -c conda-forge -c bioconda easel hmmer"
    fi
}

# RasperGade16S loads some R packages dynamically during sequence alignment
# and parsing. In particular, seqinr may not be installed automatically in
# some QIIME2/R environments, so catch this before the long R step starts.
check_raspergade_r_deps() {
    local missing
    missing=$(Rscript -e 'pkgs <- c("RasperGade16S", "seqinr", "biomformat", "tidyverse"); miss <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; cat(paste(miss, collapse=" "))' 2>/dev/null || true)
    if [[ -n "$missing" ]]; then
        warn "Missing required R package(s) for RasperGade16S: $missing"
        err "Install missing R dependencies in the active environment. For conda/mamba, try: mamba install -c conda-forge -c bioconda r-seqinr r-biomformat r-tidyverse"
    fi
}

# run_r SCRIPT_PATH — execute an R script; pass env vars; tee log alongside script
run_r() {
    local script="$1"
    step "R" "$script"
    OUTPUT_DIR="$OUTPUT_DIR" METADATA="$METADATA" DEPTH="$DEPTH" \
    QUANT_COL="$QUANT_COL" THREADS="$THREADS" GROUP_COL="$GROUP_COL" \
    BRANCHES="$BRANCHES" REP_SEQS="$RASP_REP_SEQS_FASTA" \
    NORMALIZE_GCN="${NORMALIZE_GCN:-FALSE}" \
        Rscript "$script" 2>&1 | tee "${script%.R}.log"
}

# ============================================================
# MANIFEST CREATION (same logic as original workflow)
# ============================================================
MANIFEST="$OUTPUT_DIR/manifest_v2.csv"

build_manifest() {
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"

    local first_line
    first_line=$(sed -n '1p' "$SAMPLE_MAP" || true)
    local DELIM=$'\t'
    [[ "$first_line" == *","* && "$first_line" != *$'\t'* ]] && DELIM=","

    while IFS="$DELIM" read -r col1 col2 _rest || [[ -n "$col1" ]]; do
        [[ "$col1" =~ ^(#|sample|Sample|sample-id) ]] && continue

        local sample_prefix sample_id
        sample_prefix=$(echo "$col1" | xargs)
        sample_id=$(echo "${col2:-$col1}" | xargs)
        [[ -z "$sample_prefix" || -z "$sample_id" ]] && { warn "Skipping malformed line: $col1"; continue; }

        local f1 f2
        f1=$(find "$FASTQ_DIR" -type f -name "${sample_prefix}_1.fq.gz" | head -n1 || true)
        f2=$(find "$FASTQ_DIR" -type f -name "${sample_prefix}_2.fq.gz" | head -n1 || true)

        if [[ -n "$f1" && -n "$f2" ]]; then
            printf "%s\t%s\t%s\n" "$sample_id" "$f1" "$f2" >> "$MANIFEST"
        else
            warn "Missing FASTQ files for prefix '${sample_prefix}' — skipped."
        fi
    done < <(tail -n +2 "$SAMPLE_MAP" | sed '/^[[:space:]]*$/d')

    info "Manifest written to $MANIFEST"
}

# ============================================================
# ── SHARED STEPS 1–9 ─────────────────────────────────────────
# ============================================================

# 1. Import
step 1 "Build manifest and import paired-end reads"
if checkpoint 1 "$ART/imported_reads_demux.qza"; then
    build_manifest
    qrun tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST" \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "$ART/imported_reads_demux.qza"
fi

# 2. Merge pairs
step 2 "Merging paired-end reads"
if checkpoint 2 "$ART/demux-merged.qza"; then
    qrun vsearch merge-pairs \
        --i-demultiplexed-seqs "$ART/imported_reads_demux.qza" \
        --o-merged-sequences   "$ART/demux-merged.qza" \
        --o-unmerged-sequences "$ART/dmux-unmerged.qza" \
        --p-threads "$THREADS"
fi

# 3. Demux summarize
step 3 "Summarizing demultiplexed data"
if checkpoint 3 "$ART/demux-merged.qzv"; then
    qrun demux summarize \
        --i-data "$ART/demux-merged.qza" \
        --o-visualization "$ART/demux-merged.qzv"
fi

# 4. Quality filter
step 4 "Quality filtering (min-quality 30)"
if checkpoint 4 "$ART/demux-merged.filtered.qza"; then
    qrun quality-filter q-score \
        --i-demux "$ART/demux-merged.qza" \
        --o-filtered-sequences "$ART/demux-merged.filtered.qza" \
        --o-filter-stats       "$ART/demux-merged.filter_stats.qza" \
        --p-min-quality 30
fi

# 5. Tabulate filter stats
step 5 "Tabulating filter stats"
if checkpoint 5 "$ART/demux-merged.filter_stats.qzv"; then
    qrun metadata tabulate \
        --m-input-file "$ART/demux-merged.filter_stats.qza" \
        --o-visualization "$ART/demux-merged.filter_stats.qzv"
fi

# 6. Deblur denoising
step 6 "Deblur denoising (trim-length 400)"
if checkpoint 6 "$ART/table.qza"; then
    qrun deblur denoise-16S \
        --i-demultiplexed-seqs "$ART/demux-merged.filtered.qza" \
        --p-trim-length 400 \
        --p-sample-stats \
        --o-representative-sequences "$ART/rep-seqs.qza" \
        --o-table                    "$ART/table.qza" \
        --o-stats                    "$ART/deblur-stats.qza"
fi

# 7. Feature table summary
step 7 "Summarizing feature table"
if checkpoint 7 "$ART/table.qzv"; then
    qrun feature-table summarize \
        --i-table "$ART/table.qza" \
        --o-visualization "$ART/table.qzv" \
        --m-sample-metadata-file "$METADATA"
fi

# 8. Deblur stats visualize
step 8 "Visualizing deblur stats"
if checkpoint 8 "$ART/deblur-stats.qzv"; then
    qrun deblur visualize-stats \
        --i-deblur-stats "$ART/deblur-stats.qza" \
        --o-visualization "$ART/deblur-stats.qzv"
fi

# 9. Phylogenetic tree
step 9 "Constructing phylogenetic tree (MAFFT + FastTree)"
if checkpoint 9 "$ART/rooted-tree.qza"; then
    qrun phylogeny align-to-tree-mafft-fasttree \
        --i-sequences      "$ART/rep-seqs.qza" \
        --o-alignment      "$ART/aligned-rep-seqs.qza" \
        --o-masked-alignment "$ART/masked-aligned-rep-seqs.qza" \
        --o-tree           "$ART/unrooted-tree.qza" \
        --o-rooted-tree    "$ART/rooted-tree.qza" \
        --p-n-threads "$THREADS"
fi

# ============================================================
# ── BRANCH A: ORIGINAL (steps 10–19) ─────────────────────────
# ============================================================

if branch_active original; then
    bhead "A — Original pipeline (steps 10–19)"

    # 10. Core-metrics-phylogenetic  ← RAREFACTION HERE
    # NOTE: --output-dir is only accepted when directory does not yet exist.
    #       The checkpoint uses a sentinel file inside the directory instead.
    step 10 "[A] Core-metrics-phylogenetic (RAREFACTION at depth $DEPTH)"
    if checkpoint 10 "$ART/diversities/weighted_unifrac_distance_matrix.qza"; then
        # Remove partial directory if it exists and we are forcing
        [[ "$FORCE" == "true" && -d "$ART/diversities" ]] && rm -rf "$ART/diversities"
        qrun diversity core-metrics-phylogenetic \
            --i-phylogeny "$ART/rooted-tree.qza" \
            --i-table     "$ART/table.qza" \
            --p-sampling-depth "$DEPTH" \
            --m-metadata-file "$METADATA" \
            --output-dir "$ART/diversities"
        # ── Rarefaction checkpoint: report samples dropped ──────────────
        info "Rarefaction depth: $DEPTH"
        python3 - "$ART/table.qza" "$DEPTH" <<'PY'
import sys, zipfile, json, io
qza_path, depth = sys.argv[1], int(sys.argv[2])
with zipfile.ZipFile(qza_path) as z:
    biom_name = [n for n in z.namelist() if n.endswith("feature-table.biom")]
    if not biom_name:
        print("  [rarefaction check] Could not locate feature-table.biom inside QZA.")
        sys.exit(0)
    import biom
    with z.open(biom_name[0]) as f:
        tbl = biom.parse_table(io.TextIOWrapper(f))
depths = {sid: int(tbl.data(sid, axis="sample").sum()) for sid in tbl.ids(axis="sample")}
dropped = [(s,d) for s,d in depths.items() if d < depth]
print(f"  Samples at depth >= {depth}: {len(depths)-len(dropped)}/{len(depths)}")
if dropped:
    print(f"  DROPPED at depth {depth}:")
    for s,d in sorted(dropped, key=lambda x: x[1]):
        print(f"    {s}: {d} reads")
else:
    print("  No samples dropped.")
PY
    fi

    # 11. Taxonomy classification
    step 11 "[A] Assigning taxonomy"
    if checkpoint 11 "$TAX/taxonomy.qza"; then
        qrun feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_FILE" \
            --i-reads      "$ART/rep-seqs.qza" \
            --o-classification "$TAX/taxonomy.qza" \
            --p-n-jobs "$THREADS"
    fi

    # 12. Taxonomy visualize
    step 12 "[A] Visualizing taxonomy"
    if checkpoint 12 "$TAX/taxonomy.qzv"; then
        qrun metadata tabulate \
            --m-input-file "$TAX/taxonomy.qza" \
            --o-visualization "$TAX/taxonomy.qzv"
    fi

    # 13. Taxa barplot
    step 13 "[A] Creating taxonomy barplot"
    if checkpoint 13 "$TAX/taxa-bar-plots.qzv"; then
        qrun taxa barplot \
            --i-table    "$ART/table.qza" \
            --i-taxonomy "$TAX/taxonomy.qza" \
            --m-metadata-file "$METADATA" \
            --o-visualization "$TAX/taxa-bar-plots.qzv"
    fi

    # 14. PERMANOVA tests
    step 14 "[A] Running PERMANOVA tests"
    DIST_NAMES=(unweighted_unifrac_distance_matrix weighted_unifrac_distance_matrix \
                jaccard_distance_matrix bray_curtis_distance_matrix)
    for dbase in "${DIST_NAMES[@]}"; do
        dm_qza="$ART/diversities/${dbase}.qza"
        out_viz="$STAT/permanova_${dbase}.qzv"
        if checkpoint 14 "$out_viz" && [[ -f "$dm_qza" ]]; then
            qrun diversity beta-group-significance \
                --i-distance-matrix "$dm_qza" \
                --m-metadata-file "$METADATA" \
                --m-metadata-column "$GROUP_COL" \
                --o-visualization "$out_viz" \
                --p-pairwise
        elif [[ ! -f "$dm_qza" ]]; then
            warn "Distance matrix not found: $dm_qza (skipping)"
        fi
    done

    # 15. Alpha diversity significance
    step 15 "[A] Running alpha diversity significance tests"
    ALPHA_METRICS=(faith_pd observed_features shannon evenness)
    for metric in "${ALPHA_METRICS[@]}"; do
        vec_qza="$ART/diversities/${metric}_vector.qza"
        out_qzv="$STAT/alpha_significance_${metric}.qzv"
        if checkpoint 15 "$out_qzv" && [[ -f "$vec_qza" ]]; then
            qrun diversity alpha-group-significance \
                --i-alpha-diversity "$vec_qza" \
                --m-metadata-file "$METADATA" \
                --o-visualization "$out_qzv"
        elif [[ ! -f "$vec_qza" ]]; then
            warn "Alpha vector not found: $vec_qza (skipping)"
        fi
    done

    # 16. Emperor PCoA plots
    step 16 "[A] Generating Emperor PCoA plots"
    for dbase in "${DIST_NAMES[@]}"; do
        pcoa_base=${dbase%_distance_matrix}
        pcoa_qza="$ART/diversities/${pcoa_base}_pcoa_results.qza"
        out_qzv="$ART/diversities/${pcoa_base}_emperor.qzv"
        if checkpoint 16 "$out_qzv" && [[ -f "$pcoa_qza" ]]; then
            qrun emperor plot \
                --i-pcoa "$pcoa_qza" \
                --m-metadata-file "$METADATA" \
                --o-visualization "$out_qzv"
        elif [[ ! -f "$pcoa_qza" ]]; then
            warn "PCoA results not found for ${pcoa_base} (skipping)"
        fi
    done

    # 17. Summary report
    # step 17 "[A] Generating summary report"
    # if checkpoint 17 "$STAT/summary_report.txt"; then
    #     {   echo "QIIME2 Statistical Tests Summary"
    #         echo "================================"
    #         echo "Analysis performed on: $(date --iso-8601=seconds)"
    #         echo "Metadata column used for grouping: $GROUP_COL"
    #         echo "Sampling depth: $DEPTH"
    #         echo
    #         echo "PERMANOVA Tests:"
    #         printf -- "- %s\n" "Unweighted UniFrac" "Weighted UniFrac" "Jaccard" "Bray-Curtis"
    #         echo
    #         echo "Alpha Diversity Tests:"
    #         printf -- "- %s\n" "faith_pd" "observed_features" "shannon" "evenness"
    #         echo
    #         echo "View .qzv files at: https://view.qiime2.org/"
    #     } > "$STAT/summary_report.txt"
    # fi

    # info "===== Branch A complete ====="
    # info "Artifacts in: $ART"
    # info "Taxonomy in:  $TAX"
    # info "Stats in:     $STAT"

    # # 18. R: ancombc
    # step 18 "[A] Running R downstream analysis: ancombc.R"
    # if checkpoint 18 "$OUTPUT_DIR/ancombc_done.flag"; then
    #     Rscript "$SCRIPT_DIR/ancombc.R" "$OUTPUT_DIR" "$METADATA"
    #     touch "$OUTPUT_DIR/ancombc_done.flag"
    # fi

    # # 19. R: plot_manager
    # step 19 "[A] Running R downstream analysis: plot_manager.R"
    # if checkpoint 19 "$OUTPUT_DIR/plot_manager_done.flag"; then
    #     Rscript "$SCRIPT_DIR/plot_manager.R" "$OUTPUT_DIR" "$METADATA"
    #     touch "$OUTPUT_DIR/plot_manager_done.flag"
    # fi
fi   # end Branch A

# ============================================================
# ── BRANCH B: RasperGade16S (steps 20–27) ────────────────────
# ============================================================

if branch_active raspergade; then
    bhead "B — RasperGade16S copy-number correction (steps 20–27)"

    # 20. Taxonomy — reuse Branch A artifact if already present
    step 20 "[B] Taxonomy classification (reused from Branch A if present)"
    if checkpoint 20 "$TAX/taxonomy.qza"; then
        qrun feature-classifier classify-sklearn \
            --i-classifier "$CLASSIFIER_FILE" \
            --i-reads      "$ART/rep-seqs.qza" \
            --o-classification "$TAX/taxonomy.qza" \
            --p-n-jobs "$THREADS"
    fi

    # 21. Export Deblur table, taxonomy, and representative sequences for R
    step 21 "[B] Exporting Deblur table, taxonomy, and rep-seqs FASTA"
    if checkpoint 21 "$ROBJ/raspergade/step21_export_complete.flag"; then
        mkdir -p "$ROBJ/raspergade" "$ROBJ/raspergade/taxonomy" "$ROBJ/raspergade/rep-seqs"

        qrun tools export \
            --input-path "$ART/table.qza" \
            --output-path "$ROBJ/raspergade/"
        biom convert \
            -i "$ROBJ/raspergade/feature-table.biom" \
            -o "$ROBJ/raspergade/feature-table.tsv" \
            --to-tsv

        qrun tools export \
            --input-path "$TAX/taxonomy.qza" \
            --output-path "$ROBJ/raspergade/taxonomy/"

        qrun tools export \
            --input-path "$ART/rep-seqs.qza" \
            --output-path "$ROBJ/raspergade/rep-seqs/"

        [[ -s "$ROBJ/raspergade/feature-table.biom" ]] || err "Missing exported BIOM table."
        [[ -s "$ROBJ/raspergade/taxonomy/taxonomy.tsv" ]] || err "Missing exported taxonomy.tsv."
        [[ -s "$RASP_REP_SEQS_FASTA" ]] || err "Missing exported representative-sequence FASTA: $RASP_REP_SEQS_FASTA"

        touch "$ROBJ/raspergade/step21_export_complete.flag"
        info "Exported BIOM:      $ROBJ/raspergade/feature-table.biom"
        info "Exported taxonomy:  $ROBJ/raspergade/taxonomy/taxonomy.tsv"
        info "Exported rep-seqs:  $RASP_REP_SEQS_FASTA"
    fi

    # 22. RasperGade16S R correction
    step 22 "[B] RasperGade16S copy-number correction (R)"
    RASP_R="$SCPTS/step22_raspergade.R"
    cat > "$RASP_R" <<'REOF'
suppressPackageStartupMessages({
    library(RasperGade16S)
    library(seqinr)
    library(biomformat)
    library(tidyverse)
})

# ─────────────────────────────────────────────────────────────
# Inputs
# ─────────────────────────────────────────────────────────────
base_dir <- Sys.getenv("OUTPUT_DIR", unset = ".")
out_dir  <- file.path(base_dir, "r_objects", "raspergade")

biom_path <- file.path(out_dir, "feature-table.biom")
tax_path  <- file.path(out_dir, "taxonomy", "taxonomy.tsv")
repseq_candidates <- c(
    Sys.getenv("REP_SEQS", unset = NA_character_),
    file.path(out_dir, "rep-seqs", "dna-sequences.fasta"),
    file.path(out_dir, "dna-sequences.fasta"),
    file.path(out_dir, "rep-seqs.fasta")
)
repseq_fasta <- repseq_candidates[!is.na(repseq_candidates) & file.exists(repseq_candidates)][1]

normalize_gcn <- tolower(Sys.getenv("NORMALIZE_GCN", unset = "FALSE")) %in%
    c("1", "true", "t", "yes", "y")
depth_thresh <- as.integer(Sys.getenv("DEPTH", unset = "74470"))

if (!file.exists(biom_path)) stop("Missing BIOM file: ", biom_path)
if (!file.exists(tax_path))  stop("Missing taxonomy file: ", tax_path)
if (is.na(repseq_fasta)) {
    stop(
        "No representative-sequence FASTA found. RasperGade16S needs ASV sequences to predict GCN.\n",
        "Expected one of:\n  ", paste(repseq_candidates[!is.na(repseq_candidates)], collapse = "\n  ")
    )
}

# ─────────────────────────────────────────────────────────────
# Read input table and taxonomy
# ─────────────────────────────────────────────────────────────
biom_in <- read_biom(biom_path)
otu_mat <- as(biom_data(biom_in), "matrix")   # ASVs × samples
storage.mode(otu_mat) <- "numeric"

tax_raw <- read_tsv(
    tax_path,
    col_names = c("feature_id", "taxon", "confidence"),
    skip = 1,
    show_col_types = FALSE
)

cat("── Input ──────────────────────────────────────────────────\n")
cat("BIOM:     ", biom_path, "\n", sep = "")
cat("Rep-seqs: ", repseq_fasta, "\n", sep = "")
cat("ASVs:", nrow(otu_mat), "  Samples:", ncol(otu_mat), "\n")
cat("Taxonomy rows:", nrow(tax_raw), "\n")
cat("Normalize after GCN correction:", normalize_gcn, "\n")

parse_silva <- function(x) {
    parts <- stringr::str_split(x, ";\\s*")[[1]]
    clean <- sub("^[dkpcofgs]__", "", parts)
    length(clean) <- 7
    setNames(as.list(clean), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
}

tax_df <- tax_raw |>
    rowwise() |>
    mutate(p = list(parse_silva(taxon))) |>
    unnest_wider(p) |>
    select(ASV = feature_id, Kingdom:Species) |>
    ungroup() |>
    as.data.frame()

shared_tax <- intersect(rownames(otu_mat), tax_df$ASV)
cat("ASVs matched to taxonomy:", length(shared_tax), "\n")
if (length(shared_tax) < nrow(otu_mat)) {
    warning(nrow(otu_mat) - length(shared_tax), " ASVs not in taxonomy; they will be dropped.")
}

otu_sub <- otu_mat[shared_tax, , drop = FALSE]
tax_sub <- tax_df[match(shared_tax, tax_df$ASV), , drop = FALSE]
rownames(tax_sub) <- tax_sub$ASV

# ─────────────────────────────────────────────────────────────
# Predict or load 16S rRNA gene copy number by RasperGade16S
# ─────────────────────────────────────────────────────────────
extract_gcn <- function(pred) {
    if (is.numeric(pred) && !is.null(names(pred))) return(pred)

    if (is.list(pred)) {
        for (nm in c("GCN", "gcn", "copy_number", "copyNumber")) {
            if (!is.null(pred[[nm]]) && is.numeric(pred[[nm]]) && !is.null(names(pred[[nm]]))) {
                return(pred[[nm]])
            }
        }
        if (!is.null(pred$tab)) {
            tab <- as.data.frame(pred$tab)
        } else {
            dfs <- pred[vapply(pred, is.data.frame, logical(1))]
            if (length(dfs) == 0) stop("Could not find a GCN table in the RasperGade16S prediction object.")
            tab <- as.data.frame(dfs[[1]])
        }
    } else if (is.data.frame(pred)) {
        tab <- pred
    } else {
        stop("Unsupported RasperGade16S prediction object class: ", paste(class(pred), collapse = ", "))
    }

    id_col  <- intersect(c("label", "ASV", "feature_id", "FeatureID", "id"), names(tab))[1]
    gcn_col <- intersect(c("GCN", "gcn", "x", "copy_number", "copyNumber"), names(tab))[1]
    if (is.na(id_col) || is.na(gcn_col)) {
        stop("Could not identify ID/GCN columns in prediction table. Columns are: ",
             paste(names(tab), collapse = ", "))
    }

    gcn <- as.numeric(tab[[gcn_col]])
    names(gcn) <- as.character(tab[[id_col]])
    gcn
}

gcn_rds <- Sys.getenv("GCN_RDS", unset = "")
gcn_tsv <- Sys.getenv("GCN_TSV", unset = "")

if (nzchar(gcn_rds) && file.exists(gcn_rds)) {
    cat("── Loading existing GCN RDS ───────────────────────────────\n")
    pred <- readRDS(gcn_rds)
    gcn <- extract_gcn(pred)
} else if (nzchar(gcn_tsv) && file.exists(gcn_tsv)) {
    cat("── Loading existing GCN TSV ───────────────────────────────\n")
    gcn_tbl <- read_tsv(gcn_tsv, show_col_types = FALSE)
    id_col  <- intersect(c("ASV", "feature_id", "FeatureID", "label", "id"), names(gcn_tbl))[1]
    gcn_col <- intersect(c("GCN", "gcn", "x", "copy_number", "copyNumber"), names(gcn_tbl))[1]
    if (is.na(id_col) || is.na(gcn_col)) {
        stop("GCN_TSV must contain an ID column and a GCN column. Found columns: ",
             paste(names(gcn_tbl), collapse = ", "))
    }
    gcn <- as.numeric(gcn_tbl[[gcn_col]])
    names(gcn) <- as.character(gcn_tbl[[id_col]])
} else {
    cat("── Predicting 16S GCN with RasperGade16S ──────────────────\n")
    cli_needed <- c("hmmalign", "esl-reformat")
    cli_missing <- cli_needed[!nzchar(Sys.which(cli_needed))]
    if (length(cli_missing) > 0) {
        stop(
            "RasperGade16S requires external HMMER/Easel command-line tools, but these are missing from PATH: ",
            paste(cli_missing, collapse = ", "),
            "\nInstall them in the active conda environment, e.g.: ",
            "mamba install -c conda-forge -c bioconda easel hmmer",
            "\nThen check: which hmmalign && which esl-reformat"
        )
    }
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(out_dir)
    pred <- RasperGade16S::predict_16SGCN_from_sequences(seqs = normalizePath(repseq_fasta))
    saveRDS(pred, file.path(out_dir, "raspergade16s-prediction.rds"))
    gcn <- extract_gcn(pred)
}

gcn <- gcn[!is.na(names(gcn)) & !is.na(gcn) & gcn > 0]
matched <- intersect(rownames(otu_sub), names(gcn))
cat("ASVs with predicted GCN:", length(matched), "\n")

missing_gcn <- setdiff(rownames(otu_sub), names(gcn))
if (length(missing_gcn) > 0) {
    warning(length(missing_gcn), " ASVs lack GCN predictions and will be dropped.")
    write_tsv(
        tibble(ASV = missing_gcn) |>
            left_join(as_tibble(tax_sub), by = "ASV"),
        file.path(out_dir, "asvs-missing-raspergade-gcn.tsv")
    )
}

if (length(matched) == 0) {
    stop(
        "No overlap between BIOM feature IDs and RasperGade16S FASTA/prediction IDs.\n",
        "Check that FASTA headers from rep-seqs.qza are the same ASV IDs as BIOM row names."
    )
}

otu_sub <- otu_sub[matched, , drop = FALSE]
tax_sub <- tax_sub[matched, , drop = FALSE]
gcn <- gcn[matched]

write_tsv(
    tibble(ASV = names(gcn), GCN = as.numeric(gcn)) |>
        left_join(as_tibble(tax_sub), by = "ASV"),
    file.path(out_dir, "raspergade16s-gcn.tsv")
)

cat("GCN summary:\n")
print(summary(as.numeric(gcn)))

# ─────────────────────────────────────────────────────────────
# Copy-number correction
# ─────────────────────────────────────────────────────────────
# RasperGade16S provides the predicted ASV-level GCN. The abundance correction
# is count / predicted_GCN for each ASV. If NORMALIZE_GCN=TRUE, each sample is
# rescaled back to its original library size after correction.
cat("── 16S GCN abundance correction ───────────────────────────\n")
corrected <- sweep(otu_sub, 1, gcn[rownames(otu_sub)], "/")

if (normalize_gcn) {
    original_depths <- colSums(otu_sub)
    corrected_depths <- colSums(corrected)
    scale_factors <- original_depths / corrected_depths
    scale_factors[!is.finite(scale_factors)] <- 1
    corrected <- sweep(corrected, 2, scale_factors, "*")
}

corrected <- corrected[rowSums(corrected) > 0, , drop = FALSE]
storage.mode(corrected) <- "numeric"

cat("Samples: before =", ncol(otu_sub), " after =", ncol(corrected), "\n")
cat("ASVs:    before =", nrow(otu_sub), " after =", nrow(corrected), "\n")
cat("Max corrected value:", max(corrected, na.rm = TRUE), "\n")

round_preserve_depth <- function(mat, target_depths = NULL) {
    mat <- as.matrix(mat)
    if (is.null(target_depths)) target_depths <- round(colSums(mat))
    target_depths <- as.integer(round(target_depths[colnames(mat)]))

    out <- sapply(seq_len(ncol(mat)), function(j) {
        x <- mat[, j]
        flo <- floor(x)
        target <- target_depths[j]
        add_n <- target - sum(flo)
        if (add_n > 0) {
            idx <- order(x - flo, decreasing = TRUE)[seq_len(min(add_n, length(x)))]
            flo[idx] <- flo[idx] + 1L
        } else if (add_n < 0) {
            idx <- which(flo > 0)
            if (length(idx) > 0) {
                idx <- idx[order(x[idx] - flo[idx], decreasing = FALSE)]
                idx <- idx[seq_len(min(abs(add_n), length(idx)))]
                flo[idx] <- flo[idx] - 1L
            }
        }
        flo
    })
    rownames(out) <- rownames(mat)
    colnames(out) <- colnames(mat)
    storage.mode(out) <- "integer"
    out
}

integer_depth_target <- round(colSums(corrected))
corrected_integer <- round_preserve_depth(corrected, integer_depth_target)
corrected_integer <- corrected_integer[rowSums(corrected_integer) > 0, , drop = FALSE]

cat("\n── Per-sample library sizes after correction ──────────────\n")
cat("Fractional corrected table:\n")
print(summary(as.numeric(colSums(corrected))))
cat("Integer corrected table:\n")
print(summary(as.numeric(colSums(corrected_integer))))

at_risk <- sum(colSums(corrected_integer) < depth_thresh)
cat(sprintf("Samples below rarefaction depth (%d): %d\n", depth_thresh, at_risk))
if (at_risk > 0) {
    cat("At-risk samples:\n")
    print(sort(colSums(corrected_integer)[colSums(corrected_integer) < depth_thresh]))
    cat("Recommendation: lower rarefaction depth for the GCN-corrected branch, or do CLR/TSS analyses in R.\n")
} else {
    cat("All samples above", depth_thresh, "— this rarefaction depth is safe for the integer corrected table.\n")
}

# ─────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────
out_tsv_fractional <- file.path(out_dir, "corrected-table-fractional.tsv")
out_tsv_integer    <- file.path(out_dir, "corrected-table.tsv")
out_rds_fractional <- file.path(out_dir, "corrected-table-fractional.rds")
out_rds_integer    <- file.path(out_dir, "corrected-table.rds")

write_tsv(as.data.frame(corrected) |> rownames_to_column("#OTU ID"), out_tsv_fractional)
write_tsv(as.data.frame(corrected_integer) |> rownames_to_column("#OTU ID"), out_tsv_integer)
saveRDS(corrected, out_rds_fractional)
saveRDS(corrected_integer, out_rds_integer)

cat("\nOutputs:\n")
cat("  Fractional corrected TSV: ", out_tsv_fractional, "\n", sep = "")
cat("  Integer corrected TSV:    ", out_tsv_integer, "\n", sep = "")
cat("  Fractional corrected RDS: ", out_rds_fractional, "\n", sep = "")
cat("  Integer corrected RDS:    ", out_rds_integer, "\n", sep = "")
cat("  ASV-level GCN table:      ", file.path(out_dir, "raspergade16s-gcn.tsv"), "\n", sep = "")
REOF
    if checkpoint 22 "$ROBJ/raspergade/corrected-table.rds"; then
        check_raspergade_cli_deps
        check_raspergade_r_deps
        run_r "$RASP_R"
    fi

    # 23. Reimport corrected table into QIIME 2
    step 23 "[B] Reimporting integer GCN-corrected table into QIIME 2"
    if checkpoint 23 "$RASP/table-corrected.qza"; then
        # Step 23a: Step 22 already writes an integer pseudo-count table for QIIME 2.
        [[ -s "$ROBJ/raspergade/corrected-table.tsv" ]] || err "Missing corrected integer TSV from step 22."
        [[ -s "$ROBJ/raspergade/corrected-table-fractional.tsv" ]] || err "Missing corrected fractional TSV from step 22."
        
        # Step 23b: Convert TSV to BIOM
        info "Converting TSV to BIOM..."
        biom convert \
            -i "$ROBJ/raspergade/corrected-table.tsv" \
            -o "$ROBJ/raspergade/corrected-table.biom" \
            --table-type="OTU table" --to-hdf5
        
        biom summarize-table \
            -i "$ROBJ/raspergade/corrected-table.biom" \
            | tee "$ROBJ/raspergade/corrected-table-summary.txt"
        
        # Step 23c: Import into QIIME 2
        info "Importing to QIIME 2..."
        qrun tools import \
            --type 'FeatureTable[Frequency]' \
            --input-path "$ROBJ/raspergade/corrected-table.biom" \
            --input-format BIOMV210Format \
            --output-path "$RASP/table-corrected.qza"
        
        qrun feature-table summarize \
            --i-table "$RASP/table-corrected.qza" \
            --m-sample-metadata-file "$METADATA" \
            --o-visualization "$RASP/table-corrected-summary.qzv"
        
        info "Corrected table is now rounded to integers and ready for diversity analysis."
    fi

    export OMP_NUM_THREADS="$THREADS"
    export UNIFRAC_USE_GPU=N

    # 24. Rarefaction depth decision (R) — always prints; no checkpoint on output file
    step 24 "[B] Rarefaction depth decision — comparing original vs corrected library sizes"
    RAREF_R="$SCPTS/step24_rarefaction_decision.R"
    cat > "$RAREF_R" <<'REOF'
suppressPackageStartupMessages({ library(biomformat); library(tidyverse) })

orig_biom  <- read_biom(file.path(Sys.getenv("OUTPUT_DIR"), "r_objects/raspergade/feature-table.biom"))
corr_biom  <- read_biom(file.path(Sys.getenv("OUTPUT_DIR"), "r_objects/raspergade/corrected-table.biom"))
orig_d <- colSums(as(biom_data(orig_biom), "matrix"))
corr_d <- colSums(as(biom_data(corr_biom), "matrix"))
thresh <- as.integer(Sys.getenv("DEPTH", "74470"))

common  <- intersect(names(orig_d), names(corr_d))
tbl <- tibble(sample = common,
              original  = as.integer(orig_d[common]),
              corrected = as.integer(corr_d[common]),
              drop_orig = original  < thresh,
              drop_corr = corrected < thresh)

cat(sprintf("Original rarefaction depth (workflow.sh): %d\n", thresh))
cat(sprintf("Samples dropped in Branch A:  %d / %d\n", sum(tbl$drop_orig), nrow(tbl)))
cat(sprintf("Samples dropped in Branch B:  %d / %d\n", sum(tbl$drop_corr), nrow(tbl)))

new_drops <- filter(tbl, !drop_orig & drop_corr)
if (nrow(new_drops) > 0) {
    cat("\nSamples newly at risk in Branch B:\n")
    print(select(new_drops, sample, original, corrected))
    safe_depth <- floor(min(tbl$corrected[!tbl$drop_orig]))
    cat(sprintf("\nRECOMMENDATION: use -d %d for Branch B (deepest safe depth).\n", safe_depth))
    cat("Or skip rarefaction and normalise with CLR/TSS in R.\n")
} else {
    cat(sprintf("\nAll samples safe at depth %d — no change needed for Branch B.\n", thresh))
}

ggplot(pivot_longer(tbl, c(original,corrected), names_to="branch", values_to="depth"),
       aes(x=branch, y=depth, group=sample)) +
    geom_line(alpha=.4, colour="grey40") +
    geom_point(aes(colour=depth < thresh), size=3) +
    geom_hline(yintercept=thresh, linetype="dashed", colour="red") +
    scale_colour_manual(values=c("FALSE"="steelblue","TRUE"="tomato"),
                        name=paste("Below depth =", thresh)) +
    scale_y_continuous(labels=scales::comma) +
    labs(title="Library sizes: original vs RasperGade16S corrected",
         y="Library size (reads)", x="Branch") +
    theme_bw()
ggsave(file.path(Sys.getenv("OUTPUT_DIR"), "scripts/step24_rarefaction_decision.pdf"),
       width=7, height=5)
cat("Plot saved to scripts/step24_rarefaction_decision.pdf\n")
REOF
    # This step always runs if we are at or past step 24 (no expensive output to checkpoint)
    if [[ $_STEP -lt $START_STEP ]]; then
        checkpoint 24 ""   # just advance counter; no artifact check needed
    else
        checkpoint 24 "" || true   # advance; run unconditionally when reached
        run_r "$RAREF_R"
    fi

    # 25. Core-metrics on corrected table  ← RAREFACTION
    # NOTE: inspect step 24 output before setting depth here.
    #       Override with -d NEW_DEPTH if step 24 recommended a lower value.
    step 25 "[B] Core-metrics-phylogenetic on corrected table (RAREFACTION at depth $DEPTH)"
    if checkpoint 25 "$RASP/diversities/weighted_unifrac_distance_matrix.qza"; then
        [[ "$FORCE" == "true" && -d "$RASP/diversities" ]] && rm -rf "$RASP/diversities"
        qrun diversity core-metrics-phylogenetic \
            --i-phylogeny "$ART/rooted-tree.qza" \
            --i-table     "$RASP/table-corrected.qza" \
            --p-sampling-depth "$DEPTH" \
            --m-metadata-file "$METADATA" \
            --output-dir "$RASP/diversities"
    fi

    # 26. Export final Branch B artifacts for downstream R comparison
    step 26 "[B] Exporting final corrected artifacts for R"
    if checkpoint 26 "$ROBJ/raspergade/final/feature-table.biom"; then
        mkdir -p "$ROBJ/raspergade/final"
        qrun tools export \
            --input-path "$RASP/table-corrected.qza" \
            --output-path "$ROBJ/raspergade/final/"
        biom convert \
            -i "$ROBJ/raspergade/final/feature-table.biom" \
            -o "$ROBJ/raspergade/final/feature-table.tsv" \
            --to-tsv
        qrun tools export \
            --input-path "$ART/rooted-tree.qza" \
            --output-path "$ROBJ/raspergade/final/tree/"
        info "Branch B final export: $ROBJ/raspergade/final/"
    fi

    # 27. Build phyloseq — raspergade branch
    step 27 "[B] Building phyloseq object — raspergade branch (R)"
    PHYB_R="$SCPTS/step27_phyloseq_raspergade.R"
    cat > "$PHYB_R" <<'REOF'
suppressPackageStartupMessages({
    library(phyloseq); library(biomformat); library(ape); library(tidyverse)
})

out_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "r_objects/raspergade")
meta    <- read_tsv(Sys.getenv("METADATA"), show_col_types=FALSE) |>
           column_to_rownames("sample-id")
biom_in <- read_biom(file.path(out_dir, "final/feature-table.biom"))
otu_mat <- as(biom_data(biom_in), "matrix")
tree    <- read.tree(file.path(out_dir, "final/tree/tree.nwk"))
tax_raw <- read_tsv(file.path(out_dir, "taxonomy/taxonomy.tsv"),
                    col_names=c("fid","taxon","conf"), skip=1, show_col_types=FALSE)

parse_silva <- function(x) {
    p <- str_split(x,";\\s*")[[1]]; p <- sub("^[dkpcofgs]__","",p); length(p) <- 7
    setNames(p, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
}
tax_mat <- do.call(rbind, lapply(tax_raw$taxon, parse_silva))
rownames(tax_mat) <- tax_raw$fid

shared <- intersect(colnames(otu_mat), rownames(meta))
cat("Samples matched to metadata:", length(shared), "/", ncol(otu_mat), "\n")

ps <- phyloseq(
    otu_table(otu_mat[, shared], taxa_are_rows=TRUE),
    sample_data(meta[shared, ]),
    tax_table(tax_mat[rownames(otu_mat), ]),
    phy_tree(prune_taxa(rownames(otu_mat), tree))
)
cat("Phyloseq — samples:", nsamples(ps), " ASVs:", ntaxa(ps),
    " min depth:", min(sample_sums(ps)), "\n")

saveRDS(ps, file.path(out_dir, "phyloseq_raspergade.rds"))
cat("Saved:", file.path(out_dir, "phyloseq_raspergade.rds"), "\n")
REOF
    if checkpoint 27 "$ROBJ/raspergade/phyloseq_raspergade.rds"; then
        run_r "$PHYB_R"
    fi

fi   # end Branch B

# ============================================================
# ── BRANCH C: decontam frequency (steps 28–33) ───────────────
# ============================================================

if branch_active decontam; then
    bhead "C — decontam frequency correction (steps 28–33)"

    # 28. Export Deblur table for R
    step 28 "[C] Exporting Deblur table to BIOM/TSV"
    if checkpoint 28 "$ROBJ/decontam/feature-table.biom"; then
        qrun tools export \
            --input-path "$ART/table.qza" \
            --output-path "$ROBJ/decontam/"
        biom convert \
            -i "$ROBJ/decontam/feature-table.biom" \
            -o "$ROBJ/decontam/feature-table.tsv" \
            --to-tsv
        info "Exported: $ROBJ/decontam/feature-table.biom"
    fi

    # 29. decontam frequency correction (R)
    step 29 "[C] decontam frequency method (R)"
    DECO_R="$SCPTS/step29_decontam.R"
    cat > "$DECO_R" <<'REOF'
suppressPackageStartupMessages({
    library(decontam); library(phyloseq); library(biomformat); library(tidyverse)
})

out_dir   <- file.path(Sys.getenv("OUTPUT_DIR"), "r_objects/decontam")
quant_col <- Sys.getenv("QUANT_COL")

if (nchar(quant_col) == 0)
    stop("QUANT_COL env variable is empty. Pass -q COLUMN_NAME to the workflow.")

meta    <- read_tsv(Sys.getenv("METADATA"), show_col_types=FALSE) |>
           column_to_rownames("sample-id")
biom_in <- read_biom(file.path(out_dir, "feature-table.biom"))
otu_mat <- as(biom_data(biom_in), "matrix")

if (!quant_col %in% colnames(meta))
    stop("Column '", quant_col, "' not found in metadata. Available columns:\n  ",
         paste(colnames(meta), collapse=", "))

shared <- intersect(colnames(otu_mat), rownames(meta))
otu_sub  <- otu_mat[, shared]
meta_sub <- meta[shared, ]

cat("── decontam input ─────────────────────────────────────────\n")
cat("Samples:", ncol(otu_sub), "  ASVs:", nrow(otu_sub), "\n")
cat("Concentration column:", quant_col, "\n")
cat("Samples with concentration data:", sum(!is.na(meta_sub[[quant_col]])), "\n")

ps <- phyloseq(otu_table(otu_sub, taxa_are_rows=TRUE), sample_data(meta_sub))

# Run decontam frequency method
contamdf <- isContaminant(ps, method="frequency",
                          conc=quant_col, threshold=0.1)

cat("\n── decontam results ───────────────────────────────────────\n")
cat("Total ASVs tested:        ", nrow(contamdf), "\n")
cat("Contaminants identified:  ", sum(contamdf$contaminant, na.rm=TRUE), "\n")
cat("Non-contaminants retained:", sum(!contamdf$contaminant, na.rm=TRUE), "\n")

write_tsv(rownames_to_column(contamdf, "ASV_ID") |> arrange(p),
          file.path(out_dir, "contaminant_scores.tsv"))

keep <- rownames(contamdf)[!contamdf$contaminant]
otu_clean <- otu_sub[keep, ]
write_tsv(rownames_to_column(as.data.frame(otu_clean), "#OTU ID"),
          file.path(out_dir, "decontam-table.tsv"))
saveRDS(otu_clean, file.path(out_dir, "decontam-table.rds"))
cat("Saved: decontam-table.tsv, decontam-table.rds\n")
REOF
    if checkpoint 29 "$ROBJ/decontam/decontam-table.rds"; then
        run_r "$DECO_R"
    fi

    # 30. Reimport decontam table
    step 30 "[C] Reimporting decontam-filtered table into QIIME 2"
    if checkpoint 30 "$DECO/table-decontam.qza"; then
        biom convert \
            -i "$ROBJ/decontam/decontam-table.tsv" \
            -o "$ROBJ/decontam/decontam-table.biom" \
            --table-type="OTU table" --to-hdf5
        biom summarize-table \
            -i "$ROBJ/decontam/decontam-table.biom" \
            | tee "$ROBJ/decontam/decontam-table-summary.txt"
        qrun tools import \
            --type 'FeatureTable[Frequency]' \
            --input-path "$ROBJ/decontam/decontam-table.biom" \
            --input-format BIOMV210Format \
            --output-path "$DECO/table-decontam.qza"
        qrun feature-table summarize \
            --i-table "$DECO/table-decontam.qza" \
            --m-sample-metadata-file "$METADATA" \
            --o-visualization "$DECO/table-decontam-summary.qzv"
    fi

    # 31. Core-metrics on decontam table  ← RAREFACTION
    step 31 "[C] Core-metrics-phylogenetic on decontam table (RAREFACTION at depth $DEPTH)"
    if checkpoint 31 "$DECO/diversities/weighted_unifrac_distance_matrix.qza"; then
        [[ "$FORCE" == "true" && -d "$DECO/diversities" ]] && rm -rf "$DECO/diversities"
        qrun diversity core-metrics-phylogenetic \
            --i-phylogeny "$ART/rooted-tree.qza" \
            --i-table     "$DECO/table-decontam.qza" \
            --p-sampling-depth "$DEPTH" \
            --m-metadata-file "$METADATA" \
            --output-dir "$DECO/diversities"
    fi

    # 32. Export final Branch C artifacts
    step 32 "[C] Exporting final decontam artifacts for R"
    if checkpoint 32 "$ROBJ/decontam/final/feature-table.biom"; then
        mkdir -p "$ROBJ/decontam/final"
        qrun tools export \
            --input-path "$DECO/table-decontam.qza" \
            --output-path "$ROBJ/decontam/final/"
        biom convert \
            -i "$ROBJ/decontam/final/feature-table.biom" \
            -o "$ROBJ/decontam/final/feature-table.tsv" \
            --to-tsv
        info "Branch C final export: $ROBJ/decontam/final/"
    fi

    # 33. Build phyloseq — decontam branch
    step 33 "[C] Building phyloseq object — decontam branch (R)"
    PHYC_R="$SCPTS/step33_phyloseq_decontam.R"
    cat > "$PHYC_R" <<'REOF'
suppressPackageStartupMessages({
    library(phyloseq); library(biomformat); library(ape); library(tidyverse)
})

out_dir <- file.path(Sys.getenv("OUTPUT_DIR"), "r_objects/decontam")
meta    <- read_tsv(Sys.getenv("METADATA"), show_col_types=FALSE) |>
           column_to_rownames("sample-id")
biom_in <- read_biom(file.path(out_dir, "final/feature-table.biom"))
otu_mat <- as(biom_data(biom_in), "matrix")
tree    <- read.tree(file.path(Sys.getenv("OUTPUT_DIR"),
                     "r_objects/raspergade/final/tree/tree.nwk"))
tax_raw <- read_tsv(file.path(Sys.getenv("OUTPUT_DIR"),
                    "r_objects/raspergade/taxonomy/taxonomy.tsv"),
                    col_names=c("fid","taxon","conf"), skip=1, show_col_types=FALSE)

parse_silva <- function(x) {
    p <- str_split(x,";\\s*")[[1]]; p <- sub("^[dkpcofgs]__","",p); length(p) <- 7
    setNames(p, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
}
tax_mat <- do.call(rbind, lapply(tax_raw$taxon, parse_silva))
rownames(tax_mat) <- tax_raw$fid

shared <- intersect(colnames(otu_mat), rownames(meta))
ps <- phyloseq(
    otu_table(otu_mat[, shared], taxa_are_rows=TRUE),
    sample_data(meta[shared, ]),
    tax_table(tax_mat[rownames(otu_mat), ]),
    phy_tree(prune_taxa(rownames(otu_mat), tree))
)
cat("Phyloseq — samples:", nsamples(ps), " ASVs:", ntaxa(ps),
    " min depth:", min(sample_sums(ps)), "\n")

saveRDS(ps, file.path(out_dir, "phyloseq_decontam.rds"))
cat("Saved:", file.path(out_dir, "phyloseq_decontam.rds"), "\n")
REOF
    if checkpoint 33 "$ROBJ/decontam/phyloseq_decontam.rds"; then
        run_r "$PHYC_R"
    fi

fi   # end Branch C

# ============================================================
# ── STEP 34: Cross-branch QC summary (all active branches) ───
# ============================================================

step 34 "Cross-branch QC table and phyloseq list"
QC_R="$SCPTS/step34_cross_branch_qc.R"
cat > "$QC_R" <<'REOF'
suppressPackageStartupMessages({
    library(phyloseq); library(biomformat); library(tidyverse)
})

out_dir   <- Sys.getenv("OUTPUT_DIR")
robj_dir  <- file.path(out_dir, "r_objects")
branches  <- Sys.getenv("BRANCHES")

# Discover which phyloseq RDS files exist
rds_map <- list(
    original   = file.path(robj_dir, "original/phyloseq_original.rds"),
    raspergade = file.path(robj_dir, "raspergade/phyloseq_raspergade.rds"),
    decontam   = file.path(robj_dir, "decontam/phyloseq_decontam.rds")
)

# Build original phyloseq from Branch A exports if RDS not present
if (!file.exists(rds_map$original)) {
    biom_orig <- file.path(robj_dir, "original/feature-table.biom")
    if (file.exists(biom_orig)) {
        meta  <- read_tsv(Sys.getenv("METADATA"), show_col_types=FALSE) |>
                 column_to_rownames("sample-id")
        mat   <- as(biom_data(read_biom(biom_orig)), "matrix")
        shr   <- intersect(colnames(mat), rownames(meta))
        ps    <- phyloseq(otu_table(mat[,shr], taxa_are_rows=TRUE),
                          sample_data(meta[shr,]))
        saveRDS(ps, rds_map$original)
        cat("Built original phyloseq from exported BIOM.\n")
    }
}

active <- Filter(file.exists, rds_map)
if (length(active) == 0) {
    cat("No phyloseq RDS files found yet. Run branch-specific steps first.\n")
    quit(status=0)
}

ps_list <- lapply(active, readRDS)
names(ps_list) <- names(active)

qc <- map_dfr(names(ps_list), function(b) {
    ps <- ps_list[[b]]
    tibble(branch       = b,
           n_samples    = nsamples(ps),
           n_asvs       = ntaxa(ps),
           min_depth    = min(sample_sums(ps)),
           median_depth = median(sample_sums(ps)),
           max_depth    = max(sample_sums(ps)))
})

cat("\n══ CROSS-BRANCH QC TABLE ══════════════════════════════════\n")
print(qc, n=Inf)
cat("═══════════════════════════════════════════════════════════\n\n")

saveRDS(ps_list, file.path(robj_dir, "phyloseq_all_branches.rds"))
write_tsv(qc,    file.path(robj_dir, "cross_branch_qc.tsv"))
cat("Saved:\n")
cat(" ", file.path(robj_dir, "phyloseq_all_branches.rds"), "\n")
cat(" ", file.path(robj_dir, "cross_branch_qc.tsv"), "\n")
cat("\nReady for downstream R analysis (LinDA, UniFrac, etc.).\n")
REOF
if checkpoint 34 ""; then
    BRANCHES="$BRANCHES" run_r "$QC_R"
fi

# ============================================================
# Done
# ============================================================
info "===== Revision workflow complete ====="
info "Branches run:    $BRANCHES"
info "Artifacts:       $ART"
info "Branch B output: $RASP"
info "Branch C output: $DECO"
info "R objects:       $ROBJ"
info "R scripts saved: $SCPTS"
info "Phyloseq list:   $ROBJ/phyloseq_all_branches.rds"
info "To view .qzv files: qiime tools view FILE.qzv"
exit 0
