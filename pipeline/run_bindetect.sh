#!/usr/bin/env bash
# ============================================================================
# run_bindetect.sh
# ----------------------------------------------------------------------------
# BINDetect for the SIX cell types that have both an injured and an uninjured
# arm: injured vs uninjured, one invocation per cell type.
#
# *** WHY ONE INVOCATION PER CELL TYPE ***
# BINDetect compares ALL-AGAINST-ALL by default. Handing it all twelve bigwigs
# would produce 66 pairwise comparisons, almost all of them meaningless
# (Microglia_I vs Astrocytes_U), and would force a single shared peak set on
# cell types whose peaks were called separately. One run per cell type keeps
# each comparison within its own peak set, which is the only set where both
# arms are jointly valid.
#
# *** THE FOUR SINGLE-ARM GROUPS -- SETTLED 2026-08-20, SPLIT 2/2 ***
# All four lack a second arm, but for different reasons and with different
# consequences. The decision was to RUN two and DROP two.
#
# RUN, occupancy only (no log2fc, no differential claim):
#   Ependymal_I    75.2M fragments -- the 3rd-deepest injured arm in the whole
#                  dataset, deeper than Astrocytes_I, Neurons_D_I and OPCs_I.
#   Macrophages_I  51.8M fragments -- comparable to Astrocytes_I.
# Both name a real CELL STATE ("reactive ependymal cells"), which is what makes
# an occupancy statement interpretable and comparable to the paper. Zamboni
# Fig. 4e reports AP-1 partnering with NFI and RFX in ependymal cells, from
# models TRAINED ON INJURED SAMPLES -- i.e. exactly what Ependymal_I is. So the
# TF-COMB co-occurrence run on it tests a published claim directly, and it does
# NOT inherit the bad-library problem, which lives in the UNINJURED arm.
#
# DROPPED from the biological analysis:
#   Perivascular   3.1M   } injured+uninjured POOLED into one BAM and
#   Endothelial    1.8M   } unsplittable. Their signal is AVERAGED ACROSS THE
# VERY INJURY AXIS the project is about, and no re-run recovers it -- the
# information was destroyed when the BAM was built. "Endothelial cells, injured
# and uninjured blended" is not a cell state, so it is comparable to nothing in
# the paper, whose per-cell-type motif analyses all used injured-trained models.
# NOTE the disqualifier is the POOLING, not depth: Perivascular (3.1M) is
# actually deeper than OPCs_U (2.6M), which we treat as first-class.
# Endothelial is retained as the pipeline's technical canary -- it was the
# cell type validated end-to-end on 2026-08-11 -- but not as a result.
#
# *** WHAT A SINGLE-CONDITION RUN ACTUALLY MEANS ***
# Not "self-referential", which was the initial worry. bindetect.py:462-557
# derives the bound threshold from BACKGROUND regions, not from the motif
# sites: it fits a lognormal to background footprint scores, mirrors the left
# side to build a null free of the bound tail, fits a normal to that, and takes
# threshold = norm.ppf(1 - bound_pvalue). So "bound" means "exceeds what
# unbound background chromatin in this same sample produces at p < 0.001" --
# foreground vs background, a real baseline. bindetect.py:426 shows the single-
# condition path is explicitly anticipated: len(cond_names) == 1 selects
# constant normalization instead of quantile.
# What genuinely does not exist with one arm is the CONTRAST. Also, sensitivity
# is depth-dependent (shallower -> wider background -> higher threshold -> fewer
# bound sites), so bound-site COUNTS are not comparable across samples.
#
# *** SIGNAL ORDER SETS THE SIGN OF log2fc ***
# bindetect_functions.py:429 computes
#     log2fc = log2((cond1_score + pseudo) / (cond2_score + pseudo))
# so the FIRST --signals entry is the numerator. Injured is passed first, so
#     positive log2fc = more bound AFTER injury.
# Column names in the output follow <cond1>_<cond2>_log2fc, i.e.
# `<ct>_I_<ct>_U_log2fc`.
#
# PEAKS. Both arms are scored over the cell type's own peak set,
# <ct>_peaks_no_blacklist_chr.bed -- peaks were called per cell type, not per
# condition, so this is the set on which both arms are jointly defined.
#
# *** MOTIF SET: JASPAR2020 CORE vertebrates, non-redundant, 746 motifs ***
# meta/motifs/JASPAR2020_CORE_vertebrates_non-redundant.meme
#
# NOT the sibling repo's meta/motifs.meme.txt (2,193 motifs). That file was
# assembled for the ChromBPNet track's TF-MoDISco MATCHING and is the wrong
# input for BINDetect SCANNING, for two reasons:
#   1. It contains 8 TN5_* and 6 DNASE_* sequencing-artifact motifs. They exist
#      so MoDISco can recognise Tn5 insertion bias rather than a real factor.
#      BINDetect has no such concept and would report them as transcription
#      factors, with differential scores, in <ct>_results.txt. A correctness
#      problem, not noise. JASPAR2020 CORE contains none of them (verified).
#   2. At 2,193 motifs the per-TF phase projected to ~12h per cell type
#      (~36h for six). 746 motifs is ~3x less work.
#
# WHY JASPAR2020 SPECIFICALLY: it is what Zamboni et al. used for their
# conventional motif enrichment (Methods, "the 'CORE' vertebrate collection in
# JASPAR2020 database", via Signac AddMotifs()), so our differential binding
# results are directly comparable to the paper's supervised analysis. The
# non-redundant collection matches Signac's default, which keeps only the
# latest version of each matrix.
#
# PROVENANCE: JASPAR no longer serves the 2020 release (jaspar.elixir.no hosts
# 2022 onward). Obtained from the MEME Suite motif database archive
# motif_databases.12.23.tgz -> motif_databases/JASPAR/, which ships it verbatim.
#
# --skip-excel is on: 746 motifs would otherwise mean 746 .xlsx files per cell
# type. The .txt equivalents are written regardless.
#
# Completion is tested by TOBIAS's own final log line, "Finished BINDetect
# run", never by whether an output file exists -- see run_tobias_chain.sh for
# what that mistake cost us on 2026-08-13.
#
# Always invoked through pipeline/tobias_fork.py -- TOBIAS cannot run directly
# on macOS (its multiprocessing assumes fork; spawn dies on an unpicklable
# logger).
#
# *** THE DAR SECONDARY ANALYSIS (TASKSETS=dars) ***
# A paper-matched comparison. Zamboni et al. measured motif enrichment WITHIN
# differentially accessible regions (IRENs); our primary analysis measures
# differential binding over the cell type's FULL peak set. Different
# denominators, so only rank agreement is meaningful between them.
#
# THE FULL PEAK SET IS STILL PASSED TO --peaks. This is not optional and it is
# not a style choice -- BINDetect's own help says so:
#   "--output-peaks ... This will limit all analysis to the regions in
#    --output-peaks. NOTE: --peaks must still be set to the full peak set!"
# The reason is visible in the source: background is sampled from within the
# --peaks regions (bindetect_functions.py:266+ scan_and_score draws
# reglen/200 random positions per region), and that background is what
# calibrates (a) the bound/unbound threshold, (b) the cross-condition quantile
# normalization and (c) the differential log2fc null. Passing only DARs would
# calibrate all three on regions pre-selected for differing between conditions
# -- in particular it would quantile-normalize injured against uninjured using
# only regions chosen because they differ, suppressing the signal. The subset
# is applied later, in process_tfbs (bindetect_functions.py:381-385), i.e.
# AFTER those estimates are made. So: statistics on the full set, reporting on
# the subset. Exactly what we want.
#
# INTERPRETING DAR RESULTS -- READ THIS BEFORE QUOTING THEM.
# prepare_ml_dars.R selects DARs with `avg_log2FC > 0.5`, i.e. ONE-SIDED,
# regions OPENING after injury. Footprint scores scale with accessibility, so
# within a set selected for gaining accessibility, "binding went up" is partly
# built in. These results are therefore valid as a RANKING among opening
# regions -- which TFs dominate the injury-opened enhancers -- and are NOT
# evidence of direction. That is the same claim shape the paper makes.
# Two further caveats, traced 2026-08-20: the DARs were called on the AUTHORS'
# unfiltered cell set (config ml.merged_obj, pre-PEAKQC), and the threshold is
# on the UNADJUSTED p-value (dars_max_p applied as `p_val < 0.05`), so the set
# is generous rather than stringent.
#
# Usage (from the tobias repo root):
#   bash pipeline/run_bindetect.sh                  # primary: pair + solo
#   TASKSETS=dars bash pipeline/run_bindetect.sh    # DAR secondary analysis
#   PARALLEL=1 CORES=12 bash pipeline/run_bindetect.sh
# Logs: logs/bindetect_<ct>.log
# Output: tobias/bindetect/<ct>/
# ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."

PY="${PY:-$HOME/miniconda3/envs/tobias/bin/python}"

# *** bedtools MUST BE ON PATH FOR --output-peaks ***
# BINDetect subsets to --output-peaks via pybedtools BedTool.intersect
# (bindetect_functions.py:383-385), which shells out to the bedtools binary.
# The tobias env does not ship it, so the DAR runs die mid-flight with:
#   NotImplementedError: "intersectBed" does not appear to be installed or on
#   the path, so this method is disabled.
# The failure surfaces from a child process ~2.5 min in, AFTER motif scanning,
# so it is not caught by any up-front check. Hit 2026-08-20.
BEDTOOLS_BIN="${BEDTOOLS_BIN:-$HOME/miniconda3/envs/bedtools/bin}"
if [ -x "$BEDTOOLS_BIN/bedtools" ]; then
    export PATH="$BEDTOOLS_BIN:$PATH"
fi
CB="${CB:-../TF_cobinding_multiome_chrombpnet}"
PARALLEL="${PARALLEL:-2}"
CORES="${CORES:-6}"

GENOME="$CB/meta/mm10.fa"
MOTIFS="${MOTIFS:-meta/motifs/JASPAR2020_CORE_vertebrates_non-redundant.meme}"
SRC_PEAKDIR="$CB/ML/filtered_peaks"
PEAKDIR="meta/peaks_tobias"

mkdir -p tobias/bindetect logs "$PEAKDIR"

# *** PEAKS MUST BE FILTERED TO THE CONTIGS OUR BIGWIGS ACTUALLY CONTAIN ***
# BINDetect hard-errors if a peak sits on a chromosome absent from the signal:
#   ERROR Chromosome for region "chrM 139 1579 ..." is not found in list of
#         available chromosomes ([chr1 ... chr9, chrX, chrY])
# Hit on 2026-08-20: Astrocytes carries 6 chrM peaks and Ependymal 3. Our
# pseudobulk BAMs (and hence the bigwigs) cover chr1-19, chrX, chrY only.
# chrM ATAC peaks are artifacts and are conventionally excluded anyway, and 9
# peaks out of ~264,000 across the two files is negligible.
#
# Filtering is applied UNIFORMLY to every cell type, not just the two that
# fail, so all cell types are treated identically. For the six files with no
# odd contigs the output is byte-identical to the source, so Microglia and
# Oligodendrocytes -- which completed before this step existed -- are
# unaffected and do not need re-running.
prepare_peaks() {
    local src filtered base n_src n_out
    for src in "$SRC_PEAKDIR"/*_peaks_no_blacklist_chr.bed; do
        base=$(basename "$src")
        filtered="${PEAKDIR}/${base}"
        if [ -f "$filtered" ] && [ "$filtered" -nt "$src" ]; then continue; fi
        awk 'BEGIN{FS=OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/' "$src" > "$filtered"
        n_src=$(wc -l < "$src" | tr -d ' ')
        n_out=$(wc -l < "$filtered" | tr -d ' ')
        if [ "$n_src" != "$n_out" ]; then
            echo "$(date +%H:%M:%S) | PEAKS ${base}: dropped $((n_src - n_out)) off-contig peaks ($n_src -> $n_out)"
        fi
    done
}
prepare_peaks

# Cell types with BOTH arms -> differential run. See header.
PAIRED="Microglia Oligodendrocytes Neurons_V Astrocytes Neurons_D OPCs"

# Injured-only groups -> occupancy-only run, ONE signal. See header.
# Perivascular and Endothelial are deliberately NOT here.
SOLO="Ependymal_I Macrophages_I"

finished() {   # finished <logfile>
    [ -f "$1" ] && grep -q "Finished BINDetect run" "$1"
}

run_one() {
    local ct="$1"
    local outdir="tobias/bindetect/${ct}"
    local peaks="${PEAKDIR}/${ct}_peaks_no_blacklist_chr.bed"
    local inj="tobias/footprints/${ct}_I_footprints.bw"
    local uninj="tobias/footprints/${ct}_U_footprints.bw"
    local log="logs/bindetect_${ct}.log"

    for f in "$peaks" "$inj" "$uninj"; do
        if [ ! -f "$f" ]; then
            echo "$(date +%H:%M:%S) | FAIL  ${ct} -- missing ${f}"; return 1
        fi
    done

    if finished "$log"; then
        echo "$(date +%H:%M:%S) | SKIP  ${ct} (completed earlier)"; return 0
    fi

    echo "$(date +%H:%M:%S) | BINDETECT ${ct}"
    mkdir -p "$outdir"
    # Injured first -- sets positive log2fc = up after injury. See header.
    if ! "$PY" pipeline/tobias_fork.py BINDetect \
            --signals "$inj" "$uninj" \
            --cond-names "${ct}_I" "${ct}_U" \
            --motifs "$MOTIFS" --genome "$GENOME" --peaks "$peaks" \
            --outdir "$outdir" --prefix "$ct" \
            --skip-excel --cores "$CORES" \
            > "$log" 2>&1; then
        echo "$(date +%H:%M:%S) | FAIL  ${ct} -- ${log}"
        return 1
    fi
    echo "$(date +%H:%M:%S) | DONE  ${ct}"
}
# Occupancy-only: ONE signal, ONE --cond-names. No log2fc columns are produced;
# the output carries <ct>_threshold and <ct>_bound instead. The per-TF
# beds/<TF>_<ct>_bound.bed files are still written, which is what TF-COMB's
# TFBS_from_TOBIAS(bindetect_path, condition) consumes.
run_one_solo() {
    local ct="$1"                    # e.g. Ependymal_I
    local base="${ct%_I}"            # peaks are per cell type, not per condition
    local outdir="tobias/bindetect/${ct}"
    local peaks="${PEAKDIR}/${base}_peaks_no_blacklist_chr.bed"
    local sig="tobias/footprints/${ct}_footprints.bw"
    local log="logs/bindetect_${ct}.log"

    for f in "$peaks" "$sig"; do
        if [ ! -f "$f" ]; then
            echo "$(date +%H:%M:%S) | FAIL  ${ct} -- missing ${f}"; return 1
        fi
    done

    if finished "$log"; then
        echo "$(date +%H:%M:%S) | SKIP  ${ct} (completed earlier)"; return 0
    fi

    echo "$(date +%H:%M:%S) | BINDETECT ${ct} (occupancy only, single signal)"
    mkdir -p "$outdir"
    if ! "$PY" pipeline/tobias_fork.py BINDetect \
            --signals "$sig" \
            --cond-names "${ct}" \
            --motifs "$MOTIFS" --genome "$GENOME" --peaks "$peaks" \
            --outdir "$outdir" --prefix "$ct" \
            --skip-excel --cores "$CORES" \
            > "$log" 2>&1; then
        echo "$(date +%H:%M:%S) | FAIL  ${ct} -- ${log}"
        return 1
    fi
    echo "$(date +%H:%M:%S) | DONE  ${ct}"
}

# Paper-matched secondary analysis: identical to run_one EXCEPT that results
# are limited to the cell type's injury-opened DARs via --output-peaks, and
# output goes to a SEPARATE directory so the primary results are untouched.
run_one_dars() {
    local ct="$1"
    if ! command -v bedtools >/dev/null 2>&1; then
        echo "$(date +%H:%M:%S) | FAIL  ${ct} (dars) -- bedtools not on PATH; --output-peaks needs it"
        return 1
    fi
    local outdir="tobias/bindetect_dars/${ct}"
    local peaks="${PEAKDIR}/${ct}_peaks_no_blacklist_chr.bed"   # FULL set -- see header
    local dars="${CB}/ML/celltype_dars/${ct}_I_dars.bed"
    local inj="tobias/footprints/${ct}_I_footprints.bw"
    local uninj="tobias/footprints/${ct}_U_footprints.bw"
    local log="logs/bindetect_dars_${ct}.log"

    for f in "$peaks" "$dars" "$inj" "$uninj"; do
        if [ ! -f "$f" ]; then
            echo "$(date +%H:%M:%S) | FAIL  ${ct} (dars) -- missing ${f}"; return 1
        fi
    done

    if finished "$log"; then
        echo "$(date +%H:%M:%S) | SKIP  ${ct} (dars, completed earlier)"; return 0
    fi

    echo "$(date +%H:%M:%S) | BINDETECT ${ct} (DAR subset: $(wc -l < "$dars" | tr -d ' ') regions)"
    mkdir -p "$outdir"
    if ! "$PY" pipeline/tobias_fork.py BINDetect \
            --signals "$inj" "$uninj" \
            --cond-names "${ct}_I" "${ct}_U" \
            --motifs "$MOTIFS" --genome "$GENOME" \
            --peaks "$peaks" --output-peaks "$dars" \
            --outdir "$outdir" --prefix "$ct" \
            --skip-excel --cores "$CORES" \
            > "$log" 2>&1; then
        echo "$(date +%H:%M:%S) | FAIL  ${ct} (dars) -- ${log}"
        return 1
    fi
    echo "$(date +%H:%M:%S) | DONE  ${ct} (dars)"
}

# Dispatch on a mode prefix so all kinds share one scheduler.
run_task() {
    case "$1" in
        pair:*) run_one "${1#pair:}" ;;
        solo:*) run_one_solo "${1#solo:}" ;;
        dars:*) run_one_dars "${1#dars:}" ;;
        *) echo "unknown task: $1"; return 1 ;;
    esac
}

export -f run_one run_one_solo run_one_dars run_task finished
export PY CB CORES GENOME MOTIFS PEAKDIR

TASKSETS="${TASKSETS:-pair solo}"
TASKS=""
for set in $TASKSETS; do
    case "$set" in
        pair) for ct in $PAIRED; do TASKS="$TASKS pair:$ct"; done ;;
        solo) for ct in $SOLO;   do TASKS="$TASKS solo:$ct"; done ;;
        dars) for ct in $PAIRED; do TASKS="$TASKS dars:$ct"; done ;;
        *) echo "unknown taskset: $set"; exit 1 ;;
    esac
done

echo "$(date +%H:%M:%S) | tasksets [${TASKSETS}]: $(echo $TASKS | wc -w | tr -d ' ') runs, ${PARALLEL} abreast at ${CORES} cores"
echo "$(date +%H:%M:%S) | motifs: ${MOTIFS} ($(grep -c '^MOTIF' "$MOTIFS") entries)"
status=0
printf '%s\n' $TASKS | xargs -P "$PARALLEL" -I{} bash -c 'run_task "$1"' _ {} || status=$?

echo "$(date +%H:%M:%S) | BINDETECT FINISHED (xargs status ${status})"
exit "$status"
