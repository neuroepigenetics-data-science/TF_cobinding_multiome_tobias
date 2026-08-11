#!/usr/bin/env bash
# ============================================================================
# count_unique_fragments.sh
# ----------------------------------------------------------------------------
# Per-barcode count of UNIQUE FRAGMENTS (= lines) in each sample's
# atac_fragments.tsv.gz.
#
# WHY THIS EXISTS. PEAKQC's `n_fragments` is NOT a fragment count -- it sums
# column 5 of the fragments file, which is the number of read pairs supporting
# each fragment (i.e. it includes PCR duplicates). Measured on Endothelial:
#     BED lines (unique fragments)   1,808,474   <- what a synthesised BAM holds
#     sum of column 5               5,188,752
#     PEAKQC n_fragments            5,197,437   <- matches column 5, not lines
# The synthesised BAM gets ONE READ PAIR PER LINE (expanding by the count column
# would re-introduce the PCR duplicates CellRanger already removed), so the line
# count is the depth TOBIAS will actually see. Using PEAKQC's number overstates
# every cell type by ~2.7x and would corrupt the viability tiering.
#
# Per-cell counts are needed, not a global scale factor: duplicate rate varies
# by cell and by sample, and the tight calls (Ependymal_U) are exactly where a
# global correction would be least trustworthy.
#
# Checkpointed per sample -- re-run to resume. Comment lines ('#') are skipped;
# column 4 is the barcode.
#
# Usage (from the tobias repo root):
#   bash pipeline/count_unique_fragments.sh [n_parallel]
# Output:
#   qc/fragcounts/<sample>.tsv    barcode <TAB> n_unique_fragments
# ============================================================================
set -u -o pipefail

CHROMBPNET="../TF_cobinding_multiome_chrombpnet"
OUTDIR="qc/fragcounts"
PAR="${1:-4}"

mkdir -p "$OUTDIR"

count_one() {
    local frag="$1"
    local sample
    sample=$(basename "$(dirname "$(dirname "$frag")")")
    local out="$OUTDIR/$sample.tsv"
    local tmp="$out.partial"

    if [[ -s "$out" ]]; then
        echo "$(date +%H:%M:%S) | $sample: checkpoint exists, skipping"
        return 0
    fi

    echo "$(date +%H:%M:%S) | $sample: counting"
    # Write to .partial and rename only on success, so an interrupted run
    # never leaves a truncated checkpoint that the skip-test would trust.
    if bgzip -c -d -@ 2 "$frag" \
        | awk -F'\t' '$1 !~ /^#/ {c[$4]++} END {for (b in c) printf "%s\t%d\n", b, c[b]}' \
        > "$tmp"; then
        mv "$tmp" "$out"
        echo "$(date +%H:%M:%S) | $sample: done, $(wc -l < "$out" | tr -d ' ') barcodes"
    else
        rm -f "$tmp"
        echo "$(date +%H:%M:%S) | $sample: FAILED"
        return 1
    fi
}
export -f count_one
export CHROMBPNET OUTDIR

# Smallest first, so a failure surfaces early (same convention as run_peakqc.py).
ls -S -r "$CHROMBPNET"/data/*/outs/atac_fragments.tsv.gz \
    | xargs -P "$PAR" -I{} bash -c 'count_one "$@"' _ {}

echo "$(date +%H:%M:%S) | ALL DONE: $(ls "$OUTDIR"/*.tsv 2>/dev/null | wc -l | tr -d ' ') / 21 samples"
