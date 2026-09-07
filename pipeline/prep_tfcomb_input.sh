#!/usr/bin/env bash
# ============================================================================
# prep_tfcomb_input.sh -- concatenate BINDetect bound sites into TF-COMB input
# ----------------------------------------------------------------------------
# TF-COMB's market-basket model wants ONE bed of all TFBS, with the TF name in
# column 4; BINDetect writes one bed per motif per arm. This flattens the
# 746-motif tree for a given <celltype>_<arm> into a single sorted 6-column bed.
#
# *** ALL 746 MOTIFS, ALWAYS ***
# Never prepare a partial motif set. TF-COMB ranks each pair against the
# background of all other pairs, so dropping motifs does not merely remove rows
# -- it changes the scores of the rows that remain, and therefore the ranks that
# Tier 3 is scored on.
#
# *** BOUND SITES ONLY ***
# The `_unbound.bed` and `_all.bed` companions are deliberately not used: the
# co-occurrence question is about sites the footprinting called BOUND. (`_all.bed`
# is also the file PlotAggregate reads and the one absent from the file server.)
#
# Columns are cut to BED6 (chrom/start/end/name/score/strand) -- BINDetect's
# 17-column output carries the peak record it overlapped, which TF-COMB ignores
# and which would triple the file size.
#
# Usage:  bash pipeline/prep_tfcomb_input.sh <celltype> <arm> [outdir]
#   e.g.  bash pipeline/prep_tfcomb_input.sh Astrocytes I
#         bash pipeline/prep_tfcomb_input.sh Ependymal_I I     # single-arm group
# ============================================================================
set -euo pipefail

CT="${1:?usage: prep_tfcomb_input.sh <celltype> <arm> [outdir]}"
ARM="${2:?usage: prep_tfcomb_input.sh <celltype> <arm> [outdir]}"
OUT="${3:-tobias/tfcomb/input}"
BD="${BINDETECT:-tobias/bindetect}"

# Ependymal_I / Macrophages_I already carry the arm in the group name; the bed
# files are <motif>_<group>_bound.bed, so the tag is the group name itself.
case "$CT" in
    *_I|*_U) TAG="$CT" ;;
    *)       TAG="${CT}_${ARM}" ;;
esac

mkdir -p "$OUT"
DEST="$OUT/${TAG}.bed"

# macOS ships bash 3.2, which has no `mapfile`; a file list is used instead
# (and is what xargs wants anyway).
LIST="$(mktemp -t tfcomb_beds)"
trap 'rm -f "$LIST"' EXIT
find "$BD/$CT" -name "*_${TAG}_bound.bed" | sort > "$LIST"
N=$(wc -l < "$LIST" | tr -d ' ')
if [ "$N" -ne 746 ]; then
    echo "ERROR: found $N bound beds for $TAG, expected 746. Refusing to build a"
    echo "       partial motif set -- it would change every rank Tier 3 scores."
    exit 1
fi

echo "==> $TAG: flattening $N motif beds"
xargs cat < "$LIST" \
  | cut -f1-6 \
  | sort -k1,1 -k2,2n \
  > "$DEST"

SITES=$(wc -l < "$DEST")
TFS=$(cut -f4 "$DEST" | sort -u | wc -l)
echo "    $DEST  sites=$SITES  unique TFs=$TFS"
[ "$TFS" -eq 746 ] || { echo "ERROR: $TFS distinct TF names in output, expected 746"; exit 1; }
