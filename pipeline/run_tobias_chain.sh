#!/usr/bin/env bash
# ============================================================================
# run_tobias_chain.sh
# ----------------------------------------------------------------------------
# ATACorrect -> validate -> ScoreBigwig, for every BAM in bam/.
#
# BINDetect is NOT run here. It is all-against-all by default, so it needs ONE
# invocation per cell type pairing the two arms, and the injured-only and
# aggregate-only groups have no pairing at all. Run it by hand once the
# footprint bigwigs are in place.
#
# *** THE VALIDATION GATE ***
# ScoreBigwig only runs if pipeline/check_atacbias.py passes on that run's
# *_AtacBias.pickle. Tn5's bias motif is palindromic about the insertion point,
# so the learned motif says where ATACorrect actually thought the cut site was.
# Offset 0 = correct. Offset +-9 = the 9 bp double shift, which produces
# plausible-looking output with every footprint smeared and NO error message.
# This is the failure this whole track was most at risk of; the gate costs a
# second per run.
#
# *** WHY COMPLETION IS TESTED VIA THE LOG, NOT THE OUTPUT FILE ***
# ATACorrect creates <tag>_corrected.bw EARLY and grows it throughout the run --
# measured 2026-08-13: Oligodendrocytes_I_corrected.bw was already 37 MB while
# the run was only 46% done. So `[ -s "$corrected" ]` is NOT a completion test:
# an interrupted run leaves a non-empty but TRUNCATED bigwig, and a re-run would
# skip ATACorrect and hand truncated data to the gate and to ScoreBigwig.
# TOBIAS writes "Finished ATACorrect run" / "Finished ScoreBigwig run" as its
# last log line only on success, so that is what we test. It is also backward
# compatible with runs completed before this fix.
#
# PARAMETERS. Everything is at TOBIAS defaults except --read_shift 0 0, because
# CellRanger has already applied the +4/-5 Tn5 shift and ATACorrect would apply
# it a second time. TOBIAS_snakemake's own example config runs ATACorrect,
# ScoreBigwig and BINDetect with empty parameter strings, so defaults are the
# authors' recommendation. See 2026-08-11.txt section 4.
# --score footprint is passed explicitly: it IS the default, but this is the
# parameter the whole analysis turns on, so it should not be implicit.
#
# SPLIT. TOBIAS chunks regions for multiprocessing by REGION COUNT, not by reads
# or bp (utils/regions.py, RegionList.chunks). Non-peak regions -- where bias is
# estimated -- run from tiny inter-peak gaps to multi-Mb gene deserts, so with
# the default --split 100 one chunk can carry a wildly disproportionate read
# load and set the finish time on its own. Raising SPLIT gives finer chunks and
# better load balancing. It changes only how work is divided, never the result.
# Left at the TOBIAS default of 100 unless overridden, since the benefit here is
# inferred from reading chunks() and has not been measured.
#
# Always invoked through pipeline/tobias_fork.py -- TOBIAS cannot run directly
# on macOS (its multiprocessing assumes fork; spawn dies on an unpicklable
# logger).
#
# SCHEDULING. 2 jobs abreast at --cores 6 on 14 cores. Ordered largest-first,
# so Microglia_I (255M fragments) occupies one slot while the other churns
# through the rest -- which also means validation results from several cell
# types arrive early rather than after everything finishes.
#
# Usage (from the tobias repo root):
#   bash pipeline/run_tobias_chain.sh
#   PARALLEL=3 CORES=4 SPLIT=1000 bash pipeline/run_tobias_chain.sh
# Logs: logs/atacorrect_<tag>.log, logs/scorebigwig_<tag>.log
# ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."

PY="${PY:-$HOME/miniconda3/envs/tobias/bin/python}"
CB="${CB:-../TF_cobinding_multiome_chrombpnet}"
PARALLEL="${PARALLEL:-2}"
CORES="${CORES:-6}"
SPLIT="${SPLIT:-100}"

GENOME="$CB/meta/mm10.fa"
BLACKLIST="$CB/meta/mm10-blacklist.v2.bed"
PEAKDIR="$CB/ML/filtered_peaks"

mkdir -p tobias/atacorrect tobias/footprints logs

# Largest first (fragment counts from bam/<tag>.stats.txt).
JOBS="Microglia_I Oligodendrocytes_I Neurons_V_I Ependymal_I Oligodendrocytes_U
Astrocytes_I Macrophages_I Neurons_D_I Neurons_V_U Astrocytes_U Microglia_U
Neurons_D_U OPCs_I Perivascular OPCs_U Endothelial"

# TOBIAS writes this as its final line only on a successful run.
finished() {   # finished <logfile> <tool>
    [ -f "$1" ] && grep -q "Finished $2 run" "$1"
}

run_one() {
    local tag="$1"
    local ct="${tag%_I}"; ct="${ct%_U}"          # peaks are per cell type
    local outdir="tobias/atacorrect/${tag}"
    local peaks="${PEAKDIR}/${ct}_peaks_no_blacklist_chr.bed"
    local corrected="${outdir}/${tag}_corrected.bw"
    local fp="tobias/footprints/${tag}_footprints.bw"
    local aclog="logs/atacorrect_${tag}.log"
    local sblog="logs/scorebigwig_${tag}.log"

    if [ ! -f "$peaks" ]; then
        echo "$(date +%H:%M:%S) | FAIL  ${tag} -- no peaks at ${peaks}"; return 1
    fi

    if finished "$aclog" ATACorrect && [ -s "$corrected" ]; then
        echo "$(date +%H:%M:%S) | SKIP  ${tag} atacorrect (completed earlier)"
    else
        echo "$(date +%H:%M:%S) | ATACORRECT ${tag}"
        mkdir -p "$outdir"
        if ! "$PY" pipeline/tobias_fork.py ATACorrect \
                --bam "bam/${tag}.bam" --genome "$GENOME" --peaks "$peaks" \
                --blacklist "$BLACKLIST" --read_shift 0 0 \
                --outdir "$outdir" --prefix "$tag" \
                --cores "$CORES" --split "$SPLIT" \
                > "$aclog" 2>&1; then
            echo "$(date +%H:%M:%S) | FAIL  ${tag} atacorrect -- ${aclog}"
            return 1
        fi
    fi

    # ---- the gate -------------------------------------------------------
    if ! "$PY" pipeline/check_atacbias.py "${outdir}/${tag}_AtacBias.pickle"; then
        echo "$(date +%H:%M:%S) | GATE FAILED ${tag} -- bias motif is not centred on"
        echo "             our cut sites. NOT running ScoreBigwig. Investigate before"
        echo "             trusting anything downstream."
        return 1
    fi

    if finished "$sblog" ScoreBigwig && [ -s "$fp" ]; then
        echo "$(date +%H:%M:%S) | SKIP  ${tag} scorebigwig (completed earlier)"
        return 0
    fi
    echo "$(date +%H:%M:%S) | SCOREBIGWIG ${tag}"
    if ! "$PY" pipeline/tobias_fork.py ScoreBigwig \
            --signal "$corrected" --regions "$peaks" --output "$fp" \
            --score footprint --cores "$CORES" --split "$SPLIT" \
            > "$sblog" 2>&1; then
        echo "$(date +%H:%M:%S) | FAIL  ${tag} scorebigwig -- ${sblog}"
        return 1
    fi
    echo "$(date +%H:%M:%S) | DONE  ${tag}"
}
export -f run_one finished
export PY CB CORES SPLIT GENOME BLACKLIST PEAKDIR

echo "$(date +%H:%M:%S) | $(echo $JOBS | wc -w | tr -d ' ') cell types, ${PARALLEL} abreast at ${CORES} cores, split ${SPLIT}"
status=0
printf '%s\n' $JOBS | xargs -P "$PARALLEL" -I{} bash -c 'run_one "$1"' _ {} || status=$?

echo "$(date +%H:%M:%S) | CHAIN FINISHED (xargs status ${status})"
echo
echo "=== validation summary ==="
"$PY" pipeline/check_atacbias.py --quiet tobias/atacorrect/*/*_AtacBias.pickle || true
exit "$status"
