#!/usr/bin/env bash
# ============================================================================
# build_all_bams.sh
# ----------------------------------------------------------------------------
# Build every pseudobulk BAM the TOBIAS run needs, 4 jobs abreast.
#
# THE CELL-TYPE LIST, settled 2026-08-13 from the tiering in 2026-08-11.txt
# section 3 (numbers are unique fragments in PEAKQC-passing cells):
#
#   DIFFERENTIAL, both arms -> one BINDetect per cell type
#     Microglia         255.7M / 28.1M
#     Oligodendrocytes  159.5M / 59.6M
#     Neurons_V          84.3M / 33.9M
#     Astrocytes         53.2M / 29.5M
#     Neurons_D          49.5M / 18.7M
#     OPCs               15.9M /  2.6M   asymmetric but workable
#
#   INJURED-ONLY AGGREGATE -- no contrast possible
#     Ependymal_I        75.2M   U is 3.4M, and 714 of its 755 cells come from
#                                FT_U_4, the 2nd-worst library in the dataset.
#                                That arm is one bad sample, not a shallow arm.
#     Macrophages_I      51.8M   U is 22 cells.
#
#   POOLED AGGREGATE ONLY -- too shallow to split
#     Perivascular        3.1M
#     Endothelial         1.8M   ALREADY BUILT, not in this list
#
# FT_7dpi_4 and FT_U_4 are handled by per-cell PEAKQC filtering (fld_score
# >= 100), NOT dropped wholesale -- so the fragment counts above stand as
# computed in qc/injury_viability_unique.csv, with no recount needed.
#
# SCHEDULING. Jobs run largest-first across 4 slots (longest-processing-time,
# which minimises makespan). Measured on this machine 2026-08-13:
#   scan  1.95M fragment lines/s -> ~15 min for one pass over all 21 samples
#   write 92,752 fragments/s     -> the real bottleneck, ~2.8h of work total
# Each job re-scans all 21 GB; that is redundant but it costs cores, not wall
# time, and it keeps the empirically validated script byte-for-byte unchanged.
# Expect ~2h wall clock and ~35 GB of final BAMs (604 GB free at launch).
#
# --threads 3 per job so 4 concurrent samtools sorts stay under 14 cores.
#
# Resumable: a job whose .bam already exists is skipped, so a crash or a Ctrl-C
# costs one BAM rather than the run.
#
# Usage (from the tobias repo root):
#   bash pipeline/build_all_bams.sh
#   PARALLEL=6 bash pipeline/build_all_bams.sh     # more slots
# Logs: logs/bam_<tag>.log, one per job.
# ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."

PY="${PY:-$HOME/miniconda3/envs/tobias/bin/python}"
PARALLEL="${PARALLEL:-4}"
THREADS="${THREADS:-3}"

mkdir -p bam logs

# celltype:injury  ("-" means pool both arms). Largest first.
JOBS=$(cat <<'EOF'
Microglia:I
Oligodendrocytes:I
Neurons_V:I
Ependymal:I
Oligodendrocytes:U
Astrocytes:I
Macrophages:I
Neurons_D:I
Neurons_V:U
Astrocytes:U
Microglia:U
Neurons_D:U
OPCs:I
Perivascular:-
OPCs:U
EOF
)

run_one() {
    local ct="$1" inj="$2"
    local tag="$ct" args=(--celltype "$ct")
    if [ "$inj" != "-" ]; then
        tag="${ct}_${inj}"
        args+=(--injury "$inj")
    fi
    if [ -s "bam/${tag}.bam" ]; then
        echo "$(date +%H:%M:%S) | SKIP ${tag} (bam/${tag}.bam exists)"
        return 0
    fi
    echo "$(date +%H:%M:%S) | START ${tag}"
    if "$PY" pipeline/fragments_to_bam.py "${args[@]}" --threads "$THREADS" \
            > "logs/bam_${tag}.log" 2>&1; then
        echo "$(date +%H:%M:%S) | DONE  ${tag}  ($(du -h "bam/${tag}.bam" | cut -f1))"
    else
        echo "$(date +%H:%M:%S) | FAIL  ${tag}  -- see logs/bam_${tag}.log"
        return 1
    fi
}
export -f run_one
export PY THREADS

echo "$(date +%H:%M:%S) | launching $(echo "$JOBS" | wc -l | tr -d ' ') jobs, ${PARALLEL} abreast"
# A failing job must not abort the summary -- xargs returns 123 if any child
# failed, and the surviving BAMs still need to be listed.
status=0
printf '%s\n' "$JOBS" \
    | xargs -P "$PARALLEL" -I{} bash -c 'run_one "${1%%:*}" "${1##*:}"' _ {} || status=$?

echo "$(date +%H:%M:%S) | ALL JOBS FINISHED (xargs status ${status})"
ls -lh bam/*.bam
exit "$status"
