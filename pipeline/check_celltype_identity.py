#!/usr/bin/env python
# ============================================================================
# check_celltype_identity.py
# ----------------------------------------------------------------------------
# Does each pseudobulk BAM actually contain the cell type its filename claims?
#
# WHY THIS EXISTS. On 2026-08-24 the Tier 4 cell-identity check failed (2 of 5
# lineage TFs landed in the wrong cell type). The failure signature -- injury
# signal fine everywhere, cell identity wrong -- is EXACTLY what a swapped BAM
# or cells leaking between cell types would produce, so it had to be excluded
# before anything else could be believed.
#
# THE TEST. For each cell type, take peaks present in ITS peak set and in no
# other cell type's ("unique peaks"), then count reads from every BAM in every
# cell type's unique peaks, normalised per million mapped reads. If the BAMs
# are correctly labelled the diagonal dominates.
#
# Measured 2026-08-24 (300 unique peaks per cell type, injured arms):
#     BAM \ peaks       Astro   Oligo   Micro  NeurD    own/mean(other)
#     Astrocytes_I      362.4    44.8    54.1   36.5           8.03
#     Oligodendro_I      50.7   269.2    54.3   33.8           5.82
#     Microglia_I        42.9    43.2   445.8   36.3          10.93
#     Neurons_D_I        56.5    52.3    58.6  249.5           4.47
# 4.5-11x diagonal dominance => BAMs are correctly labelled, and the Tier 4
# failure is NOT a data-handling bug. See 2026-08-24.txt section 5.
#
# Requires bedtools on PATH (env `bedtools` has v2.31.1).
#
# Usage from the repo root:
#   python pipeline/check_celltype_identity.py [--n 300] [--out DIR] CT [CT ...]
# ============================================================================
import argparse, os, subprocess, sys, tempfile
import numpy as np
import pysam

PEAKDIR = "meta/peaks_tobias"


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, text=True, **kw)


def unique_peaks(ct, others, workdir, n):
    """Peaks in ct's set that overlap no other cell type's peaks; subsampled to n."""
    others_cat = os.path.join(workdir, f"others_{ct}.bed")
    with open(others_cat, "w") as fh:
        for o in others:
            sh(f"cut -f1-3 {PEAKDIR}/{o}_peaks_no_blacklist_chr.bed", stdout=fh)
    merged = os.path.join(workdir, f"om_{ct}.bed")
    sh(f"sort -k1,1 -k2,2n {others_cat} | bedtools merge -i - > {merged}")
    uniq = os.path.join(workdir, f"uniqall_{ct}.bed")
    sh(f"cut -f1-3 {PEAKDIR}/{ct}_peaks_no_blacklist_chr.bed | sort -k1,1 -k2,2n "
       f"| bedtools intersect -v -a - -b {merged} > {uniq}")
    rows = [l.split()[:3] for l in open(uniq)]
    if not rows:
        sys.exit(f"{ct}: no unique peaks -- cannot run the identity test")
    step = max(1, len(rows) // n)
    picked = rows[::step][:n]
    return [(c, int(a), int(b)) for c, a, b in picked], len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("celltypes", nargs="+", help="cell types, e.g. Astrocytes Microglia")
    ap.add_argument("--arm", default="_I", help="BAM suffix (default _I)")
    ap.add_argument("--n", type=int, default=300, help="unique peaks sampled per cell type")
    ap.add_argument("--out", default=None, help="keep intermediate BEDs here")
    args = ap.parse_args()

    cts = args.celltypes
    workdir = args.out or tempfile.mkdtemp(prefix="idcheck_")
    os.makedirs(workdir, exist_ok=True)

    R = {}
    for ct in cts:
        R[ct], total = unique_peaks(ct, [o for o in cts if o != ct], workdir, args.n)
        print(f"{ct:<20} unique peaks {total:>7,}   sampled {len(R[ct])}", file=sys.stderr)

    print(f"\nreads per million, in each cell type's UNIQUE peaks ({args.n} each)")
    header = "BAM \\ peaks"
    print(f"{header:<20}" + "".join(f"{c[:13]:>14}" for c in cts) + f"{'own/mean(oth)':>15}")
    ok = True
    for bam_ct in cts:
        path = f"bam/{bam_ct}{args.arm}.bam"
        if not os.path.exists(path):
            sys.exit(f"missing {path}")
        bf = pysam.AlignmentFile(path, "rb")
        tot = sum(s.mapped for s in bf.get_index_statistics())
        vals = [sum(bf.count(c, a, b) for c, a, b in R[pk]) / tot * 1e6 for pk in cts]
        bf.close()
        own = vals[cts.index(bam_ct)]
        oth = np.mean([v for i, v in enumerate(vals) if cts[i] != bam_ct])
        ratio = own / oth if oth else float("inf")
        ok &= own == max(vals)
        print(f"{bam_ct + args.arm:<20}" + "".join(f"{v:>14.1f}" for v in vals) + f"{ratio:>15.2f}")

    print(f"\ndiagonal dominant for every cell type: {'YES' if ok else 'NO -- INVESTIGATE'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
