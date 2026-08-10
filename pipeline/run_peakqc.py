#!/usr/bin/env python
# ============================================================================
# run_peakqc.py
# ----------------------------------------------------------------------------
# PEAKQC fragment-length-distribution (FLD) scoring for every cell, run
# PER SAMPLE and then concatenated.
#
# WHY PER SAMPLE. PEAKQC scores cells individually, so each fragment must map to
# exactly one cell. The pooled per-cell-type BEDs in the chrombpnet repo cannot
# do that -- SplitFragments wrote bare 10x barcodes with no sample prefix, so
# cells from different samples sharing a barcode get merged (measured: Microglia
# 17,151 unique barcodes for 17,341 cells). Barcodes ARE unique within a sample
# -- verified, 0 duplicate (sample, barcode) pairs across all 67,072 cells --
# so each sample is scored against its own atac_fragments.tsv.gz.
#
# WHY IT MATTERS FOR TOBIAS. Footprint depth is a signal-to-noise problem, not a
# read-count problem. Cells with poor nucleosomal periodicity contribute
# essentially uniform background cut sites, which raise coverage while flattening
# the footprint. FLD scoring is a different axis from the standard QC the authors
# already applied (nFrags / TSS / doublets), so it is additive, not a redo.
#
# Only the 67,072 cells in the authors' final object are scored; barcodes not in
# meta/cell_metadata.csv are ignored.
#
# Usage (from the tobias repo root):
#   python pipeline/run_peakqc.py --pilot                 # smallest sample only
#   python pipeline/run_peakqc.py --samples FT_U_4 ALL_U_3
#   python pipeline/run_peakqc.py                         # all 21
#
# Outputs:
#   qc/peakqc_<sample>.csv     per-sample per-cell metrics (checkpoint)
#   meta/peakqc_scores.csv     all samples concatenated, keyed by `cell`
#   qc/plots/<sample>_*.png    density + overview plots
# ============================================================================
import argparse
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd

from peakqc.fld_scoring import add_fld_metrics

CHROMBPNET = "../TF_cobinding_multiome_chrombpnet"
METADATA = "meta/cell_metadata.csv"
QC_DIR = "qc"
PLOT_DIR = "qc/plots"
OUT = "meta/peakqc_scores.csv"


def log(*a):
    print(time.strftime("%H:%M:%S"), "|", *a, flush=True)


def fragments_path(sample):
    return os.path.join(CHROMBPNET, "data", sample, "outs", "atac_fragments.tsv.gz")


def run_sample(sample, md, n_threads, plot):
    """Score one sample's cells. Returns a DataFrame indexed by `cell`."""
    sub = md[md["sample"] == sample]
    frag = fragments_path(sample)
    if not os.path.exists(frag):
        log(f"  {sample}: MISSING {frag} -- skipping")
        return None
    if sub.empty:
        log(f"  {sample}: no cells in metadata -- skipping")
        return None

    # obs index must be the bare barcode as it appears in column 4 of the
    # fragments file; `cell` (sample-prefixed) is kept as a column so results
    # can be joined back unambiguously.
    obs = pd.DataFrame(
        {"cell": sub["cell"].values, "celltype": sub["celltype"].values},
        index=pd.Index(sub["barcode"].values, name="barcode"),
    )
    adata = ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs)

    before = set(adata.obs.columns)
    log(f"  {sample}: {len(obs):,} cells, scoring against {os.path.basename(frag)}")
    t0 = time.time()
    add_fld_metrics(
        adata,
        fragments=frag,
        barcode_col=None,          # use the obs index
        n_threads=n_threads,
        plot=plot,
        save_density=os.path.join(PLOT_DIR, f"{sample}_density.png") if plot else None,
        save_overview=os.path.join(PLOT_DIR, f"{sample}_overview.png") if plot else None,
    )
    added = [c for c in adata.obs.columns if c not in before]
    log(f"  {sample}: done in {time.time()-t0:.0f}s; added columns: {added}")

    out = adata.obs.copy()
    out["sample"] = sample
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*", default=None, help="subset of samples")
    ap.add_argument("--pilot", action="store_true", help="smallest sample only")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if not os.path.exists(METADATA):
        sys.exit(f"missing {METADATA} -- run pipeline/extract_cell_metadata.R first")
    md = pd.read_csv(METADATA, usecols=["cell", "barcode", "sample", "celltype"])
    log(f"metadata: {len(md):,} cells, {md['sample'].nunique()} samples")

    counts = md["sample"].value_counts()
    if args.pilot:
        samples = [counts.idxmin()]
        log(f"PILOT: smallest sample = {samples[0]} ({counts.min():,} cells)")
    elif args.samples:
        samples = args.samples
    else:
        samples = list(counts.sort_values().index)   # small first, fail fast

    os.makedirs(QC_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)

    frames, failed = [], []
    for i, s in enumerate(samples, 1):
        log(f"[{i}/{len(samples)}] {s}")
        ckpt = os.path.join(QC_DIR, f"peakqc_{s}.csv")
        if os.path.exists(ckpt):
            log(f"  {s}: checkpoint exists, loading")
            frames.append(pd.read_csv(ckpt))
            continue
        try:
            df = run_sample(s, md, args.threads, not args.no_plot)
        except Exception as e:                      # one bad sample must not kill the run
            log(f"  {s}: FAILED -- {type(e).__name__}: {e}")
            failed.append(s)
            continue
        if df is None:
            failed.append(s)
            continue
        df.to_csv(ckpt, index=True)
        frames.append(df.reset_index())

    if not frames:
        sys.exit("no samples produced results")

    allres = pd.concat(frames, ignore_index=True)
    os.makedirs("meta", exist_ok=True)
    allres.to_csv(OUT, index=False)
    log(f"wrote {OUT}: {len(allres):,} cells from {len(frames)} samples")
    if failed:
        log(f"FAILED samples ({len(failed)}): {failed}")
        sys.exit(1)
    log("DONE")


if __name__ == "__main__":
    main()
