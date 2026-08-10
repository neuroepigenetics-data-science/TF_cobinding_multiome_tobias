#!/usr/bin/env python
# ============================================================================
# analyze_peakqc.py
# ----------------------------------------------------------------------------
# Turn per-cell FLD scores into the decision they exist to inform: which
# cell types can support TOBIAS footprinting, and at what threshold.
#
# The point is EFFECTIVE depth. The raw per-cell-type fragment counts (Microglia
# 284M, Ependymal 101M, ...) assume every cell contributes usable signal. Cells
# with no nucleosomal periodicity contribute near-uniform background cut sites,
# which raise coverage while flattening the footprint. So the number that
# actually matters is fragments contributed by cells that PASS, not total
# fragments. This script computes both and shows the gap.
#
# It also picks a threshold from the data rather than applying PEAKQC's
# rule-of-thumb 100 blindly: the pilot showed a sharply bimodal distribution
# with an empty band between ~1 and ~50, so the natural cut is the gap, and
# where that gap sits is a property of this dataset.
#
# Usage (from the tobias repo root):
#   python pipeline/analyze_peakqc.py [--threshold N]
#
# Reads:  meta/peakqc_scores.csv, meta/cell_metadata.csv
# Writes: qc/peakqc_summary.txt, qc/celltype_viability.csv
# ============================================================================
import argparse
import os

import numpy as np
import pandas as pd

SCORES = "meta/peakqc_scores.csv"
METADATA = "meta/cell_metadata.csv"


def find_gap(scores, lo=0.1, hi=500, nbins=200):
    """Largest empty band in log-space -- the natural cut in a bimodal set.

    Returns (gap_low, gap_high, suggested_threshold) or None if unimodal.
    """
    s = scores[scores > 0]
    if len(s) < 50:
        return None
    edges = np.logspace(np.log10(lo), np.log10(hi), nbins)
    counts, _ = np.histogram(s, bins=edges)
    best_len = best_i = 0
    cur_len = cur_i = 0
    for i, c in enumerate(counts):
        if c == 0:
            if cur_len == 0:
                cur_i = i
            cur_len += 1
            if cur_len > best_len:
                best_len, best_i = cur_len, cur_i
        else:
            cur_len = 0
    if best_len < 5:                       # no meaningful gap -> unimodal
        return None
    g_lo, g_hi = edges[best_i], edges[best_i + best_len]
    return g_lo, g_hi, float(np.sqrt(g_lo * g_hi))   # geometric midpoint


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=None,
                    help="override; default is the data-driven gap midpoint")
    args = ap.parse_args()

    if not os.path.exists(SCORES):
        raise SystemExit(f"missing {SCORES} -- run pipeline/run_peakqc.py first")
    d = pd.read_csv(SCORES)
    md = pd.read_csv(METADATA, usecols=["cell", "experiment", "condition", "cluster_ids_injury"])
    d = d.merge(md, on="cell", how="left")
    out = []

    def emit(s=""):
        print(s)
        out.append(s)

    emit(f"PEAKQC summary — {len(d):,} cells, {d['sample'].nunique()} samples")
    emit("=" * 78)

    # --- threshold ----------------------------------------------------------
    gap = find_gap(d["fld_score"].values)
    if args.threshold is not None:
        thr = args.threshold
        emit(f"\nthreshold: {thr} (user-specified)")
    elif gap:
        thr = gap[2]
        emit(f"\ndistribution is BIMODAL: empty band {gap[0]:.2f} – {gap[1]:.1f}")
        emit(f"threshold: {thr:.1f} (geometric midpoint of the gap)")
        emit(f"  for comparison, PEAKQC's rule of thumb is 100")
    else:
        thr = 100.0
        emit("\ndistribution is NOT clearly bimodal; falling back to PEAKQC's 100")
    d["pass"] = d["fld_score"] >= thr
    emit(f"\noverall pass: {d['pass'].sum():,} / {len(d):,} ({100*d['pass'].mean():.1f}%)")

    # --- is it an artefact of depth? ---------------------------------------
    emit("\n" + "-" * 78)
    emit("Is the split explained by depth?  (if yes, it is not a quality signal)")
    for lbl, g in [("fail", d[~d["pass"]]), ("pass", d[d["pass"]])]:
        emit(f"  {lbl}: n_fragments median={g['n_fragments'].median():>9,.0f}   "
             f"mean_fragment_size median={g['mean_fragment_size'].median():.1f}")
    rho = d[["n_fragments", "fld_score"]].corr(method="spearman").iloc[0, 1]
    emit(f"  spearman(n_fragments, fld_score) = {rho:.3f}"
         f"   {'-> largely independent of depth' if abs(rho) < 0.3 else '-> DEPTH-CONFOUNDED, be careful'}")

    # --- per experiment: is the bimodality sample-specific? -----------------
    emit("\n" + "-" * 78)
    emit("Pass rate by experiment / condition  (is this specific to sorted samples?)")
    for key in ["experiment", "condition"]:
        if key in d.columns and d[key].notna().any():
            t = d.groupby(key)["pass"].agg(["count", "sum", "mean"])
            t.columns = ["cells", "passing", "frac"]
            emit(f"\n  by {key}:")
            for k, r in t.iterrows():
                emit(f"    {str(k):<14} {int(r.cells):>7,} cells   {100*r.frac:>5.1f}% pass")

    # --- the payload: effective depth per cell type -------------------------
    emit("\n" + "=" * 78)
    emit("EFFECTIVE DEPTH per cell type  (fragments from PASSING cells only)")
    emit("=" * 78)
    rows = []
    for ct, g in d.groupby("celltype"):
        raw = g["n_fragments"].sum()
        eff = g.loc[g["pass"], "n_fragments"].sum()
        rows.append(dict(celltype=ct, cells=len(g), passing=int(g["pass"].sum()),
                         pass_frac=g["pass"].mean(), raw_fragments=int(raw),
                         effective_fragments=int(eff),
                         retained=eff / raw if raw else 0,
                         median_score=g["fld_score"].median()))
    t = pd.DataFrame(rows).sort_values("effective_fragments", ascending=False)
    emit(f"\n{'celltype':<18}{'cells':>7}{'pass':>7}{'pass%':>7}"
         f"{'raw frags':>14}{'effective':>14}{'kept':>7}")
    for _, r in t.iterrows():
        emit(f"{r.celltype:<18}{r.cells:>7,}{r.passing:>7,}{100*r.pass_frac:>6.0f}%"
             f"{r.raw_fragments:>14,}{r.effective_fragments:>14,}{100*r.retained:>6.0f}%")

    # --- and per cell type x injury, which is what BINDetect compares -------
    if "cluster_ids_injury" in d.columns and d["cluster_ids_injury"].notna().any():
        emit("\n" + "=" * 78)
        emit("EFFECTIVE DEPTH per cell type x injury  (BINDetect's actual comparison)")
        emit("=" * 78)
        rows = []
        for grp, g in d.groupby("cluster_ids_injury"):
            eff = g.loc[g["pass"], "n_fragments"].sum()
            rows.append(dict(group=grp, cells=len(g), effective=int(eff),
                             raw=int(g["n_fragments"].sum())))
        t2 = pd.DataFrame(rows).sort_values("effective", ascending=False)
        emit(f"\n{'group':<22}{'cells':>8}{'raw frags':>14}{'effective':>14}{'kept':>7}")
        for _, r in t2.iterrows():
            k = 100 * r.effective / r.raw if r.raw else 0
            emit(f"{r.group:<22}{r.cells:>8,}{r.raw:>14,}{r.effective:>14,}{k:>6.0f}%")

    os.makedirs("qc", exist_ok=True)
    t.to_csv("qc/celltype_viability.csv", index=False)
    with open("qc/peakqc_summary.txt", "w") as fh:
        fh.write("\n".join(out) + "\n")
    emit("\nwrote qc/peakqc_summary.txt and qc/celltype_viability.csv")


if __name__ == "__main__":
    main()
