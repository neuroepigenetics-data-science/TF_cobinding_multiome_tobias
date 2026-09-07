#!/usr/bin/env python
# ============================================================================
# diagnose_tier3.py -- POST-HOC diagnostics for the Tier 3 result
# ----------------------------------------------------------------------------
# *** THIS DOES NOT RESCORE TIER 3. ***
# The pre-registered verdict is computed and written by run_tfcomb_tier3.py and
# stands exactly as recorded. Everything here is a post-hoc diagnostic asking
# WHY the matrix came out as it did, and every number below is labelled as such.
# Per VALIDATION_CRITERIA.md's own rule, a criterion already judged against
# results is never rewritten -- a critique is appended instead.
#
# THE OBSERVATION THAT PROMPTED THIS. In the scored matrix, NFI+SOX9 is the best
# column in 3 of 4 rows and the two NFI-containing columns are the top two
# almost everywhere. That is a COLUMN effect, not the row effect the tier tests.
#
# THE STRUCTURAL SUSPICION. The pre-registered partner sets pool families of
# very unequal size, and the statistic is a median over the POOLED pairs:
#
#     NFI+SOX9    = 4 NFI motifs + 1 SOX9   -> 140 + 35 pairs   (80% NFI)
#     NFI+RFX     = 4 NFI + 4 RFX           -> 140 + 140 pairs  (50% NFI)
#     SOX10+TCF4  = 1 + 1                   ->  35 +  35 pairs
#     ELK3+CEBPA  = 1 + 1                   ->  35 +  35 pairs
#
# So "NFI+SOX9" is four-fifths an NFI measurement, and a column can win on the
# strength of one promiscuous family while the lineage-specific member of the
# pair contributes almost nothing. That is the same class of defect --
# a statistic responding to family size rather than biology -- that demoted
# Tier 1; see the motif-family-size-bias note. Decomposing per family is the
# direct test of it.
#
# THE SECOND QUESTION: SPECIFICITY. Diagonal dominance is a WITHIN-row test by
# design. The complementary view is cross-row: does a lineage family rank best
# in ITS OWN cell type? That is only legitimate among the three DIFFERENTIAL
# rows -- Ependymal_I is a single-arm occupancy quantity and is excluded from
# every cross-row comparison here, exactly as the criterion's limitation 4 says.
#
# Usage: ~/venvs/tfcomb/bin/python pipeline/diagnose_tier3.py [--dir tobias/tfcomb]
# ============================================================================
import argparse
import json
import os

import numpy as np
import pandas as pd

DIFFERENTIAL = ["Astrocytes", "Oligodendrocytes", "Microglia"]   # cross-row safe
SINGLE_ARM = ["Ependymal_I"]
LINEAGE = ["NFI", "SOX9", "RFX", "SOX10", "TCF4", "ELK3", "CEBPA"]
# Which cell type each lineage family is the paper's expectation for.
EXPECTED_FOR = {"SOX9": "Astrocytes", "NFI": "Astrocytes",
                "RFX": "Ependymal_I", "SOX10": "Oligodendrocytes",
                "TCF4": "Oligodendrocytes", "ELK3": "Microglia",
                "CEBPA": "Microglia"}


def med_rank(M, idx, a, b):
    ai = [idx[t] for t in a if t in idx]
    bi = [idx[t] for t in b if t in idx]
    if not ai or not bi:
        return np.nan, 0
    v = M[np.ix_(ai, bi)].ravel()
    v = v[~np.isnan(v)]
    return (float(np.median(v)) if len(v) else np.nan), int(len(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="tobias/tfcomb")
    ap.add_argument("--dist-root", default="tobias/bindetect",
                    help="BINDetect output root, for the <ct>_distances.txt matrices")
    args = ap.parse_args()

    cache = os.path.join(args.dir, "rank_matrices")
    meta = json.load(open(os.path.join(cache, "index.json")))
    tf_names, fams = meta["tf_names"], meta["families"]
    total = meta["total_pairs"]
    idx = {t: i for i, t in enumerate(tf_names)}
    ap1 = fams["AP-1"]
    cts = [c for c in DIFFERENTIAL + SINGLE_ARM]

    mats = {ct: np.load(os.path.join(cache, f"{ct}.npy")) for ct in cts}

    # ---- 1. per-family decomposition ------------------------------------
    rows = {}
    for ct in cts:
        rows[ct] = {}
        for fam in LINEAGE:
            m, n = med_rank(mats[ct], idx, ap1, fams[fam])
            rows[ct][fam] = m
    F = pd.DataFrame(rows).T[LINEAGE]

    print("=" * 82)
    print("DIAGNOSTIC 1 -- median rank of AP-1 <-> FAMILY pairs, family by family")
    print("=" * 82)
    print("(post-hoc; the scored Tier 3 matrix pooled these into 4 columns)")
    print()
    disp = F.copy()
    disp.insert(0, "pairs_total", [total[c] for c in F.index])
    print(disp.to_string(float_format=lambda v: f"{v:,.0f}"))
    print()
    print("as PERCENTILE of all pairs in that cell type (lower = stronger):")
    P = F.div(pd.Series({c: total[c] for c in F.index}), axis=0) * 100
    print(P.to_string(float_format=lambda v: f"{v:.1f}%"))

    # ---- 2. is the NFI+SOX9 column just NFI? ----------------------------
    print()
    print("=" * 82)
    print("DIAGNOSTIC 2 -- does the winning column reduce to NFI alone?")
    print("=" * 82)
    for ct in cts:
        nfi, _ = med_rank(mats[ct], idx, ap1, fams["NFI"])
        sox9, _ = med_rank(mats[ct], idx, ap1, fams["SOX9"])
        pooled, _ = med_rank(mats[ct], idx, ap1, fams["NFI"] + fams["SOX9"])
        print(f"  {ct:17s} NFI={nfi:9,.0f}  SOX9={sox9:9,.0f}  pooled(NFI+SOX9)={pooled:9,.0f}"
              f"   pooled-NFI gap={pooled - nfi:+,.0f}")
    print("  -> a small gap means the pooled column is carrying NFI's signal, not SOX9's.")

    # ---- 3. specificity: does a family rank best in its own cell type? --
    print()
    print("=" * 82)
    print("DIAGNOSTIC 3 -- cross-row specificity, DIFFERENTIAL rows only")
    print("=" * 82)
    print("(Ependymal_I excluded: single-arm occupancy is a different quantity")
    print(" and is never comparable across rows -- criterion limitation 4)")
    print()
    Fd = F.loc[DIFFERENTIAL]
    print(Fd.to_string(float_format=lambda v: f"{v:,.0f}"))
    print()
    hits = []
    for fam in LINEAGE:
        exp = EXPECTED_FOR[fam]
        if exp not in DIFFERENTIAL:
            print(f"  {fam:6s} expected in {exp:17s} -- not testable across rows (single-arm)")
            continue
        best = Fd[fam].idxmin()
        ok = best == exp
        hits.append(ok)
        print(f"  {'OK  ' if ok else 'MISS'} {fam:6s} expected strongest in {exp:17s}"
              f" actually strongest in {best}")
    print(f"\n  {sum(hits)}/{len(hits)} lineage families rank strongest in their own cell type.")
    print("  NOTE: this is a POST-HOC view, not a pre-registered criterion. It is")
    print("  reported because it is the question Tier 3 was trying to ask, and")
    print("  because a 'nothing is specific' answer and a 'the statistic was wrong'")
    print("  answer have different consequences for the pipeline.")

    # ---- 4. the null, applied per INDIVIDUAL family ----------------------
    # Same instrument as the pre-registered secondary criterion, but without the
    # pooling. This matters because a random motif paired with AP-1 already sits
    # near the top of the table -- AP-1 pairs are globally high-ranked, which is
    # the co-binding face of "AP-1 binds basically everywhere". Raw percentiles
    # therefore flatter every family, and only the null says what is real.
    print()
    print("=" * 82)
    print("DIAGNOSTIC 4 -- per-family permutation null (10,000 draws, size-matched)")
    print("=" * 82)
    rng = np.random.default_rng(1)
    ap1_set = set(ap1)
    pool = [t for t in tf_names if t not in ap1_set]
    draws_n = 10000
    recs = []
    for ct in cts:
        M_ct = mats[ct]
        ap1_i = [idx[t] for t in ap1 if t in idx]
        pool_i = np.array([idx[t] for t in pool if t in idx])
        sub = M_ct[ap1_i, :]
        for fam in LINEAGE:
            obs, _ = med_rank(M_ct, idx, ap1, fams[fam])
            n = len(fams[fam])
            d = np.empty(draws_n)
            for i in range(draws_n):
                pick = rng.choice(pool_i, size=n, replace=False)
                v = sub[:, pick].ravel()
                v = v[~np.isnan(v)]
                d[i] = np.median(v) if len(v) else np.nan
            pval = float((d <= obs).sum() + 1) / (draws_n + 1)
            recs.append({"cell_type": ct, "family": fam, "n_motifs": n,
                         "observed": obs, "null_median": float(np.nanmedian(d)),
                         "p": pval})
    R = pd.DataFrame(recs)
    piv = R.pivot(index="cell_type", columns="family", values="p")[LINEAGE]
    piv = piv.loc[[c for c in cts]]
    print("p-values (one-sided; p<0.05 = stronger than a size-matched random partner)")
    print(piv.to_string(float_format=lambda v: f"{v:.4f}"))
    print()
    nfam, ncell = len(LINEAGE), len(cts)
    bonf = 0.05 / (nfam * ncell)
    print(f"Bonferroni across all {nfam*ncell} family x cell-type tests: p < {bonf:.5f}")
    sig = R[R["p"] < bonf].sort_values("p")
    if len(sig):
        print("survives Bonferroni:")
        for _, r in sig.iterrows():
            tag = " <- paper's expectation" if EXPECTED_FOR.get(r["family"]) == r["cell_type"] else ""
            print(f"    {r['cell_type']:17s} {r['family']:6s} obs={r['observed']:9,.0f} "
                  f"null={r['null_median']:9,.0f} p={r['p']:.5f}{tag}")
    else:
        print("    nothing survives Bonferroni.")
    R.to_csv(os.path.join(args.dir, "tier3_diagnostic_family_null.csv"), index=False)

    # ---- 5. what ARE the top AP-1 partners? ------------------------------
    # Asked without reference to any expectation: if the pipeline sees real
    # co-binding structure, the top partners should be interpretable.
    print()
    print("=" * 82)
    print("DIAGNOSTIC 5 -- strongest AP-1 partners per cell type (no expectation used)")
    print("=" * 82)
    for ct in cts:
        M_ct = mats[ct]
        ap1_i = [idx[t] for t in ap1 if t in idx]
        sub = M_ct[ap1_i, :]
        with np.errstate(invalid="ignore"):
            best = np.nanmin(sub, axis=0)          # best rank each partner achieves
        order = np.argsort(np.where(np.isnan(best), np.inf, best))
        names = []
        for j in order:
            t = tf_names[j]
            if t in ap1_set:
                continue
            names.append(f"{t.split('_')[0]}({best[j]:,.0f})")
            if len(names) == 12:
                break
        print(f"  {ct:17s} {', '.join(names)}")

    # ---- 6. are the top partners just AP-1 look-alikes? ------------------
    # BINDetect writes a 746x746 motif DISTANCE matrix per cell type
    # (<ct>_distances.txt; 0 = identical motif, ~1 = unrelated). If the top
    # AP-1 partners sit at small distances from AP-1 itself, the co-occurrence
    # is motif redundancy -- the same DNA matched twice -- and not co-binding.
    # This is divergence 3 ("footprinting cannot resolve near-identical
    # motifs") measured directly rather than assumed.
    dist_path = os.path.join(args.dist_root, "Astrocytes", "Astrocytes_distances.txt")
    if os.path.exists(dist_path):
        print()
        print("=" * 82)
        print("DIAGNOSTIC 6 -- motif distance of the top partners from AP-1 itself")
        print("=" * 82)
        dnames = [n for n in pd.read_csv(dist_path, sep="\t", nrows=0).columns if n != "#"]
        raw = np.loadtxt(dist_path, skiprows=1)
        DM = pd.DataFrame(raw, index=dnames, columns=dnames)
        core = [c for c in dnames if c.split("_")[0] in ("FOSJUN", "FOSL1", "FOS", "JUNB", "JUN")]
        print(f"(distance to the AP-1 core motifs {[c.split('_')[0] for c in core]};")
        print(" lineage partners for reference: SOX9/SOX10/NFI/RFX/ELK3/TCF4 all sit at ~0.997-0.999)")
        print()
        for ct in cts:
            M_ct = mats[ct]
            ap1_i = [idx[t] for t in ap1 if t in idx]
            sub = M_ct[ap1_i, :]
            with np.errstate(invalid="ignore"):
                best = np.nanmin(sub, axis=0)
            order = np.argsort(np.where(np.isnan(best), np.inf, best))
            top, ds = [], []
            for j in order:
                t = tf_names[j]
                if t in ap1_set or t not in DM.index:
                    continue
                d = float(DM.loc[t, core].astype(float).min())
                top.append(f"{t.split('_')[0]}({d:.2f})")
                ds.append(d)
                if len(top) == 10:
                    break
            print(f"  {ct:17s} median dist={np.median(ds):.3f}   {', '.join(top)}")
        print()
        print("  A LOW median means that cell type's top AP-1 partners are largely")
        print("  motifs resembling AP-1 -- an artifact. A HIGH median means the top")
        print("  partners are genuinely distinct factors.")

    F.to_csv(os.path.join(args.dir, "tier3_diagnostic_family_ranks.csv"))
    print(f"\nwrote {args.dir}/tier3_diagnostic_family_ranks.csv, tier3_diagnostic_family_null.csv")


if __name__ == "__main__":
    main()
