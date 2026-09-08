#!/usr/bin/env python
# ============================================================================
# run_tfcomb_tier3.py -- TIER 3, the load-bearing discrimination test
# ----------------------------------------------------------------------------
# Scored EXACTLY as pre-registered in VALIDATION_CRITERIA.md ("TIER 3 -- THE
# LOAD-BEARING CHECK", rewritten 2026-08-24 before TF-COMB was installed, so it
# is genuine pre-registration). Nothing here was tuned after seeing a number;
# any deviation from the written criterion is listed under DIVERGENCES below and
# was decided before the first run.
#
# *** THE QUESTION ***
# The paper's mechanism is AP-1 PLUS A LINEAGE-SPECIFIC PARTNER -- the partner
# is what makes an injury-responsive enhancer cell-type specific. Tier 1 showed
# empirically that AP-1's own rank does not discriminate cell types (it is rank 1
# nearly everywhere, neurons included). So cell-type specificity shows up in the
# co-binding or nowhere, which is why this tier is load-bearing and the others
# are not.
#
# *** THE TEST: DIAGONAL DOMINANCE OF A 4x4 MATRIX ***
# Rows = cell types, columns = each cell type's expected partner set.
#
#                    NFI+SOX9   NFI+RFX   SOX10+TCF4   ELK3+CEBPA
#   Astrocytes          *
#   Ependymal                       *
#   Oligodendrocytes                             *
#   Microglia                                              *
#
# Score(cell type c, partner set P) = the MEDIAN RANK, within c's own rule
# table, of every AP-1<->P pair rule.  Lower rank = stronger co-occurrence.
# PASS = the diagonal entry is the best (lowest median rank) in its ROW for
# >= 3 of the 4 cell types.
#
# *** THE NAMED FAILURE MODE ***
# The pre-registration names this in advance and it must not be reported as a
# pass: if EVERY partner set ranks highly in EVERY cell type, that is a
# NON-DISCRIMINATING result and counts as a FAILURE, even though every expected
# partner "appears". It is the exact pathology Tier 1 turned out to have.
#
# *** DIVERGENCES FROM THE WRITTEN CRITERION -- all decided before running ***
# 1. AP-1 HAS 35 MEMBERS, NOT 33. The criterion says "AP-1 = the 33 matching
#    motifs" and pipeline/plot_differential_heatmap.py repeats the figure, but
#    the pre-registered regex (FOS|JUN|ATF3|BATF|JDP2) matches 35 of the 746,
#    and all 35 are genuine AP-1/TRE-binding bZIPs -- no false positives. The
#    "33" was a miscount, not a different family definition. The DEFINITION is
#    what was pre-registered, so all 35 are used and the slip is recorded here.
# 2. RFX IS HELD TO RFX1/2/3/4 as the criterion states, although RFX5 and RFX7
#    are also in the motif set. RFX5 (MHC-II X-box) and RFX7 are distinct
#    factors; including them would be a post-hoc widening of a pre-registered
#    set, so they are excluded even though it is the less flattering choice.
# 3. THE DIFFERENTIAL ROWS RANK BY log2fc OF COSINE, the Ependymal row by raw
#    cosine. This follows the criterion's own limitation 4, which says the three
#    two-armed rows use a DiffCombObj and Ependymal a single-condition CombObj.
#    The quantities differ, which is why diagonal dominance is a WITHIN-ROW
#    comparison; it is never compared across rows.
# 4. RULES ARE DEDUPLICATED TO UNORDERED PAIRS. TF-COMB emits both A->B and
#    B->A with identical cosine; keeping both would double-count every pair in
#    the median. This is arithmetic hygiene, not a scoring choice.
#
# *** EPENDYMAL IS KEPT, AND WHY THAT IS NOT A CONTRADICTION ***
# Supervisor ruling 6 (2026-08-26) dropped ependymal cells. That ruling closes
# the uninjured-material question and forbids any DIFFERENTIAL claim; it
# explicitly notes the occupancy results "can still feed a single-condition
# TF-COMB run". Tier 3's limitation 4 anticipated exactly this and declared the
# within-row test legitimate. Dropping the row would also remove RFX, which
# limitation 3 names as one of the three CLEANEST discriminators. So the row is
# computed and reported, carrying its caveat. The 3-of-4 threshold is unchanged.
#
# *** macOS ***
# fork start method, as for TOBIAS -- see pipeline/tobias_fork.py. Must be run
# as a file, never piped from stdin, or the workers re-import __main__ and fork
# endlessly. See pipeline/setup_tfcomb_env.sh.
#
# Usage:
#   ~/venvs/tfcomb/bin/python pipeline/run_tfcomb_tier3.py \
#       [--input tobias/tfcomb/input] [--out tobias/tfcomb] [--threads 8] \
#       [--draws 10000] [--seed 1]
# ============================================================================
import argparse
import json
import multiprocessing as mp
import os
import re
import sys
import time

# ---------------------------------------------------------------------------
# Pre-registered family definitions. AP1_RE is copied verbatim from
# pipeline/plot_differential_heatmap.py so the two analyses cannot drift apart.
# ---------------------------------------------------------------------------
AP1_RE = re.compile(r"FOS|JUN|ATF3|BATF|JDP2", re.I)   # JDP2 binds the AP-1/TRE site

FAMILY_RE = {
    "NFI":   re.compile(r"^NFI[ABCX]_"),        # NFIA/B/C/X, per the criterion
    "RFX":   re.compile(r"^RFX[1234]_"),        # RFX1/2/3/4 only -- divergence 2
    "SOX9":  re.compile(r"^SOX9_"),
    "SOX10": re.compile(r"^SOX10_"),
    "TCF4":  re.compile(r"^TCF4_"),
    "ELK3":  re.compile(r"^ELK3_"),
    "CEBPA": re.compile(r"^CEBPA_"),
}

# Columns of the 4x4 matrix: the paper's expected partner set per cell type.
PARTNER_SETS = {
    "NFI+SOX9":   ["NFI", "SOX9"],
    "NFI+RFX":    ["NFI", "RFX"],
    "SOX10+TCF4": ["SOX10", "TCF4"],
    "ELK3+CEBPA": ["ELK3", "CEBPA"],
}

# Rows: cell type -> (its diagonal column, the arms it is built from).
# A one-element arm tuple means a single-condition CombObj (no differential).
ROWS = [
    ("Astrocytes",       "NFI+SOX9",   ("Astrocytes_I", "Astrocytes_U")),
    ("Ependymal_I",      "NFI+RFX",    ("Ependymal_I",)),
    ("Oligodendrocytes", "SOX10+TCF4", ("Oligodendrocytes_I", "Oligodendrocytes_U")),
    ("Microglia",        "ELK3+CEBPA", ("Microglia_I", "Microglia_U")),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def resolve_families(tf_names):
    """Map each family label to its member motifs, present in this motif set."""
    fams = {"AP-1": sorted(t for t in tf_names if AP1_RE.search(t))}
    for label, rx in FAMILY_RE.items():
        fams[label] = sorted(t for t in tf_names if rx.match(t))
    return fams


def build_combobj(bed, prefix, threads, verbosity=0):
    """One arm: bound sites -> co-occurrence counts -> market basket rules."""
    from tfcomb import CombObj
    C = CombObj(verbosity=verbosity)
    C.TFBS_from_bed(bed)
    C.set_prefix(prefix)
    log(f"    {prefix}: {len(C.TFBS):,} sites")
    C.count_within(min_dist=0, max_dist=100, min_overlap=0, max_overlap=0,
                   threads=threads)
    C.market_basket(measure="cosine", threads=threads)
    log(f"    {prefix}: {C.rules.shape[0]:,} rules")
    return C


def ranked_table(rules, measure_col):
    """Rank rules by `measure_col` descending; rank 1 = strongest.

    Deduplicated to unordered pairs (divergence 4): TF-COMB emits A->B and B->A
    with identical measures, and keeping both double-counts every pair.
    """
    import pandas as pd
    df = rules.copy()
    key = [tuple(sorted(p)) for p in zip(df["TF1"], df["TF2"])]
    df["_pair"] = key
    df = df[df["TF1"] != df["TF2"]]                       # drop self-pairs
    df = df.drop_duplicates(subset="_pair", keep="first")
    df = df.sort_values(measure_col, ascending=False).reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)
    return df


def rank_matrix(df, tf_names):
    """Symmetric (n_tf x n_tf) matrix of pair ranks, NaN where a pair is absent.

    Built once per cell type so the 10,000-draw permutation null is a matrix
    slice per draw instead of a full-table scan; the naive .isin() version was
    40,000 scans of a ~280k-row table.
    """
    import numpy as np
    idx = {t: i for i, t in enumerate(tf_names)}
    n = len(tf_names)
    M = np.full((n, n), np.nan, dtype=np.float64)
    i = df["TF1"].map(idx).to_numpy()
    j = df["TF2"].map(idx).to_numpy()
    r = df["rank"].to_numpy(dtype=np.float64)
    ok = ~(np.isnan(i.astype(float)) | np.isnan(j.astype(float)))
    i, j, r = i[ok].astype(int), j[ok].astype(int), r[ok]
    M[i, j] = r
    M[j, i] = r                      # rules were deduped to unordered pairs
    return M, idx


def pair_ranks(M, idx, set_a, set_b):
    """Ranks of every pair joining a motif in set_a to one in set_b."""
    import numpy as np
    ai = [idx[t] for t in set_a if t in idx]
    bi = [idx[t] for t in set_b if t in idx]
    if not ai or not bi:
        return np.array([])
    sub = M[np.ix_(ai, bi)].ravel()
    return sub[~np.isnan(sub)]


def main():
    import numpy as np
    import pandas as pd

    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="tobias/tfcomb/input")
    ap.add_argument("--out", default="tobias/tfcomb")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--draws", type=int, default=10000,
                    help="permutation draws for the secondary null")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # -- build every arm once ------------------------------------------------
    arms = sorted({a for _, _, tup in ROWS for a in tup})
    log(f"building {len(arms)} arms: {', '.join(arms)}")
    objs = {}
    for arm in arms:
        bed = os.path.join(args.input, f"{arm}.bed")
        if not os.path.exists(bed):
            sys.exit(f"ERROR: missing input {bed} -- run pipeline/prep_tfcomb_input.sh")
        objs[arm] = build_combobj(bed, arm, args.threads)

    # -- families, resolved against the real motif names ---------------------
    tf_names = sorted({s.name for s in objs[arms[0]].TFBS})
    fams = resolve_families(tf_names)
    log(f"motifs in set: {len(tf_names)}")
    for label in ["AP-1"] + list(FAMILY_RE):
        log(f"    {label:6s} n={len(fams[label]):2d}  {', '.join(fams[label][:4])}"
            + (" ..." if len(fams[label]) > 4 else ""))
    if not fams["AP-1"]:
        sys.exit("ERROR: no AP-1 motifs resolved")
    ap1 = fams["AP-1"]
    ap1_set = set(ap1)
    # Null draws sample from non-AP-1 motifs only: an AP-1<->AP-1 pair is a
    # within-family pair and would not be a fair stand-in for a partner.
    pool = [t for t in tf_names if t not in ap1_set]

    # -- one ranked rule table per ROW --------------------------------------
    tables, measures, mats = {}, {}, {}
    for ct, _diag, tup in ROWS:
        if len(tup) == 2:
            inj, unj = tup
            log(f"{ct}: differential {inj} vs {unj}")
            D = objs[inj].compare(objs[unj], measure="cosine")
            rules = D.rules
            cand = [c for c in rules.columns if c.endswith("_log2fc")]
            if len(cand) != 1:
                sys.exit(f"ERROR: expected 1 log2fc column, got {cand}")
            col = cand[0]
        else:
            log(f"{ct}: single-arm occupancy (no differential claim)")
            rules = objs[tup[0]].rules
            col = "cosine"
        t = ranked_table(rules, col)
        tables[ct] = t
        mats[ct] = rank_matrix(t, tf_names)
        measures[ct] = col
        log(f"    ranked on '{col}', {len(t):,} unordered pairs")

    # -- the 4x4 matrix ------------------------------------------------------
    matrix, detail = {}, []
    for ct, diag, _tup in ROWS:
        M_ct, idx_ct = mats[ct]
        row = {}
        for colname, famlabels in PARTNER_SETS.items():
            members = [m for f in famlabels for m in fams[f]]
            r = pair_ranks(M_ct, idx_ct, ap1, members)
            row[colname] = float(np.median(r)) if len(r) else float("nan")
            detail.append({"cell_type": ct, "partner_set": colname,
                           "n_partner_motifs": len(members), "n_pairs": int(len(r)),
                           "median_rank": row[colname],
                           "is_diagonal": colname == diag})
        matrix[ct] = row

    M = pd.DataFrame(matrix).T[list(PARTNER_SETS)]
    total_pairs = {ct: len(tables[ct]) for ct in M.index}

    # -- primary criterion: diagonal dominance -------------------------------
    wins, verdicts = 0, []
    for ct, diag, _tup in ROWS:
        best = M.loc[ct].idxmin()
        ok = best == diag
        wins += ok
        verdicts.append((ct, diag, best, ok))

    # -- secondary criterion: permutation null on the diagonal ---------------
    log(f"permutation null: {args.draws:,} draws per row")
    nulls = {}
    for ct, diag, _tup in ROWS:
        M_ct, idx_ct = mats[ct]
        n = len([m for f in PARTNER_SETS[diag] for m in fams[f]])
        obs = M.loc[ct, diag]
        ap1_i = [idx_ct[t] for t in ap1 if t in idx_ct]
        pool_i = np.array([idx_ct[t] for t in pool if t in idx_ct])
        sub_ap1 = M_ct[ap1_i, :]                    # 35 x n_tf, sliced once
        draws = np.empty(args.draws)
        for i in range(args.draws):
            pick = rng.choice(pool_i, size=n, replace=False)
            v = sub_ap1[:, pick].ravel()
            v = v[~np.isnan(v)]
            draws[i] = np.median(v) if len(v) else np.nan
        # one-sided: how often does a random partner set rank at least as high?
        p = float((draws <= obs).sum() + 1) / (args.draws + 1)
        nulls[ct] = {"n_motifs": n, "observed_median_rank": float(obs),
                     "null_median": float(np.nanmedian(draws)),
                     "null_p2.5": float(np.nanpercentile(draws, 2.5)),
                     "p_value": p, "passes_p05": bool(p < 0.05)}
        log(f"    {ct:17s} obs={obs:9.1f}  null median={np.nanmedian(draws):9.1f}  p={p:.4f}")

    # -- report --------------------------------------------------------------
    print()
    print("=" * 78)
    print("TIER 3 -- median rank of AP-1 <-> partner-set pairs (lower = stronger)")
    print("=" * 78)
    disp = M.copy()
    disp.insert(0, "n_pairs_total", [total_pairs[i] for i in M.index])
    print(disp.to_string(float_format=lambda v: f"{v:,.1f}"))
    print()
    print("Diagonal = the partner set the paper expects for that cell type.")
    for ct, diag, best, ok in verdicts:
        mark = "OK  " if ok else "MISS"
        print(f"  {mark} {ct:17s} expected {diag:12s} best was {best:12s}"
              f" ({measures[ct]})")
    print()
    primary = wins >= 3
    print(f"PRIMARY  diagonal dominance: {wins}/4 rows  ->  "
          f"{'PASS' if primary else 'FAIL'} (>=3 required)")

    sec = sum(v["passes_p05"] for v in nulls.values())
    print(f"SECONDARY beat a size-matched null at p<0.05: {sec}/4 rows")

    # The named failure mode: if every set ranks similarly everywhere, the
    # result is non-discriminating and is a FAILURE regardless of the diagonal.
    spreads = [(M.loc[ct].max() - M.loc[ct].min()) / M.loc[ct].median()
               for ct in M.index]
    print(f"\nrow spread (max-min)/median: "
          + ", ".join(f"{ct}={s:.2f}" for ct, s in zip(M.index, spreads)))
    print("  -- a row spread near 0 means that row does not discriminate at all;")
    print("     the pre-registration counts a uniformly-high matrix as a FAILURE.")

    out = {
        "matrix": M.to_dict(),
        "measures": measures,
        "total_pairs": total_pairs,
        "families": {k: v for k, v in fams.items()},
        "verdicts": [{"cell_type": c, "expected": d, "best": b, "diagonal_wins": bool(o)}
                     for c, d, b, o in verdicts],
        "primary_wins": int(wins), "primary_pass": bool(primary),
        "null": nulls, "row_spread": dict(zip(M.index, map(float, spreads))),
        "params": vars(args),
    }
    with open(os.path.join(args.out, "tier3_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    pd.DataFrame(detail).to_csv(os.path.join(args.out, "tier3_matrix.csv"), index=False)
    M.to_csv(os.path.join(args.out, "tier3_median_ranks.csv"))

    # Cache the rank matrices. Rebuilding them costs ~25 min (the first arm pays
    # a long numba warm-up), and any follow-up diagnostic on a failed tier needs
    # exactly these. Saved AFTER the verdict above is computed and written, so
    # nothing downstream can influence the pre-registered result.
    cache = os.path.join(args.out, "rank_matrices")
    os.makedirs(cache, exist_ok=True)
    for ct in tables:
        M_ct, idx_ct = mats[ct]
        np.save(os.path.join(cache, f"{ct}.npy"), M_ct)
    with open(os.path.join(cache, "index.json"), "w") as fh:
        json.dump({"tf_names": tf_names, "families": fams, "measures": measures,
                   "total_pairs": total_pairs}, fh, indent=2)
    log(f"cached rank matrices to {cache}/")
    log(f"wrote {args.out}/tier3_results.json, tier3_matrix.csv, tier3_median_ranks.csv")


if __name__ == "__main__":
    mp.set_start_method("fork", force=True)
    main()
