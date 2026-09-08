#!/usr/bin/env python
# ============================================================================
# check_atacbias.py
# ----------------------------------------------------------------------------
# Validate an ATACorrect run by asking where the LEARNED Tn5 bias motif is
# centred. Run this on every *_AtacBias.pickle before trusting any footprint.
# It takes ~1 second per file and it is the only thing standing between us and
# a silently wrong result.
#
# WHY IT WORKS. Tn5 binds as a dimer and its sequence-preference motif is
# palindromic about the insertion point. So the motif itself tells us which
# coordinate ATACorrect actually treated as the cut site:
#     best offset  0   -> we placed the cut sites correctly
#     best offset +-9  -> the 9 bp double shift (--read_shift left at 4 -5)
#     best offset +-1  -> the BED `end` off-by-one
# The forward and reverse matrices are learned from DISJOINT read sets, so
# their agreement is independent evidence, not an artefact of the same data.
#
# This is a script-ified version of the ad-hoc check run on 2026-08-11 that
# validated the BAM synthesis; the method and the reference numbers are in
# 2026-08-11.txt section 4. Reference values for Endothelial:
#     symmetry at offset 0 = 0.9775   (next best -5 = 0.5243)
#     forward/reverse correlation = 0.999926 at lag 0
#
# Usage (from the tobias repo root):
#   python pipeline/check_atacbias.py tobias/atacorrect/*/*_AtacBias.pickle
#   python pipeline/check_atacbias.py --quiet tobias/atacorrect/Microglia_I/*.pickle
# Exit status is nonzero if ANY file fails, so it can gate a pipeline.
# ============================================================================
import argparse
import glob
import os
import pickle
import sys

import numpy as np

# Endothelial scored 0.9775 with 1.8M fragments. Deeper cell types should do at
# least as well; these are floors, not targets.
MIN_SYMMETRY = 0.90
MIN_STRAND_CORR = 0.99


def profile(matrix):
    """Per-position bias signal: sum of |pssm| over bases. 1-D, length 25."""
    return np.abs(np.array(matrix.pssm)).sum(axis=0)


def strand_agreement(f, r):
    """Correlation of the two strand profiles, and the lag that aligns them."""
    def corr_at(s):
        a = f[max(0, s):len(f) + min(0, s)]
        b = r[max(0, -s):len(r) + min(0, -s)]
        return np.corrcoef(a, b)[0, 1]
    best_lag = max(range(-9, 10), key=corr_at)
    return np.corrcoef(f, r)[0, 1], np.max(np.abs(f - r)), best_lag


def symmetry_scores(f):
    """For each candidate centre, how mirror-symmetric is the profile about it."""
    k = len(f) // 2
    out = {}
    for off in range(-9, 10):
        c = k + off
        n = min(c, len(f) - 1 - c)
        if n < 6:                       # too little overlap to mean anything
            continue
        left = f[c - n:c][::-1]
        right = f[c + 1:c + 1 + n]
        out[off] = np.corrcoef(left, right)[0, 1]
    return out


def check(path, quiet=False):
    bias = pickle.load(open(path, "rb")).bias
    f, r = profile(bias["forward"]), profile(bias["reverse"])

    corr, maxdiff, lag = strand_agreement(f, r)
    scores = symmetry_scores(f)
    best_off = max(scores, key=scores.get)
    sym0 = scores.get(0, float("nan"))

    problems = []
    if best_off != 0:
        problems.append(f"motif is centred at offset {best_off:+d}, not 0")
    if lag != 0:
        problems.append(f"strands disagree by {lag:+d} bp")
    if sym0 < MIN_SYMMETRY:
        problems.append(f"symmetry at 0 is {sym0:.4f} < {MIN_SYMMETRY}")
    if corr < MIN_STRAND_CORR:
        problems.append(f"strand correlation {corr:.6f} < {MIN_STRAND_CORR}")

    name = os.path.basename(path).replace("_AtacBias.pickle", "")
    verdict = "PASS" if not problems else "FAIL"
    print(f"{verdict}  {name:22s} symmetry@0={sym0:.4f}  best_offset={best_off:+d}  "
          f"strand_corr={corr:.6f}  lag={lag:+d}")

    if not quiet or problems:
        top = sorted(scores, key=lambda o: -scores[o])[:4]
        print("        symmetry by offset: "
              + "  ".join(f"{o:+d}:{scores[o]:.4f}" for o in top))
        print(f"        strand max abs diff: {maxdiff:.4f}")
    for p in problems:
        print(f"        !! {p}")
    if problems:
        print("        >> DO NOT TRUST THIS RUN. Offset +-9 means the 9 bp double")
        print("           shift: check that ATACorrect got --read_shift 0 0.")
    return not problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pickles", nargs="+", help="*_AtacBias.pickle paths (globs ok)")
    ap.add_argument("--quiet", action="store_true",
                    help="one line per file unless it fails")
    args = ap.parse_args()

    paths = sorted({p for arg in args.pickles for p in glob.glob(arg)})
    if not paths:
        sys.exit("no pickles matched")

    ok = [check(p, args.quiet) for p in paths]
    print(f"\n{sum(ok)}/{len(ok)} passed")
    sys.exit(0 if all(ok) else 1)


if __name__ == "__main__":
    main()
