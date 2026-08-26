#!/usr/bin/env python
# ============================================================================
# plot_differential_heatmap.py
# ----------------------------------------------------------------------------
# Analogue of Zamboni et al. Fig. 4a: which TF motifs are most differentially
# active after injury, across cell types.
#
# THEIRS vs OURS. Fig. 4a is motif ENRICHMENT in injury-responsive regions,
# adjusted for TF expression. This is footprint-based differential BINDING over
# each cell type's full peak set. Different quantities -- only rank agreement
# among the top factors is meaningful (VALIDATION_CRITERIA.md divergence 1).
#
# NEURONS ARE PLOTTED, AND SEPARATED BY A GAP. They are the standing negative
# control for any "glia-specific" claim (statistical policy item 5). On the
# 2026-08-20 results AP-1 is just as strong in neurons, and this figure is
# meant to show that honestly rather than crop it out.
#
# Colour is DIVERGING (blue = down after injury, red = up) with a neutral grey
# midpoint, because `change` is signed and centred on zero.
#
# Usage from the repo root:
#   python pipeline/plot_differential_heatmap.py [--top 6] [--outdir figures]
# ============================================================================
import argparse, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

GLIA    = ["Astrocytes", "Oligodendrocytes", "OPCs", "Microglia"]
NEURONS = ["Neurons_V", "Neurons_D"]
AP1     = re.compile(r"FOS|JUN|ATF3|BATF|JDP2", re.I)   # JDP2 = Jun dimerization protein 2, binds the AP-1/TRE site

# diverging ramp: blue arm / neutral / red arm, lightness-matched per arm
BLUE = ["#0d366b", "#184f95", "#1c5cab", "#2a78d6", "#6da7ec", "#9ec5f4", "#cde2fb"]
GREY = "#f0efec"
RED  = ["#fbd5d4", "#f6b1b0", "#ef8785", "#e34948", "#c2302f", "#9a2322", "#74191a"]
CMAP = LinearSegmentedColormap.from_list("bl_gy_rd", BLUE + [GREY] + RED, N=256)

INK, INK2, MUTED, RULE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"


def load(ct, root):
    d = pd.read_csv(f"{root}/{ct}/{ct}_results.txt", sep="\t")
    col = f"{ct}_I_{ct}_U_change"
    return d.set_index("name")[col].rename(ct)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=12, help="top TFs per cell type to union")
    ap.add_argument("--per-family", type=int, default=2,
                    help="max members shown per redundant family (AP-1, C/EBP, NFI, RFX)")
    ap.add_argument("--root", default="tobias/bindetect")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--tag", default="full_peak_set")
    args = ap.parse_args()

    cts = GLIA + NEURONS
    M = pd.concat([load(ct, args.root) for ct in cts], axis=1)

    # rows = union of each cell type's top-N, plus every factor the paper names
    paper = ["STAT3", "CEBPD", "TFEB", "TFEC", "BHLHE40", "HIF1A", "CTCF"]
    keep = set()
    for ct in cts:
        keep |= set(M[ct].sort_values(ascending=False).head(args.top).index)
    keep |= {p for p in paper if p in M.index}
    M = M.loc[sorted(keep, key=lambda t: -M.loc[t, GLIA].mean())]

    # COLLAPSE REDUNDANT FAMILY MEMBERS. AP-1 alone matches 33 of the 746 motifs
    # (mostly dimer matrices) and its members score near-identically, so an
    # un-collapsed top-N is simply the AP-1 block repeated. Keep the strongest
    # `--per-family` members of each redundant family so the rows show
    # biological diversity rather than motif-database structure. This is a
    # DISPLAY choice only -- it changes no score. See
    # [[motif-family-size-bias]] for why family size distorts naive rankings.
    def family(tf):
        if AP1.search(tf):                 return "AP-1"
        if re.match(r"CEBP", tf, re.I):    return "C/EBP"
        if re.match(r"NFI[ABCX]$", tf):    return "NFI"
        if re.match(r"RFX", tf, re.I):     return "RFX"
        return None
    seen, rows = {}, []
    for tf in M.index:
        f = family(tf)
        if f is None:
            rows.append(tf); continue
        seen[f] = seen.get(f, 0) + 1
        if seen[f] <= args.per_family:
            rows.append(tf)
    M = M.loc[rows]

    lim = float(np.nanmax(np.abs(M.values)))
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)

    nrow, ncol = M.shape
    fig, (axg, axn) = plt.subplots(
        1, 2, figsize=(7.4, 0.30 * nrow + 2.0),
        gridspec_kw={"width_ratios": [len(GLIA), len(NEURONS)], "wspace": 0.10})

    for ax, cols, title in ((axg, GLIA, "Glia"), (axn, NEURONS, "Neurons  (negative control)")):
        sub = M[cols]
        im = ax.imshow(sub.values, cmap=CMAP, norm=norm, aspect="auto")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=38,
                           ha="right", fontsize=8.5, color=INK2)
        ax.set_title(title, fontsize=9, color=INK2, pad=8, loc="left")
        ax.set_yticks(range(nrow))
        # selective direct labels: only cells with a substantial effect
        for i in range(nrow):
            for j in range(len(cols)):
                v = sub.values[i, j]
                if abs(v) >= 0.45:
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.4,
                            color="#ffffff" if abs(v) > lim * 0.55 else INK)
        ax.set_xticks(np.arange(-.5, len(cols), 1), minor=True)
        ax.set_yticks(np.arange(-.5, nrow, 1), minor=True)
        ax.grid(which="minor", color="#fcfcfb", linewidth=1.4)
        ax.tick_params(which="minor", length=0)
        for s in ax.spines.values():
            s.set_visible(False)

    axg.set_yticklabels(M.index, fontsize=7.4,
                        color=INK)
    for lbl, tf in zip(axg.get_yticklabels(), M.index):
        if AP1.search(tf):
            lbl.set_color("#c2302f"); lbl.set_fontweight("bold")
    axn.set_yticklabels([])
    axn.tick_params(axis="y", length=0)

    cb = fig.colorbar(im, ax=(axg, axn), fraction=0.030, pad=0.03)
    cb.set_label("differential binding score\n(+ = more bound after injury)",
                 fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7.5, color=MUTED, labelcolor=INK2)
    cb.outline.set_visible(False)

    # place the title block in absolute inches so it never collides as nrow varies
    H = fig.get_figheight()
    fig.subplots_adjust(top=1 - 0.95 / H)
    fig.suptitle("Differential TF binding after spinal cord injury",
                 fontsize=12.5, color=INK, x=0.055, ha="left", y=1 - 0.22 / H)
    fig.text(0.055, 1 - 0.55 / H,
             "AP-1 family in red.  Analogue of Zamboni et al. Fig. 4a — see script header "
             "for why the quantities differ.",
             fontsize=7.8, color=MUTED, ha="left", va="top")

    os.makedirs(args.outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        p = f"{args.outdir}/differential_heatmap_{args.tag}.{ext}"
        fig.savefig(p, bbox_inches="tight", dpi=200, facecolor="#fcfcfb")
        print("wrote", p)


if __name__ == "__main__":
    main()
