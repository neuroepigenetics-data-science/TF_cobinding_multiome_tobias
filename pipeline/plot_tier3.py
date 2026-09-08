#!/usr/bin/env python
# ============================================================================
# plot_tier3.py -- the two figures that carry the Tier 3 argument
# ----------------------------------------------------------------------------
# Two figures, not a gallery. Each one exists to make a specific claim legible
# that a table of numbers does not.
#
# FIGURE 1 (the important one) -- "why Tier 3 failed, and why microglia didn't".
# For every non-AP-1 motif: how similar its motif is to AP-1 (x) against how
# strongly it co-occurs with AP-1 (y). The diagnosis is a SHAPE: in Astrocytes
# and Oligodendrocytes the strongest "partners" pile up on the LEFT -- motifs
# that resemble AP-1, i.e. the same DNA matched twice. In Microglia the top of
# the plot sits on the RIGHT, genuinely distinct factors. One figure therefore
# shows both the confound and the single clean result, which is why it is worth
# more than a plot of the failed matrix.
#
# FIGURE 2 -- the verdict itself. Median rank of each partner set per cell type,
# with the diagonal (the paper's expectation) marked. Encoded as TWO classes
# (diagonal vs the rest) rather than four coloured series, because the only
# question asked is "does the diagonal win its row?" -- and four categorical
# hues would exceed the validated all-pairs cap for this chart form.
#
# Colour: categorical slots 1 and 2 of the validated palette (blue #2a78d6,
# orange #eb6834); bulk points are neutral ink, which is context rather than a
# series. Surface #fcfcfb and the ink ramp match the existing project figures
# (pipeline/plot_differential_heatmap.py), so the set reads as one system.
# Identity is never colour-alone: every highlighted motif is direct-labelled and
# both figures carry a legend.
#
# Usage: ~/venvs/tfcomb/bin/python pipeline/plot_tier3.py
# ============================================================================
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SURFACE = "#fcfcfb"
INK, INK2, MUTED, RULE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE, ORANGE = "#2a78d6", "#eb6834"          # validated categorical slots 1, 2
DOT = "#c9c7c0"                               # bulk points: context, not a series

CELLS = ["Astrocytes", "Oligodendrocytes", "Microglia", "Ependymal_I"]
LINEAGE = {"SOX9": "SOX9", "SOX10": "SOX10", "NFIA": "NFI", "NFIB": "NFI",
           "NFIC": "NFI", "NFIX": "NFI", "RFX1": "RFX", "RFX2": "RFX",
           "RFX3": "RFX", "RFX4": "RFX", "ELK3": "ELK3", "CEBPA": "CEBPA",
           "TCF4": "TCF4"}
AP1_LIKE_CUT = 0.70      # for the shaded band only; never used to score anything


def load(dirpath, dist_root):
    meta = json.load(open(os.path.join(dirpath, "rank_matrices", "index.json")))
    tf_names, fams = meta["tf_names"], meta["families"]
    idx = {t: i for i, t in enumerate(tf_names)}
    ap1 = [t for t in fams["AP-1"] if t in idx]
    out = {}
    for ct in CELLS:
        M = np.load(os.path.join(dirpath, "rank_matrices", f"{ct}.npy"))
        sub = M[[idx[t] for t in ap1], :]
        with np.errstate(invalid="ignore"):
            best = np.nanmin(sub, axis=0)
        dp = os.path.join(dist_root, ct, f"{ct}_distances.txt")
        dn = [n for n in pd.read_csv(dp, sep="\t", nrows=0).columns if n != "#"]
        DM = pd.DataFrame(np.loadtxt(dp, skiprows=1), index=dn, columns=dn)
        core = [c for c in dn if c.split("_")[0] in ("FOSJUN", "FOSL1", "FOS", "JUNB", "JUN")]
        rows = []
        for j, t in enumerate(tf_names):
            if t in set(ap1) or t not in DM.index or np.isnan(best[j]):
                continue
            rows.append({"tf": t, "short": t.split("_")[0],
                         "dist": float(DM.loc[t, core].astype(float).min()),
                         "rank": float(best[j])})
        out[ct] = pd.DataFrame(rows)
    return out, meta


def place_labels(ax, items, min_gap_pt=10.5, half_width_pt=52.0):
    """Annotate with real collision avoidance, in dpi-independent units.

    Two earlier attempts failed in instructive ways. Absolute figure pixels
    break because savefig re-renders at dpi=200 with a tight bbox. A fixed
    ladder keyed to rank order breaks because two labels can be adjacent in
    rank yet land on top of each other, and two on opposite sides can collide
    in the middle of the panel.

    So: convert each point to AXES FRACTION (a ratio, therefore dpi-invariant),
    scale by the axes size in points, greedily push each label down until it
    clears every label already placed, and emit the result as an offset in
    points. Both sides are packed together, so left- and right-anchored labels
    cannot overlap either.
    """
    fig = ax.figure
    bb = ax.get_position()
    h_pt = bb.height * fig.get_figheight() * 72.0
    w_pt = bb.width * fig.get_figwidth() * 72.0
    inv = ax.transAxes.inverted()

    placed = []          # (x_centre_pt, y_pt)
    for x, y, text, colour, side in items:
        fx, fy = inv.transform(ax.transData.transform((x, y)))
        px, py = fx * w_pt, fy * h_pt
        # label box centre sits to one side of the point
        cx = px + (half_width_pt * 0.5 + 7) * (1 if side == "left" else -1)
        ty = py
        for _ in range(80):
            if not any(abs(ty - qy) < min_gap_pt and abs(cx - qx) < half_width_pt
                       for qx, qy in placed):
                break
            ty -= min_gap_pt        # push downward on screen (y grows upward)
        placed.append((cx, ty))
        ax.annotate(text, (x, y), fontsize=7.6, color=colour, fontweight="bold",
                    zorder=5, xytext=(7 if side == "left" else -7, ty - py),
                    textcoords="offset points",
                    ha="left" if side == "left" else "right", va="center")


def figure1(data, meta, outdir):
    """Colour the TOP hits by whether their motif resembles AP-1.

    An earlier draft highlighted the expected lineage partners instead. That
    version could not show its own claim: the factors that make Microglia
    interesting (TEF, DBP, NFIL3, CEBPG) were unlabelled grey dots, because they
    are neither AP-1-like nor on the paper's expected list. Labelling the top
    hits themselves is what makes the contrast between the panels visible.
    """
    TOPN = 8
    label_jobs = []
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0), facecolor=SURFACE)
    for ax, ct in zip(axes.ravel(), CELLS):
        d = data[ct].sort_values("rank")
        ax.set_facecolor(SURFACE)
        ax.axvspan(0, AP1_LIKE_CUT, color=RULE, alpha=0.5, lw=0, zorder=0)
        ax.scatter(d["dist"], d["rank"], s=9, c=DOT, lw=0, zorder=1)

        # the paper's expected partners: open rings, no labels (Figure 2 scores them)
        lin = d[d["short"].isin(LINEAGE)]
        ax.scatter(lin["dist"], lin["rank"], s=46, facecolors="none",
                   edgecolors=INK2, lw=1.1, zorder=2)

        top = d.head(TOPN)
        n_like = int((top["dist"] < AP1_LIKE_CUT).sum())
        items = []
        for _, r in top.iterrows():
            like = r["dist"] < AP1_LIKE_CUT
            c = ORANGE if like else BLUE
            ax.scatter(r["dist"], r["rank"], s=44, c=c, lw=0.9,
                       edgecolors=SURFACE, zorder=4)
            items.append((r["dist"], r["rank"], r["short"], c,
                          "left" if like else "right"))
        label_jobs.append((ax, items))

        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xlim(-0.05, 1.16)
        arm = "occupancy, single arm" if ct == "Ependymal_I" else "injured vs uninjured"
        ax.set_title(f"{ct}    ", fontsize=11.5, color=INK, loc="left", pad=20,
                     fontweight="bold")
        ax.text(0.0, 1.045, f"{arm}   —   {n_like} of the top {TOPN} are AP-1 look-alikes",
                transform=ax.transAxes, fontsize=8.2,
                color=ORANGE if n_like >= TOPN / 2 else BLUE)
        ax.tick_params(labelsize=8, colors=INK2, length=3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(RULE)
        ax.grid(True, color=RULE, lw=0.6, alpha=0.7, zorder=0)
        ax.set_axisbelow(True)

    for ax in axes[1]:
        ax.set_xlabel("motif distance from AP-1   (0 = identical motif,  1 = unrelated)",
                      fontsize=8.8, color=INK2)
    for ax in axes[:, 0]:
        ax.set_ylabel("strongest co-occurrence rank with AP-1\n(log scale, top = strongest)",
                      fontsize=8.8, color=INK2)

    handles = [
        Line2D([], [], marker="o", ls="", ms=6.5, mfc=ORANGE, mec=ORANGE,
               label="top hit whose motif RESEMBLES AP-1  (< 0.7) — an artifact"),
        Line2D([], [], marker="o", ls="", ms=6.5, mfc=BLUE, mec=BLUE,
               label="top hit that is a GENUINELY DISTINCT factor  (≥ 0.7)"),
        Line2D([], [], marker="o", ls="", ms=6.5, mfc="none", mec=INK2, mew=1.1,
               label="the paper's expected lineage partner (scored in Figure 2)"),
        Line2D([], [], marker="o", ls="", ms=5, mfc=DOT, mec=DOT,
               label="all other motifs"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=8.5, labelcolor=INK2, bbox_to_anchor=(0.5, -0.012))
    fig.suptitle("Why Tier 3 failed — the strongest AP-1 “partners” are mostly motifs "
                 "that resemble AP-1 itself",
                 fontsize=13.5, color=INK, x=0.010, ha="left", y=0.998,
                 fontweight="bold")
    fig.text(0.010, 0.960,
             "Astrocytes and Oligodendrocytes fill their strongest hits with AP-1 look-alikes on "
             "the LEFT — the same stretch of DNA matched by two similar motifs, not two factors\n"
             "binding together. Microglia is the exception: its top hits are distinct factors on "
             "the RIGHT (TEF, DBP, CEBPG — a PAR-bZIP/C-EBP signature), and it is the one cell "
             "type that recovers the paper's biology.",
             fontsize=8.9, color=INK2, ha="left", va="top")
    fig.tight_layout(rect=[0, 0.050, 1, 0.930])
    for ax, items in label_jobs:
        place_labels(ax, items)
    save(fig, outdir, "tier3_motif_similarity_confound")


def figure2(dirpath, outdir):
    M = pd.read_csv(os.path.join(dirpath, "tier3_median_ranks.csv"), index_col=0)
    res = json.load(open(os.path.join(dirpath, "tier3_results.json")))
    diag = {v["cell_type"]: v["expected"] for v in res["verdicts"]}
    order = [c for c in CELLS if c in M.index]

    fig, ax = plt.subplots(figsize=(9.9, 4.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for y, ct in enumerate(order):
        row = M.loc[ct]
        ax.plot([row.min(), row.max()], [y, y], color=RULE, lw=1.6, zorder=1,
                solid_capstyle="round")
        for col in M.columns:
            is_d = col == diag[ct]
            ax.scatter(row[col], y, s=104 if is_d else 62,
                       c=BLUE if is_d else DOT, zorder=3,
                       edgecolors=SURFACE, lw=1.6)
        best = row.idxmin()
        ax.annotate(f"{diag[ct]}", (row[diag[ct]], y), fontsize=8,
                    color=BLUE, fontweight="bold", xytext=(0, 11),
                    textcoords="offset points", ha="center")
        if best != diag[ct]:
            ax.scatter(row[best], y, s=104, facecolors="none", edgecolors=ORANGE,
                       lw=1.8, zorder=4)
            ax.annotate(f"{best} wins", (row[best], y), fontsize=8, color=ORANGE,
                        xytext=(0, -17), textcoords="offset points", ha="center")

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=9.6, color=INK)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("median rank of AP-1 ↔ partner-set pairs   "
                  "(out of ~277,800 pairs; further LEFT = stronger)",
                  fontsize=8.8, color=INK2)
    ax.tick_params(labelsize=8.4, colors=INK2, length=3)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.grid(True, axis="x", color=RULE, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.margins(x=0.10, y=0.22)

    handles = [
        Line2D([], [], marker="o", ls="", ms=8, mfc=BLUE, mec=BLUE,
               label="the paper's expectation for that cell type (the diagonal)"),
        Line2D([], [], marker="o", ls="", ms=8, mfc="none", mec=ORANGE, mew=1.8,
               label="the set that actually ranked best"),
        Line2D([], [], marker="o", ls="", ms=7, mfc=DOT, mec=DOT,
               label="the other partner sets"),
    ]
    # below the axes: at "lower right" the legend sat on top of the Ependymal row
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=8.2, labelcolor=INK2, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Tier 3 as scored — the expected partner wins only 1 of 4 rows "
                 "(3 were required)",
                 fontsize=12.6, color=INK, x=0.012, ha="left", fontweight="bold")
    fig.text(0.012, 0.895,
             "Ependymal_I is a single-arm occupancy quantity and is compared only within its own row.",
             fontsize=8.4, color=MUTED, ha="left", va="top")
    fig.tight_layout(rect=[0, 0.055, 1, 0.895])
    save(fig, outdir, "tier3_diagonal_dominance")


def save(fig, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(p, bbox_inches="tight", dpi=200, facecolor=SURFACE)
        print("wrote", p)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="tobias/tfcomb")
    ap.add_argument("--dist-root", default="tobias/bindetect")
    ap.add_argument("--out", default="figures")
    args = ap.parse_args()
    data, meta = load(args.dir, args.dist_root)
    figure1(data, meta, args.out)
    figure2(args.dir, args.out)


if __name__ == "__main__":
    main()
