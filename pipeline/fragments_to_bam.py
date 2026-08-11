#!/usr/bin/env python
# ============================================================================
# fragments_to_bam.py
# ----------------------------------------------------------------------------
# Build a pseudobulk BAM for one cell type (optionally one injury arm) directly
# from the 21 per-sample atac_fragments.tsv.gz files.
#
# WHY THIS EXISTS. TOBIAS ATACorrect requires --bam and has no fragment input;
# this dataset ships only fragments (no BAMs exist, and GEO does not deposit
# them). Fragments encode exactly what ATACorrect consumes -- Tn5 cut sites and
# strand -- and sequence context comes from --genome, not from read sequences.
# So this is a format adapter, not a fabrication. Published precedent for the
# same manoeuvre: IReNA-v2 (Pinlyu3/IReNA-v2, STEP4.2/4.3), which converts
# cell-type fragments to BAM and runs ATACorrect --read_shift 0 0.
#
# *** WHY NOT USE ML/celltype_fragments/*.bed.gz FROM THE SIBLING REPO ***
# Those were written by SplitFragments with BARE 10x barcodes, no sample
# prefix, so cells from different samples sharing a barcode collapse together
# (measured: Microglia 17,151 unique barcodes for 17,341 cells). Building from
# the per-sample files keyed on (sample, barcode) avoids that, and lets us drop
# PEAKQC-failing cells in the same pass. Barcodes are unique WITHIN a sample --
# verified, 0 duplicate (sample, barcode) pairs across all 67,072 cells.
#
# ---------------------------------------------------------------------------
# COORDINATES -- the part that is easy to get silently wrong
# ---------------------------------------------------------------------------
# CellRanger has ALREADY applied the Tn5 shift. From the 10x Cell Ranger ARC
# docs (our assay): "The start of the interval is moved forward by 4bp from a
# left-most alignment position and backward 5bp from the right-most alignment
# position." Meanwhile the TOBIAS FAQ says ATACorrect shifts +4/-5 internally.
# Both are true, so the default double-shifts and smears every footprint by
# 9 bp with no error message.  => ALWAYS run ATACorrect with --read_shift 0 0.
#
# TOBIAS computes the cutsite as (tobias/utils/ngs.pyx, get_cutsite), for reads
# with no soft-clipping:
#       forward: cutsite = reference_start + 1 + pos_shift
#       reverse: cutsite = reference_end   + 1 + neg_shift
# With --read_shift 0 0 this is reference_start+1 and reference_end+1 (1-based).
# We therefore place, for a fragment line (chrom, start, end) [BED, 0-based]:
#       forward read: reference_start = start   -> cutsite 0-based = start
#       reverse read: reference_end   = end     -> cutsite 0-based = end
# i.e. the two cut sites land exactly on the two fragment endpoints, symmetric,
# and separated by the fragment length.
#
# NOTE we deliberately differ from IReNA-v2 by 1 bp on the forward strand: their
# Read_fragment_to_GR loads the 0-based BED start straight into a 1-based
# IRanges start, so their forward cutsite sits 1 bp left of the insertion. That
# is immaterial for aggregate footprints and not immaterial for per-site calls.
#
# The remaining genuine ambiguity is whether the right-hand insertion is at BED
# `end` or `end-1`; the 10x wording does not settle it. TOBIAS estimates forward
# and reverse bias SEPARATELY (ReadList.split_strands), so a systematic 1 bp
# reverse offset largely absorbs into the learned reverse bias motif -- and a
# 1 bp asymmetry between the forward and reverse profiles in _atacorrect.pdf is
# exactly the signature if we chose wrong. CHECK THAT PLOT.
#
# ---------------------------------------------------------------------------
# READ CONSTRUCTION -- driven by what ATACorrect actually requires
# ---------------------------------------------------------------------------
# ReadList.from_bam filters on exactly two things:
#       read.is_unmapped == False and read.is_duplicate == False
# No MAPQ threshold, no proper-pair requirement. bias_estimation additionally
# skips reads whose leading CIGAR op is not a match, so we emit a plain <n>M.
#
# ONE READ PAIR PER FRAGMENT LINE. Column 5 of the fragments file is the number
# of read pairs supporting the fragment, i.e. it counts PCR duplicates that
# CellRanger already collapsed. Expanding by it would re-introduce them.
# (This is also why PEAKQC's n_fragments is ~2.3x the true fragment count --
# it sums column 5. See pipeline/count_unique_fragments.sh.)
#
# Reads are min(--read-len, fragment length) long and lie WHOLLY INSIDE the
# fragment, so they can never run off a contig end.
#
# Usage (from the tobias repo root):
#   python pipeline/fragments_to_bam.py --celltype Endothelial
#   python pipeline/fragments_to_bam.py --celltype Endothelial --injury U
#   python pipeline/fragments_to_bam.py --celltype Ependymal --no-filter
# Output:
#   bam/<celltype>[_<injury>].bam (+ .bai), coordinate-sorted
#   bam/<celltype>[_<injury>].stats.txt
# ============================================================================
import argparse
import os
import subprocess
import sys
import time

import pandas as pd
import pysam

CHROMBPNET = "../TF_cobinding_multiome_chrombpnet"
METADATA = "meta/cell_metadata.csv"
SCORES = "meta/peakqc_scores.csv"
CHROM_SIZES = os.path.join(CHROMBPNET, "meta", "mm10.chr.sizes")
BAM_DIR = "bam"
FLD_THRESHOLD = 100.0     # PEAKQC rule of thumb; see qc/peakqc_summary.txt


def log(*a):
    print(time.strftime("%H:%M:%S"), "|", *a, flush=True)


def fragments_path(sample):
    return os.path.join(CHROMBPNET, "data", sample, "outs", "atac_fragments.tsv.gz")


def build_keepset(celltype, injury, apply_filter):
    """(sample, barcode) pairs to include -> dict of sample -> set(barcode)."""
    md = pd.read_csv(METADATA, usecols=["cell", "barcode", "sample", "celltype",
                                        "cluster_ids_injury"])
    md = md[md["celltype"] == celltype]
    if injury:
        md = md[md["cluster_ids_injury"] == f"{celltype}_{injury}"]
    n_before = len(md)

    if apply_filter:
        sc = pd.read_csv(SCORES, usecols=["cell", "fld_score"])
        md = md.merge(sc, on="cell", how="left")
        missing = md["fld_score"].isna().sum()
        if missing:
            sys.exit(f"{missing} cells have no PEAKQC score -- refusing to guess")
        md = md[md["fld_score"] >= FLD_THRESHOLD]

    keep = {s: set(g["barcode"]) for s, g in md.groupby("sample")}
    log(f"{celltype}{'_' + injury if injury else ''}: "
        f"{len(md):,} of {n_before:,} cells kept across {len(keep)} samples"
        f"{'' if apply_filter else '  (NO PEAKQC FILTER)'}")
    return keep


def write_reads(bam, keep, read_len, tid_of):
    """Stream all samples, emit one read pair per kept fragment line."""
    n_frag = n_skip_chrom = 0
    for sample in sorted(keep):
        frag = fragments_path(sample)
        if not os.path.exists(frag):
            sys.exit(f"missing {frag}")
        barcodes = keep[sample]
        t0, n_this = time.time(), 0
        with pysam.BGZFile(frag, "rb") as fh:
            for raw in fh:
                line = raw.decode()
                if line[0] == "#":
                    continue
                chrom, start, end, bc, _count = line.rstrip("\n").split("\t")
                if bc not in barcodes:
                    continue
                tid = tid_of.get(chrom)
                if tid is None:          # non-standard contig, not in mm10.chr.sizes
                    n_skip_chrom += 1
                    continue
                start, end = int(start), int(end)
                L = min(read_len, end - start)
                if L < 1:
                    continue

                name = f"{sample}:{bc}:{n_frag}"
                tlen = end - start
                # forward mate: reference_start == fragment start
                a = pysam.AlignedSegment()
                a.query_name = name
                a.flag = 99            # paired, proper, mate reverse, first in pair
                a.reference_id = tid
                a.reference_start = start
                a.mapping_quality = 60
                a.cigartuples = [(0, L)]
                a.next_reference_id = tid
                a.next_reference_start = end - L
                a.template_length = tlen
                a.query_sequence = "N" * L
                a.query_qualities = pysam.qualitystring_to_array("I" * L)
                a.set_tag("CB", bc)
                bam.write(a)
                # reverse mate: reference_end == fragment end
                b = pysam.AlignedSegment()
                b.query_name = name
                b.flag = 147           # paired, proper, reverse, second in pair
                b.reference_id = tid
                b.reference_start = end - L
                b.mapping_quality = 60
                b.cigartuples = [(0, L)]
                b.next_reference_id = tid
                b.next_reference_start = start
                b.template_length = -tlen
                b.query_sequence = "N" * L
                b.query_qualities = pysam.qualitystring_to_array("I" * L)
                b.set_tag("CB", bc)
                bam.write(b)
                n_frag += 1
                n_this += 1
        log(f"  {sample}: {n_this:,} fragments in {time.time()-t0:.0f}s")
    return n_frag, n_skip_chrom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--injury", choices=["I", "U"], default=None,
                    help="restrict to one injury arm (default: both pooled)")
    ap.add_argument("--no-filter", action="store_true",
                    help="skip the PEAKQC fld_score >= 100 filter")
    ap.add_argument("--read-len", type=int, default=36)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    sizes = pd.read_csv(CHROM_SIZES, sep="\t", header=None, names=["chrom", "size"])
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": int(n)} for c, n in
                     zip(sizes["chrom"], sizes["size"])]}
    tid_of = {c: i for i, c in enumerate(sizes["chrom"])}

    keep = build_keepset(args.celltype, args.injury, not args.no_filter)
    if not keep:
        sys.exit("no cells selected")

    os.makedirs(BAM_DIR, exist_ok=True)
    tag = args.celltype + (f"_{args.injury}" if args.injury else "")
    final = os.path.join(BAM_DIR, f"{tag}.bam")
    unsorted_bam = os.path.join(BAM_DIR, f"{tag}.unsorted.bam")

    t0 = time.time()
    # Written unsorted: the 21 source files are each coordinate-sorted, but
    # interleaving them is not, so samtools sort does the global ordering.
    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as bam:
        n_frag, n_skip = write_reads(bam, keep, args.read_len, tid_of)
    log(f"wrote {n_frag:,} fragments ({2*n_frag:,} reads) in {time.time()-t0:.0f}s; "
        f"{n_skip:,} skipped on non-standard contigs")

    log("sorting")
    subprocess.run(["samtools", "sort", "-@", str(args.threads),
                    "-o", final, unsorted_bam], check=True)
    subprocess.run(["samtools", "index", "-@", str(args.threads), final], check=True)
    os.remove(unsorted_bam)

    stats = subprocess.run(["samtools", "flagstat", final],
                           capture_output=True, text=True, check=True).stdout
    with open(os.path.join(BAM_DIR, f"{tag}.stats.txt"), "w") as fh:
        fh.write(f"cells: {sum(len(v) for v in keep.values()):,}\n")
        fh.write(f"fragments: {n_frag:,}\nreads: {2*n_frag:,}\n")
        fh.write(f"skipped_nonstandard_contig: {n_skip:,}\n\n{stats}")
    log(f"DONE -> {final}")
    print(stats)


if __name__ == "__main__":
    main()
