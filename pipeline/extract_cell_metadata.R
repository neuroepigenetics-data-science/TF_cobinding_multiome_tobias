# ============================================================================
# extract_cell_metadata.R
# ----------------------------------------------------------------------------
# One-time export of the merged object's cell metadata to a small CSV, so that
# nothing downstream has to load the 7.2 GB object again (it needs ~18.8 GB
# resident on a 24 GB machine).
#
# WHY THIS EXISTS. PEAKQC scores cells individually, so every fragment must be
# attributable to exactly one cell. The pooled per-cell-type fragment BEDs in
# the chrombpnet repo CANNOT do that: SplitFragments wrote plain 10x barcodes
# with no sample prefix, so cells from different samples that happen to share a
# barcode are merged. Measured: Microglia has 17,151 unique barcodes for 17,341
# cells (190 collisions), Astrocytes 6,604/6,625, Ependymal 5,613/5,625.
# Barcodes ARE unique within a sample, so PEAKQC must be run per sample against
# the original data/<sample>/outs/atac_fragments.tsv.gz. That requires knowing
# which barcodes belong to which sample -- which is what this exports.
#
# Writes: meta/cell_metadata.csv  — one row per cell:
#           cell        the object's cell name (unique)
#           barcode     the bare 10x barcode as it appears in the fragments file
#           sample      sample id, matching pipeline/samples.csv
#           celltype    cluster_ids
#           condition   U / 1dpi / 3dpi / 7dpi / 28dpi
#           (plus any other metadata columns present, verbatim)
#
# Run from the tobias repo root:
#   Rscript pipeline/extract_cell_metadata.R [path/to/object.rds]
# ============================================================================

suppressPackageStartupMessages({ library(Seurat) })

args <- commandArgs(trailingOnly = TRUE)
obj_path <- if (length(args) >= 1) args[[1]] else
  "../TF_cobinding_multiome_chrombpnet/multiome_obj_clean_240326.rds"

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)

if (!file.exists(obj_path)) stop("object not found: ", obj_path)
say("reading ", obj_path, " (this takes a few minutes and ~19 GB RAM)")
obj <- readRDS(obj_path)
say("loaded: ", ncol(obj), " cells")

md <- obj@meta.data
md$cell <- rownames(md)

say("metadata columns: ", paste(colnames(md), collapse = ", "))

# --- work out how the cell name encodes the barcode -------------------------
# Reported so the mapping is visible rather than assumed. Seurat merges usually
# produce "<something>_<BARCODE>" or "<BARCODE>_<n>"; we need the bare barcode
# exactly as it appears in column 4 of atac_fragments.tsv.gz (e.g. AAACGAAAG...-1).
say("example cell names: ", paste(head(md$cell, 3), collapse = " | "))

bc <- sub(".*_([ACGT]+-[0-9]+)$", "\\1", md$cell)   # trailing barcode form
hit <- grepl("^[ACGT]+-[0-9]+$", bc)
if (!all(hit)) {
  bc2 <- sub("^([ACGT]+-[0-9]+)_.*$", "\\1", md$cell)  # leading barcode form
  hit2 <- grepl("^[ACGT]+-[0-9]+$", bc2)
  if (sum(hit2) > sum(hit)) { bc <- bc2; hit <- hit2 }
}
say("barcode parse: ", sum(hit), " / ", length(hit), " cell names matched")
if (!all(hit)) {
  say("UNPARSED examples: ", paste(head(md$cell[!hit], 5), collapse = " | "))
  warning("some cell names did not yield a barcode; inspect before using")
}
md$barcode <- bc

# --- sample column ----------------------------------------------------------
# Do not guess: report the candidates and pick the one that looks like a sample
# id, else leave it to be derived from the cell-name prefix.
cand <- intersect(c("sample", "sample_id", "orig.ident", "dataset", "library"),
                  colnames(md))
say("sample-like columns present: ", if (length(cand)) paste(cand, collapse = ", ") else "NONE")
if (length(cand)) {
  md$sample <- as.character(md[[cand[1]]])
  say("using '", cand[1], "' as sample; ", length(unique(md$sample)), " distinct values")
} else {
  md$sample <- sub("_([ACGT]+-[0-9]+)$", "", md$cell)
  say("derived sample from cell-name prefix; ", length(unique(md$sample)), " distinct values")
}
print(table(md$sample))

# --- sanity: barcode must be unique WITHIN each sample ----------------------
dup <- sum(duplicated(paste(md$sample, md$barcode)))
say("duplicate (sample, barcode) pairs: ", dup, if (dup == 0) "  <- good" else "  <- PROBLEM")

dir.create("meta", showWarnings = FALSE)
front <- c("cell", "barcode", "sample", "celltype", "condition")
if ("cluster_ids" %in% colnames(md)) md$celltype <- as.character(md$cluster_ids)
front <- front[front %in% colnames(md)]
md <- md[, c(front, setdiff(colnames(md), front))]

out <- "meta/cell_metadata.csv"
write.csv(md, out, row.names = FALSE)
say("wrote ", out, ": ", nrow(md), " cells x ", ncol(md), " columns")
say("cell types:"); print(table(md$celltype))
