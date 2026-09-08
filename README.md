# TF_cobinding_multiome_tobias

A **GPU-free route to transcription-factor occupancy and co-binding** on the Zamboni
et al. spinal-cord-injury multiome dataset, built with TOBIAS and TF-COMB.

    PEAKQC  ->  pseudobulk BAM  ->  ATACorrect  ->  ScoreBigwig  ->  BINDetect  ->  TF-COMB
    cell QC     per cell type      Tn5 debias     footprint       bound sites    co-binding

---

## Why this exists

**The purpose is pipeline validation, not discovery.** The plan is to reproduce a
*known* answer — the published result on the Zamboni dataset — so that the pipeline
can be trusted later on **the lab's own data, where no published comparison will
exist**. That is the whole logic of the project and it explains most of the design
choices, including the motif set (JASPAR2020 CORE vertebrates, matching the paper's
supervised analysis) and the pre-registered validation criteria.

It is GPU-free because the sibling track (`../TF_cobinding_multiome_chrombpnet`,
ChromBPNet — the method the paper actually used for its co-binding claim) is blocked
on GPU access.

**Source paper:** Zamboni, Llorens-Bobadilla et al., *The regulatory code of injury-
responsive enhancers enables precision cell-state targeting in the CNS*, Nature
Neuroscience 29:337–349 (2026), doi:10.1038/s41593-025-02131-w.

---

## Status

| Stage | State |
|---|---|
| PEAKQC cell filtering | done — 96.2% of 67,072 cells pass |
| Pseudobulk BAMs | done — 16, synthesised from fragments, verified lossless |
| Tn5 bias gate | done — 16/16 pass across a 140× depth range |
| ATACorrect → ScoreBigwig | done — 16 bias-corrected footprint tracks |
| BINDetect | done — 8 groups × 746 motifs, **36.2 M bound sites** |
| TF-COMB co-binding (Tier 3) | done — **FAILED as pre-registered**, cause diagnosed |
| **The lab's own data** | **not started — this is the actual goal** |

### What reproduces, and what does not

**Reproduces.** AP-1 is the dominant injury factor in glia (the paper's headline).
Secondary factors TFEC, BHLHE40, HIF1A, TFEB, CEBPD, STAT3 all positive. Aggregate
footprints recover Extended Data Fig. 6b (FOS 5.26× injured/uninjured, CTCF invariant
at 0.85×). C/EBP is sharply microglia-specific — found twice, by two independent
analyses.

**Does not reproduce.** The cell-type-specific *co-binding* mechanism (Tier 3:
1 of 4 rows, 3 required). The cause is diagnosed and is not a pipeline bug: the top
of an AP-1 co-occurrence ranking is dominated by motifs that **resemble AP-1**
(Smad2::Smad3 at motif distance 0.06, NFE2 0.24, BACH1 0.41) — the same DNA matched
by two similar matrices. The paper avoids this by learning motifs *de novo* from a
neural network, weighting by contribution score, and filtering to TFs expressed in
≥1% of cells in that cell type. We do none of the three. See `VALIDATION_CRITERIA.md`
→ "RESULTS AS JUDGED — TIER 3".

**Found something the paper missed.** Two bad libraries. Their nucleus QC used counts
and mitochondrial fraction only; the fragment-periodicity check catches what that cannot.

---

## What is in this repo — and what is not

The repo holds **code, criteria and figures only**. All data is derived and gitignored:

| Not in git | Size | How to regenerate |
|---|---|---|
| `bam/` | 35 G | `pipeline/build_all_bams.sh` |
| `tobias/atacorrect/`, `tobias/footprints/` | 32 G | `pipeline/run_tobias_chain.sh` |
| `tobias/bindetect/` | 61 G | `pipeline/run_bindetect.sh` |
| `tobias/tfcomb/` | 1 G | `pipeline/prep_tfcomb_input.sh` + `run_tfcomb_tier3.py` |
| `qc/`, `meta/peakqc_scores.csv` | 223 M | `pipeline/run_peakqc.py` (~90 min) |
| `meta/*.fa`, `meta/peaks_tobias/` | 2.7 G | downloaded / derived |

**The daily working notes (`<date>.txt`) are gitignored too.** They are the reasoning
record and they are **single-copy** — they exist only on the working machine and the
file server. Because of that, **commit messages and PR bodies are the durable home for
methods evidence in this project**; they are written long deliberately. Keep the habit.

---

## Environments

Three, and they are separate on purpose.

| Env | Used for | Build with |
|---|---|---|
| `tobias` (conda) | TOBIAS, PEAKQC, everything up to BINDetect | existing; **never repair in place** — it produced every result |
| `bedtools` (conda) | `bedtools`, required by BINDetect `--output-peaks` | existing |
| `~/venvs/tfcomb` | TF-COMB / Tier 3 | `bash pipeline/setup_tfcomb_env.sh` |
| `zamboni-r` (conda) | Seurat/Signac R work | existing |

`tfcomb` is deliberately a separate **native arm64** venv. The `tobias` conda env is
x86_64 under Rosetta, and TF-COMB cannot install there — llvmlite ships an arm64 macOS
wheel and no x86_64 one. The setup script asserts the architecture before installing;
see its header for the full diagnosis.

---

## Running it

Scripts take **overridable environment variables** rather than edits. Set them:

    PY=/path/to/python            # default $HOME/miniconda3/envs/tobias/bin/python
    CB=../TF_cobinding_multiome_chrombpnet   # genome, blacklist, peaks live here
    BEDTOOLS_BIN=/path/to/bedtools/bin
    MOTIFS=meta/motifs/JASPAR2020_CORE_vertebrates_non-redundant.meme
    PARALLEL=2  CORES=6

Order:

    python pipeline/run_peakqc.py              # cell QC          ~90 min
    bash   pipeline/build_all_bams.sh          # pseudobulk BAMs
    python pipeline/check_atacbias.py          # Tn5 bias gate -- run BEFORE the chain
    bash   pipeline/run_tobias_chain.sh        # ATACorrect -> ScoreBigwig
    bash   pipeline/run_bindetect.sh           # differential TF binding
    bash   pipeline/prep_tfcomb_input.sh <celltype> <arm>
    ~/venvs/tfcomb/bin/python pipeline/run_tfcomb_tier3.py
    ~/venvs/tfcomb/bin/python pipeline/diagnose_tier3.py
    ~/venvs/tfcomb/bin/python pipeline/plot_tier3.py

---

## Things that will bite you

- **TOBIAS cannot be invoked directly on macOS.** Its multiprocessing assumes `fork`;
  macOS defaults to `spawn` and every run dies on `logger cannot be pickled`. Always go
  through `pipeline/tobias_fork.py`. The same applies to TF-COMB — run it from a real
  `.py` file with an `if __name__ == "__main__":` guard, never a heredoc, or the workers
  re-import `__main__` and fork endlessly.
- **Do not re-apply the Tn5 shift.** The fragment files are already shifted. Applying it
  again in `fragments_to_bam.py` would silently corrupt every footprint.
- **`bedtools` must be on PATH** for BINDetect's `--output-peaks`.
- **Never build a partial motif set for TF-COMB.** It ranks each pair against the
  background of all others, so dropping motifs changes the scores of the rows that remain.
  `prep_tfcomb_input.sh` refuses anything but 746.
- **`*_all.bed` is 16 GB and PlotAggregate reads it.** It is excluded from the file-server
  copy, so remaking the aggregate figures needs it fetched first.
- **Output existing is not the same as output being complete.** Check timestamps and
  logs, not just presence — a background job's "exit 0" can be the wrapper's `echo`, not
  the work.

### Known issue

`pipeline/build_all_bams.sh:87` uses `[ -s bam/${tag}.bam ]` as its completion test. A
truncated BAM is non-empty, so an interrupted build is treated as done. Carried since
2026-08-14; more relevant now, because lab data means new BAM builds.

---

## Validation

`VALIDATION_CRITERIA.md` is the contract. Its rule: **a criterion may be specified or
tightened only before the data it judges exists**; criteria already judged are never
rewritten, a critique is appended instead. Verdicts on record: Tier 1 detected but was
demoted for family-size bias; Tier 2 passed; Tier 3 **failed**; Tier 4 was retired on
its premise by the supervisor; Tier 5 was 1 of 3 with the thresholds themselves recorded
as miscalibrated; Tier 0 (null control) is specified but deprioritised.

A revised Tier 3 and an RNA expression filter are both **specified and deliberately not
run** — writing a test after seeing the numbers makes it a hypothesis, not a test.

---

## Layout

    pipeline/   all scripts; headers carry the reasoning, read them first
    meta/       motifs, peaks, genome (largely gitignored)
    figures/    committed figures (PDF + PNG)
    VALIDATION_CRITERIA.md   pre-registered criteria and every verdict
    <date>.txt  daily working notes -- gitignored, single-copy
