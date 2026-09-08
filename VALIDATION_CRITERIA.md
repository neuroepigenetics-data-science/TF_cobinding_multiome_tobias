# Validation criteria — TOBIAS track vs. Zamboni et al.

**Written 2026-08-20 at 14:58, BEFORE any BINDetect result existed.** The run was
launched at 14:58:51; the first results are not expected for ~4h. This file is
deliberately timestamped and committed so that the criteria cannot be adjusted
after seeing the output.

## Why this file exists

The purpose of the TOBIAS track is not primarily to discover new biology. It is
to establish **whether this GPU-free pipeline reproduces a known answer**, so
that it can be trusted on our own lab data where no published comparison exists.

That goal has a hard requirement attached: *a validation only validates if it
could have failed.* If we decide what counts as agreement after seeing the
results, we learn nothing transferable — we will simply have described our own
output. So the pass/fail conditions are fixed here, in advance.

Reference: Zamboni, Llorens-Bobadilla et al., *Nat Neurosci* 29:337–349 (2026),
DOI 10.1038/s41593-025-02131-w. Local copy:
`../Literature/s41593-025-02131-w/s41593-025-02131-w.md`.

## AMENDMENT LOG

Amendments are listed here so the document's history is legible. The rule: a
criterion may be **specified or tightened only before the data it judges
exists**. Criteria already judged against results are NEVER rewritten — a
critique is appended instead, and the original verdict stands.

- **2026-08-20 ~15:05** — divergence 2 reworded ("we have no RNA" was wrong;
  this is a multiome dataset). Before any result existed.
- **2026-08-26 — SUPERVISOR DECISIONS. These change SCOPE and PRIORITY, and
  resolve two open items on scientific grounds. They do not rewrite any judged
  criterion.** See "SUPERVISOR RULINGS" below.
- **2026-08-24** — Tiers 1 and 2 have already been judged, so their criteria are
  left **verbatim** and a POST-HOC CRITIQUE is appended below them. Tier 3 is
  **rewritten**, which is legitimate because TF-COMB has not been installed or
  run. Tiers 0, 4 and 5 and the statistical policy are **new**. Tier 0 is
  genuine pre-registration (its data does not exist). Tiers 4 and 5 are
  **specified-before-inspection**, a weaker standard: they run on results that
  already exist, but the criteria below were written and committed before those
  numbers were looked at. That distinction is stated wherever it applies and
  must survive into any write-up.

---

## What we are comparing against, and what we are NOT

We are comparable to the paper's **supervised motif enrichment** (Fig. 4a),
which used *the same motif set we are using* — JASPAR2020 CORE vertebrates, via
Signac `AddMotifs()`. That is the deliberate reason for the motif-set choice.

We are **only loosely** comparable to their primary deep-learning analysis
(ChromBPNet → TF-MoDISco → tangermeme), which discovered motifs *de novo* and
then matched them to a database. Different procedure, different failure modes.

We are **not** comparable to their SCENIC+ enhancer-GRN analysis at all; that
uses the RNA modality, which this pipeline does not touch.

---

## TIER 1 — a detection check (NOT load-bearing; see critique below)

> **2026-08-24 ANNOTATION — the criterion text below is unchanged and its
> verdict stands.** It was originally headed "the load-bearing check". That was
> wrong: Tier 1 turns out to have no cell-type specificity, so it can establish
> that the pipeline *detects the dominant injury factor* but cannot establish
> that it recovered the paper's *cell-type-specific* biology. It is necessary
> but not sufficient. See POST-HOC CRITIQUE OF TIERS 1-2. The load-bearing
> check is now Tier 3.

**AP-1 must come out positive (up after injury) in glial cell types.**

| Cell type | Expectation |
|---|---|
| Astrocytes | AP-1 strongly positive |
| Oligodendrocytes | AP-1 strongly positive |
| OPCs | AP-1 positive |
| Microglia | AP-1 positive but **weakest of the four** |

AP-1 family = FOS, FOSB, FOSL1, FOSL2, JUN, JUNB, JUND, ATF3, and the dimer
matrices (FOS::JUN, FOSL1::JUNB, …).

Paper basis: AP-1 dominance in glial IRENs was found by **three independent
methods** (supervised motif enrichment, ChromBPNet/MoDISco, chromVAR). It is the
strongest and best-replicated finding in the paper.

Microglia being weakest is a **prediction, not a hedge**: the paper reports
microglia IRENs as less AP-1 reliant, driven instead by TFEC, BHLHE40 and HIF1a.
A weak microglial AP-1 is therefore a CORRECT result.

**Verdicts, fixed in advance:**

- **PASS** — AP-1 family among the top positive differential factors in ≥3 of
  {Astrocytes, Oligodendrocytes, OPCs, Microglia}, and microglia is the weakest
  of the four.
- **INVESTIGATE** — AP-1 positive but not top-ranked, or microglia not weakest.
  Report as a partial reproduction; do not present the pipeline as validated.
- **FAIL, suspect our pipeline first** — AP-1 *negative* (down after injury) in
  any glial type. Check in this order: (1) sign convention — confirm the column
  is literally `<ct>_I_<ct>_U_log2fc`; (2) peak set; (3) normalization. Do not
  conclude we have overturned a triple-replicated published result.

---

## TIER 2 — cell-type-specific secondary factors

Weaker checks; each rests on a single method in the paper, so individual misses
are tolerable. Expect these among the positive differential factors:

| Cell type | Expected secondary factors |
|---|---|
| Astrocytes | STAT3, CEBPD |
| OPCs | TFEB |
| Microglia | TFEC, BHLHE40, HIF1A ranking high |

**Not testable:** ependymal EGR1/KLF6. The paper's claim is about IRENs, i.e.
differential; we have only occupancy for Ependymal_I. Recorded here so its
absence is not later mistaken for a failure.

All 25 paper-named TFs were confirmed present in our JASPAR file before the run.

---

## POST-HOC CRITIQUE OF TIERS 1-2 — written 2026-08-24, AFTER the results

**This section is post-hoc by construction and must be labelled as such
wherever it is reported.** It does not change the Tier 1/2 verdicts above. It
records that the criteria used the wrong instrument, and — importantly — the
corrected instrument makes the Tier 1 result *stronger*, which is exactly why
substituting it silently would have been dishonest.

**1. Tier 1's statistic contradicted this document's own divergence 3.**
Divergence 3 says "evaluate at FAMILY level, never per gene", but the verdict
was read off AP-1's *best* rank — a per-motif statistic, and a max over 33
correlated family members. Measured:

    AP-1 family size          33 of 746 motifs
    AP-1 best rank            1
    random 33-motif set       best rank <= 64 in 95% of draws (median 16)

A family that size almost always has *someone* near the top. **"Best rank = 1"
is close to guaranteed and is not evidence.**

**2. The fair statistic is the family's MEDIAN rank against a size-matched
null, and the result survives it decisively:**

    cell type          AP-1 median rank    size-matched null
    Astrocytes                       19    373 (5th pct 272)   p < 1e-4
    Oligodendrocytes                 19                        p < 1e-4
    OPCs                             21                        p < 1e-4
    Microglia                        22                        p < 1e-4

**3. THE SERIOUS FINDING — Tier 1 has no specificity.**

    Neurons_V                        19                        p < 1e-4
    Neurons_D                        19                        p < 1e-4

AP-1's family median rank is 19-22 in **all six** cell types, glia and neurons
alike. A criterion that the negative control passes just as decisively cannot
license the conclusion this file was written to support. This is not a pipeline
failure — it agrees with the paper, which calls AP-1 the "generic stimulus
response element" — but the criteria were written as though AP-1's rank would
discriminate cell types, and it does not. Hence the demotion, and hence Tier 3
being rebuilt as an explicit discrimination test.

**4. The "microglia weakest" sub-clause: met, on a margin too small to
interpret.** By median rank microglia IS weakest (22 vs OPCs 21, Astrocytes 19,
Oligodendrocytes 19). But that range is compressed to 3 rank positions and does
not reproduce the paper's Gini gap (0.585 vs 0.759). Report as "technically met,
not a reproduction". The original criterion also failed to name a measure at
all, which is a defect in the pre-registration and is recorded as one.

**5. An untested confound, now tested and CLEAR.** Nothing in the original file
checked whether AP-1 tops the list merely because everything shifted. It does
not: median `change` across all 746 motifs is ~0 (-0.006 to +0.042 by cell
type), while AP-1's median `change` is 0.38-1.24 — 10-40x the global median.
This diagnostic is promoted to a standing requirement in Tier 5b.

**6. Tier 2 survives a proper null test.** The paper's secondary factors landed
at ranks 50, 53, 70, 75, 75, 120 of 746; median 72 against a uniform-null median
of 373, **p = 0.0019**. Tier 2's original wording ("ranking high") was undefined,
but the claim holds under a real test.

**7. Minor.** Tier 1's FAIL branch names the column `<ct>_I_<ct>_U_log2fc`; in
TOBIAS 0.17.5 the summary column is `<ct>_I_<ct>_U_change` (an effect size,
`bindetect.py:521`). The per-site `_log2fc` lives in the per-TF overview files.
Divergence 5 is now *confirmed* rather than suspected for the DARs (config
`ml.merged_obj`), though still only suspected for `ML/filtered_peaks` itself.

---

## TIER 3 — THE LOAD-BEARING CHECK: co-binding as a DISCRIMINATION test

**Rewritten 2026-08-24, before TF-COMB was installed or run — genuine
pre-registration.** The original wording ("expected partners appear among the
top AP-1 co-occurring pairs for >=3 of 4") is superseded because it had the same
defect Tier 1 turned out to have: no null, no defined "top", and **no
requirement that the partner be specific to the right cell type**. Asking only
"does SOX10 appear among oligodendrocyte AP-1 partners" would very likely pass
for the same non-discriminating reason AP-1 passes everywhere.

**Why this is now the load-bearing check.** The paper's mechanism is AP-1 *plus
a lineage-specific partner* — the partner is what makes an enhancer cell-type
specific. AP-1's own rank was never going to discriminate cell types, and Tier 1
showed empirically that it does not. Cell-type specificity will appear in the
co-binding analysis or nowhere.

### The statistic

For cell type *c* and TF *t*, take the AP-1<->*t* co-occurrence strength from
TF-COMB's `market_basket()` rules using its default `cosine` measure. Aggregate
at FAMILY level throughout, per divergence 3: AP-1 = the 33 matching motifs,
NFI = NFIA/B/C/X, RFX = RFX1/2/3/4, SOX9 and SOX10 as their single motifs.
Score = the **median rank** of the family's AP-1 pair rules within cell type *c*.

### Primary criterion — DIAGONAL DOMINANCE (a 4x4 confusion matrix)

Rows = cell types, columns = the paper's expected partner set for each.

| | NFI+SOX9 | NFI+RFX | SOX10+TCF4 | ELK3+CEBPA |
|---|---|---|---|---|
| **Astrocytes** | *diagonal* | | | |
| **Ependymal** | | *diagonal* | | |
| **Oligodendrocytes** | | | *diagonal* | |
| **Microglia** | | | | *diagonal* |

**PASS** = the diagonal entry is the best (lowest median rank) in its ROW for
**>=3 of the 4** cell types.
**FAIL, EXPLICITLY NAMED IN ADVANCE** = every partner set ranks highly in every
cell type. That is a *non-discriminating* result and counts as a FAILURE of
Tier 3 even though every expected partner "appears" — it is the exact pathology
Tier 1 exhibited, and it must not be reported as a pass.

### Secondary criterion — beat a null

Each diagonal partner's AP-1 co-occurrence median rank must beat a size-matched
random-motif null at p < 0.05 (same permutation procedure as the Tier 1
critique: draw *n* random motifs, take the median rank, 10,000 draws).

### Limitations stated IN ADVANCE, so they are not mistaken for pipeline failure

1. **NFI cannot separate astrocytes from ependymal** — it is the expected
   partner for BOTH rows. Those two rows are discriminated only by their
   non-shared members, SOX9 and RFX respectively.
2. **SOX9 vs SOX10 may not be separable at all.** They are SOX-family paralogues
   binding near-identical consensus sequences, and divergence 3 states plainly
   that footprinting cannot resolve such paralogues. So the
   Astrocytes-vs-Oligodendrocytes discrimination may fail *for motif-similarity
   reasons rather than pipeline reasons*. If it fails, check whether TCF4 (an
   E-box factor, unrelated to SOX) carries the oligodendrocyte row on its own.
3. **The cleanest discriminators are therefore RFX (ependymal), TCF4
   (oligodendrocyte) and ELK3/CEBPA (microglia)** — three families with no
   overlap with each other. Astrocytes is the weakest row by construction.
4. **Ependymal has no second arm**, so its row is computed from a
   single-condition `CombObj` (occupancy co-occurrence) rather than a
   `DiffCombObj`. That is a different quantity from the other three rows and
   its rank is comparable *within* the ependymal row only. The diagonal-dominance
   test is a within-row comparison, so this is legitimate — but it must be
   stated.

---

## TIER 0 — THE NULL CONTROL (genuine pre-registration; data does not yet exist)

**Added 2026-08-24. This is the strongest single control available to us, and if
it fails every other tier is void — so it should run BEFORE Tier 3.**

Nothing in the original file asked the most basic question: *does this pipeline
report a strong differential when there is no real difference?* Tiers 1-3 all
assume the answer is no. That assumption is testable and has never been tested.

### Design

Take **Astrocytes** — 16 injured samples, 53.2M injured fragments, and the cell
type with the largest real AP-1 effect (+1.394), so the contrast between "real
injury signal" and "no signal" is sharpest. Split its injured samples into two
disjoint arms, `Astrocytes_IA` and `Astrocytes_IB`, then run the full chain
(fragments -> BAM -> ATACorrect -> bias gate -> ScoreBigwig -> BINDetect) with
identical peaks and motifs, comparing IA against IB.

**Split by SAMPLE, never by random cell.** A random cell split breaks batch
structure and would understate variance, giving a falsely clean null. A
sample-level split carries real library-to-library variation, which is the noise
floor a new dataset would actually face. Balance the two arms on timepoint and
on cell count, assign with a fixed seed, and record the assignment.

### Criteria, fixed in advance

- **PASS** — AP-1 family median rank **> 150** of 746 (i.e. not even in the top
  20%; a true null expects ~373), AND fewer than 5 of 746 motifs reach
  |change| >= 0.3. For reference the smallest *real* AP-1 effect across the six
  cell types was +0.445 (OPCs).
- **INVESTIGATE** — AP-1 median rank between 100 and 150, or 5-20 motifs past
  the effect floor. Suggests a residual batch/depth artifact that inflates but
  does not manufacture the real result.
- **FAIL, AND TIERS 1-3 BECOME UNINTERPRETABLE** — AP-1 median rank < 100 in an
  injured-vs-injured comparison. That would mean the differential signal is
  driven by something other than injury, and no amount of agreement with the
  paper would rescue it.

**Implementation note:** `fragments_to_bam.py` currently selects on
`cluster_ids_injury`; this needs a keepset that selects on an explicit sample
list. That is new code and must be written before Tier 0 can run.

---

## TIER 4 — CELL-IDENTITY DISCRIMINATION (specified before inspection)

**Added 2026-08-24. These results already exist, so this is NOT pre-registration
— but the criteria below were written and committed before the numbers were
looked at.** Label it that way in any write-up.

**Why it is worth adding:** its ground truth is independent of the paper
entirely, and unlike Tier 1 it is a genuine discrimination test. It also probes
something more basic than the injury response — whether the pipeline recovers
cell identity at all. A pipeline that cannot place SPI1 in microglia has no
business being trusted on new data.

### The statistic

Within each cell type, rank all 746 motifs by **bound fraction**
(`<cond>_bound / total_tfbs`) in the injured arm — or the single arm for the
occupancy-only groups. Ranks are scale-free, which is what makes this legitimate
under divergence 4 (bound *counts* are not comparable across cell types, but
within-cell-type ranks are).

### Expected assignments

| Lineage TF(s) | Should rank highest in |
|---|---|
| SPI1, IRF8 | Microglia and/or Macrophages_I |
| SOX10, OLIG2 | Oligodendrocytes and/or OPCs |
| SOX9, NFI | Astrocytes |
| RFX1/2/3/4 | Ependymal_I |
| NEUROD1/2, MEF2C | Neurons_V and/or Neurons_D |

**Out-of-tissue negative controls** — verified present in our motif set, and
none has any business being active in spinal cord: **GATA1** (erythroid),
**HNF4A** (hepatocyte), **POU5F1/OCT4** (pluripotency).

### Criteria

- **PASS** = at least **4 of 5** lineage assignments land on the correct cell
  type or its lineage pair, AND all three out-of-tissue negatives fall outside
  the top 25% in every CNS cell type.
- **FAIL** = fewer than 3 correct, or any out-of-tissue negative in the top 10%
  anywhere. Either would mean the bound calls are not tracking real biology.

---

## TIER 5 — INVARIANT-FACTOR AND GLOBAL-SHIFT CONTROLS (specified before inspection)

**Added 2026-08-24, same standard as Tier 4 — criteria written before the
numbers were inspected.**

**5a. CTCF sensitivity.** CTCF has the canonical strongest and best-characterised
ATAC footprint of any factor. It must rank in the **top 10% by bound fraction in
every cell type**. Failing this means the footprinting lacks basic sensitivity
and nothing downstream is trustworthy.

**5b. CTCF invariance.** CTCF occupancy is largely condition-invariant, and the
paper itself treats CTCF as a *common* (non-injury) factor — Extended Data
Fig. 6b contrasts it with injury-induced FOS. So CTCF's |change| must fall
**below the median |change|** across all 746 motifs in every differential cell
type. A large CTCF differential would indicate the contrast is picking up
something global rather than injury biology.

**5c. Global-shift diagnostic** (promoted from the post-hoc critique into a
standing requirement). The **median `change` across all 746 motifs must lie
within ±0.1 of zero** in every cell type. If it does not, TOBIAS's
cross-condition quantile normalization has not done its job and every
differential rank is suspect.

**PASS = all three.** These are cheap — they read existing output and need no
re-run.

---

## STATISTICAL POLICY — binding on every tier

Added 2026-08-24, after Tier 1 was found to have used an unfair statistic.

1. **Family-level median rank against a size-matched permutation null is the
   default statistic** for any multi-motif family (AP-1, NFI, RFX, C/EBP, SOX).
   Draw *n* random motifs, take the median rank, 10,000 draws. **Never use
   "best rank"** — for a 33-member family a top-1 hit is close to guaranteed.
2. **Multiple testing.** 746 motifs x 6 cell types = 4,476 tests. Apply
   Benjamini-Hochberg across the 746 within each cell type and report q-values
   alongside ranks.
3. **Effect-size floor.** |change| < 0.1 is uninterpretable regardless of
   p-value. Do not narrate rank-300 noise.
4. **Every cross-cell-type claim must be a DISCRIMINATION test** — does the
   correct cell type win? — never a presence test. Tier 1 is the cautionary
   example: every expected factor "appeared", in every cell type, including the
   ones where the paper predicts it should not dominate.
5. **Neurons_V and Neurons_D are the standing negative control** for any claim
   of the form "this is glia-specific". Any such claim must be checked against
   them before it is made.

---

## EXPECTED DIVERGENCES — stated in advance so they are not mistaken for failures

These are known methodological differences. Each will produce apparent
disagreement with the paper that is **not** a pipeline error. They must be
stated in any write-up rather than discovered by a reviewer — and must not be
"fixed" by tuning.

1. **Different denominators.** The paper measured motif *enrichment in
   differentially accessible regions* (IRENs). We measure *footprint-based
   differential binding over the cell type's full peak set*. Different
   quantities. Only **rank agreement among top factors** is meaningful;
   score-level correlation is not expected and its absence proves nothing.

2. **They adjusted for TF expression; our pipeline currently does not.**
   *(Amended 2026-08-20 15:0x, still before any result existed — the original
   wording said "we have no RNA", which was wrong. This is a multiome dataset;
   the RNA is present. The accurate statement is that the TOBIAS pipeline
   consumes ATAC fragments only, so nothing in it is expression-aware today.)*

   The paper's enrichment was "adjusted according to the expression of the
   motif-linked transcription factor after injury", and for TF-MoDISco they
   retained only TFs "expressed in at least 1% of cells within the cluster of
   interest."

   This is **partially correctable, and the correction is a post-hoc filter on
   `<ct>_results.txt` — it does NOT require re-running BINDetect.** Choosing
   JASPAR for the scan therefore did not cost us the expression criterion; the
   two compose.

   **COMMITMENT MADE IN ADVANCE:** if we apply an expression filter, it is
   reported as a *pre-specified secondary analysis*, alongside the unfiltered
   result — never as a replacement for it. Deciding whether to filter *after*
   seeing whether AP-1 came out right would be exactly the tuning this file
   exists to prevent.

   Residual, not correctable: a TF that is expressed and footprint-bound but
   whose motif is shared with a paralogue cannot be resolved to the right gene
   (see item 3).

3. **Motif redundancy — evaluate at FAMILY level, never per gene.** Footprinting
   cannot distinguish FOS from FOSB/FOSL1/FOSL2, JUN from JUNB/JUND, or NFIA
   from NFIB/NFIC/NFIX: the motifs are near-identical. This is exactly why the
   paper writes "AP-1" and "NFI" rather than individual gene names. In our
   JASPAR set FOS appears in 21 entries and JUN in 27, because of the dimer
   matrices — so **any concentration measure (e.g. a Gini coefficient) is
   inflated relative to the paper's and is not numerically comparable.** Their
   Gini values (oligodendrocytes 0.759, astrocytes 0.698, ependymal 0.600,
   microglia 0.585) may be compared for *ordering only*, and even that is weak.

4. **Bound-site counts are not comparable across cell types.** The bound
   threshold is calibrated per sample against its own background, so a shallower
   sample yields a higher threshold and fewer bound sites. Compare *which* TFs
   are bound, never *how many* sites.

5. **The peak set predates our cell filtering.** `ML/filtered_peaks/*.bed` is
   dated 2026-08-04; PEAKQC did not exist until 2026-08-10, so the peaks were
   likely called on the authors' unfiltered cell set. Impact judged small — peaks
   define where we look, not what we measure inside — but it matters most for
   Ependymal, which draws on the affected libraries.

6. **No real BAM exists for this dataset**, so synthesised-vs-real can never be
   compared directly. Validation of the BAM synthesis is source-level reasoning
   plus internal consistency (16/16 bias-gate passes across a 140× depth range).
   State this upfront in methods; a reviewer will ask.

---

## SUPERVISOR RULINGS — 2026-08-26

Relayed from the supervisor in response to the 2026-08-25 report. Recorded
verbatim in substance, with what each one closes.

**1. THE NEURON RESULT IS NOT A PROBLEM — CLOSED.**
> "AP-1 binds basically everywhere and it pops up in many analyses."

This resolves the largest open finding from 2026-08-20. Our observation that
AP-1 is as strong in neurons as in glia is *expected biology*, not a pipeline
defect and not a contradiction of the paper. It is consistent with the paper's
own framing of AP-1 as the "generic stimulus response element", and with the
2026-08-24 finding that AP-1 sits at family median rank 19-22 in all six cell
types. **No further investigation.**

**2. AP-1 RANKING MAY NOT BE CELL-TYPE SPECIFIC ANYWAY.**
> "I'm not sure whether the AP-1 ranking would be cell-type specific... it could
> be biased by the cell numbers and sequencing depth/quality... there are more
> interesting aspects such as the co-factors."

Directionally what the Tier 1 critique concluded independently. One piece of our
own evidence to keep attached: the depth confound WAS tested on 2026-08-24 and
argues against it — correlation between the injured/uninjured depth ratio and
the AP-1 effect size is **r = -0.498**, i.e. NEGATIVE, where a depth artifact
would give a strong positive. Cell-number bias remains untested.

**3. TIER 4'S PREMISE WAS WRONG — the tier is retired, not merely failed.**
> "I would not expect that the lineage factors will have the strongest binding
> compared to other TFs... they just need to be present at important loci
> (probably many of them) but necessarily not with rank 1."

This is the substantive resolution of the Tier 4 failure, and it is a better
answer than the statistical one. Tier 4 asked whether lineage TFs *rank highest*
in their own cell type. **That is not what lineage factors do.** They need to be
present at the right loci, not to out-rank every other factor. So the tier tested
a proposition that was never true, independent of the statistic used.

The 2026-08-26 diagnostics show the statistic was *also* invalid — the
bound-fraction ranking is the same list in every cell type (mean Spearman
rho = 0.969 across all 15 cell-type pairs) and is 85-90% explained by motif GC
content (rho = -0.84 to -0.91). Both the premise and the instrument were wrong.

**TIER 4 IS RETIRED.** Its 2/5 failure stands on record but is NOT evidence
against the pipeline. Any future cell-identity check must ask about *presence at
the right loci*, not rank.

**4. THE TIER 0 NULL CONTROL IS DEPRIORITISED.**
> "Since your results look really really good (in comparison with the paper) I
> would not invest too much time in the negative control. I have the feeling
> that the approach works pretty well and you could proceed with TF-COMB."

Supervisor's call, and the cost/benefit is reasonable: TF-COMB is the actual
objective and Tier 0 is an overnight job plus new code.

**ONE TECHNICAL NOTE, RECORDED ONCE AND NOT PRESSED.** Tier 0 answers exactly
the concern raised in ruling 2 — whether depth and cell-number differences can
manufacture a differential. Its value is also *higher*, not lower, for the lab
data (ruling 5), because there no published answer exists to check against. It
remains specified above and can be run whenever wanted; nothing is lost by
deferring it.

**5. NEW SCOPE — THE LAB'S OWN DATA.**
> "Most exciting of course would be to proceed with our own data. I talked to
> Desi, she is able to identify the cells in the different groups of interest...
> She has the fragment files on her computer. You could ask her to give you first
> groups which are easy to compare and where we expect the strongest differences."

**The data are FRAGMENT FILES — exactly what this pipeline was built to consume.**
`pipeline/fragments_to_bam.py` exists precisely because the Zamboni dataset ships
no BAMs, and that converter is what makes the lab data tractable without a GPU.
The riskiest component of the project turns out to be the one that transfers.

**6. EPENDYMAL IS DROPPED — CLOSED.**
> "Don't worry about the ependymal cells. you have 6 other cell types where you
> can run the analysis."

The longest-open item in the project (since 2026-08-11) is closed by decision
rather than by resolution. Ependymal_I's occupancy results remain on disk and can
still feed a single-condition TF-COMB run, but no differential claim will be made
and the uninjured-material question is not being pursued.

---

## RESULTS AS JUDGED — Tiers 4 and 5, run 2026-08-24

Recorded immediately after running, against the criteria exactly as written
above. **No threshold was adjusted after seeing these numbers.**

**TIER 5 — 5a FAIL, 5b FAIL, 5c PASS.**
- *5a (CTCF top 10% by bound fraction):* fails in 4 of 8 groups. Actual
  percentiles 7.1%-24.4% — CTCF is in the **top quartile in all eight**, but the
  10% line was crossed only in Astrocytes, both neuron types and Ependymal_I.
  Microglia 18.4%, Macrophages_I 24.4%.
- *5b (CTCF |change| below median |change|):* fails in all 6. CTCF sits at the
  55th-86th percentile of |change|. **But in absolute terms CTCF's differential
  is 3.8-13.5x SMALLER than AP-1's** (e.g. Astrocytes: CTCF 0.118 vs AP-1 median
  1.242).
- *5c (global shift):* passes everywhere, median `change` -0.006 to +0.042.
- **Assessment:** the scientific concerns behind 5a and 5b appear satisfied —
  CTCF is detected in the top quartile everywhere, and is an order of magnitude
  less variable than AP-1. The *thresholds* were set without calibration data
  and are the likely cause of the failures. That is recorded as a defect in the
  criteria, NOT corrected retroactively. A recalibrated 5a/5b may be proposed as
  a revision; the failures above stand on record either way.

**TIER 4 — FAIL (2 of 5 lineage assignments correct; >=4 required).**

    SPI1/IRF8    -> Microglia          CORRECT
    SOX10/OLIG2  -> Oligodendrocytes   CORRECT
    SOX9/NFI     -> Ependymal_I        WRONG (expected Astrocytes)
    RFX          -> Neurons_D          WRONG (expected Ependymal_I)
    NEUROD/MEF2C -> Microglia          WRONG (expected Neurons)

    out-of-tissue negatives all clean: GATA1 75.6%, HNF4A 50.9%, POU5F1 73.9%
    (criterion: outside top 25% everywhere)

- **A design flaw in this tier, admitted:** the expectation table assigns NFI to
  Astrocytes alone, but the paper's own Fig. 4e lists **NFI for ependymal cells
  too**. Tier 3's limitation 1 states this explicitly and it was not carried
  across to Tier 4. So "SOX9/NFI -> Ependymal_I" may be *correct biology scored
  as wrong*.
- **The two that passed have distinctive motifs** (ETS/SPI1, SOX10); two of the
  three failures involve promiscuous or shared motifs — NEUROD is a generic
  bHLH E-box shared with, among others, BHLHE40, which is a genuine microglial
  injury factor. This is divergence 3 biting harder than anticipated.
- **The competing reading must not be dismissed:** if footprint occupancy does
  not track cell identity, that is a serious problem for using this pipeline on
  new data, and it would not have been caught by Tiers 1-3. The clean
  out-of-tissue negatives argue against wholesale promiscuity, but do not settle
  it.
- **Status: unresolved.** A revised Tier 4 restricted to distinctive-motif
  families, using a specificity ratio (own cell type vs mean of others) rather
  than percentile ranks, is proposed — to be run and reported as a clearly
  labelled revision, with this failure retained on record.

---

## EXECUTION ORDER

Tiers are numbered by logical dependency, not by the order they were written
(Tier 0 was added last; the amendment log explains why the file is not
reordered). Run them in this order:

1. **Tier 5** and **Tier 4** — cheap, read existing output, no re-run. Do these
   first because they are nearly free and either one failing would undermine
   everything else.
2. **Tier 0** — the null control. Needs new BAMs and a chain run (overnight).
   If it fails, Tiers 1-3 are void, so it should precede any write-up.
3. **Tier 3** — the load-bearing discrimination test. Needs `pip install tfcomb`.
4. Tiers 1 and 2 are already judged; see their verdicts and the post-hoc
   critique.

---

## Scope of this run

- **Differential (both arms):** Microglia, Oligodendrocytes, Neurons_V,
  Astrocytes, Neurons_D, OPCs.
- **Occupancy only (single arm, no log2fc):** Ependymal_I, Macrophages_I.
- **Excluded:** Perivascular, Endothelial — pooled across the injury axis and
  therefore comparable to nothing in the paper. Endothelial is retained as the
  pipeline's technical canary only.

---

## RESULTS AS JUDGED — TIER 3, run 2026-09-07

**Recorded immediately after running, against the criterion exactly as written
above. No threshold was adjusted after seeing these numbers, and the criterion
text is not edited — the critique below is appended, per the amendment rule.**

Run by `pipeline/run_tfcomb_tier3.py` (TF-COMB 1.1.1, `cosine`, `max_dist=100`),
inputs built by `pipeline/prep_tfcomb_input.sh` from all 746 motifs' bound sites,
2.3-3.5M sites per arm. Reproduced bit-identically on a second run.

### TIER 3 — FAIL (1 of 4 rows; >=3 required)

Median rank of AP-1 <-> partner-set pairs, out of ~277,800 unordered pairs:

| | NFI+SOX9 | NFI+RFX | SOX10+TCF4 | ELK3+CEBPA |
|---|---|---|---|---|
| **Astrocytes** | **11,711** | 13,973 | 19,848 | 23,716 |
| **Ependymal_I** | 27,667 | **31,900** | 45,870 | 82,109 |
| **Oligodendrocytes** | 12,663 | 14,312 | **32,753** | 69,793 |
| **Microglia** | 23,007 | 16,987 | 111,176 | **23,384** |

(bold = the diagonal, i.e. the paper's expectation for that row)

- Only **Astrocytes** puts its diagonal first. Ependymal, Oligodendrocytes and
  Microglia are all won by an off-diagonal column.
- **Secondary criterion: 1 of 4** beat a size-matched random null at p<0.05
  (Ependymal p=0.0029; Astrocytes p=0.165, Oligodendrocytes p=0.572,
  Microglia p=0.505). For Oligodendrocytes the diagonal is *worse* than random.
- The pre-registered "non-discriminating matrix" failure mode is **not** what
  happened — row spreads are 0.71-4.06, so rows do separate the columns. They
  separate them the *same way*: `NFI+SOX9` is best in 3 of 4 rows. This is a
  COLUMN effect where the tier tests for a row effect.

### Three divergences from the criterion as written, all decided before running

1. **AP-1 has 35 members, not 33.** The pre-registered regex
   (`FOS|JUN|ATF3|BATF|JDP2`) matches 35 of the 746 and all 35 are genuine
   AP-1/TRE-binding bZIPs. The definition was pre-registered; the count was a
   miscount. All 35 used.
2. **RFX held to RFX1/2/3/4** as written, though RFX5 and RFX7 are in the motif
   set — widening a pre-registered set after the fact is not allowed, even
   though including them was the option that might have helped the row.
3. **Ependymal_I retained**, on limitation 4's own terms (single-arm occupancy,
   within-row comparison only, no differential claim). Supervisor ruling 6
   dropped ependymal cells but explicitly permits a single-condition TF-COMB
   run. Dropping the row would also have removed RFX, named in limitation 3 as
   one of the three cleanest discriminators.

### POST-HOC CRITIQUE — why it failed (`pipeline/diagnose_tier3.py`)

**Not evidence that the pipeline lacks cell-type-specific co-binding.** Three
defects in the statistic are demonstrable, and specific signal survives all of
them.

**1. Pooling unequal families makes the statistic an NFI measurement.**
`NFI+SOX9` is 140 NFI pairs against 35 SOX9 pairs. Decomposed:

    Astrocytes:  NFI 11,925 | SOX9 2,436 | pooled 11,711  (pooled-NFI gap: -214)

SOX9 alone sits at the **0.9th percentile** in Astrocytes — the strongest single
family/cell-type combination in the whole table — and pooling buried it. This is
the same class of defect that demoted Tier 1: a statistic answering to family
size rather than biology.

**2. Ranking against all pairs is confounded, because AP-1 pairs are globally
top-ranked.** A *random* motif paired with AP-1 already lands near the 6th
percentile in the differential rows. Raw percentiles therefore flatter every
family, and only the permutation null is informative — which the criterion did
include, to its credit, as the secondary test.

**3. THE LARGEST EFFECT, AND A NEW FINDING: the top of the AP-1 co-occurrence
ranking is motif redundancy, not co-binding.** Using BINDetect's own
`<ct>_distances.txt` (0 = identical motif, ~1 = unrelated), the strongest
"partners" in Astrocytes and Oligodendrocytes are motifs that resemble AP-1:

    Smad2Smad3 0.06 | NFE2 0.24 | BACH1 0.41 | BACH2 0.49 | MAFK 0.53 | NFE2L1 0.58

against ~0.997-0.999 for every genuine lineage partner (SOX9, SOX10, NFI, RFX,
ELK3, TCF4). Median distance of the top 10 partners: **Astrocytes 0.553,
Oligodendrocytes 0.553** — those rows are largely reporting the same DNA matched
by two similar matrices. This is divergence 3 measured rather than assumed, and
it is a real limitation of motif-based co-binding, not a bug.

**WHAT SURVIVES, and it is not nothing.** Per-family permutation nulls
(10,000 draws, size-matched) put three of the paper's four expected partnerships
in the right cell type at nominal significance:

    Ependymal_I  NFI    obs 16,602  null 72,902  p = 0.0011   survives Bonferroni (28 tests)
    Microglia    CEBPA  obs    821  null ~23,000 p = 0.0093   0.3rd percentile
    Astrocytes   SOX9   obs  2,436  null ~16,700 p = 0.0388   0.9th percentile

and **Microglia is the one row whose top partners are genuinely distinct from
AP-1** (median distance 0.940): TEF 0.87, DBP 0.81, NFIL3 0.91, CEBPG 0.98 — a
PAR-bZIP/C/EBP signature, found without reference to any expectation. That
independently reproduces the 2026-08-26 differential-heatmap result that the
C/EBP block is sharply microglia-specific, from a different statistic on a
different quantity.

**Status: Tier 3 FAILS as pre-registered, and that verdict stands.** A revised
criterion — per-family rather than pooled, scored against the AP-1 permutation
null, and excluding partners within some motif-distance of AP-1 — is proposed,
to be run and reported as a clearly labelled revision with this failure retained
on record. It is NOT run here, because specifying it after seeing these numbers
makes it a hypothesis, not a test.
