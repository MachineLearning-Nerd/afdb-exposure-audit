# e428 results — exposure audit of two published "unseen" evaluations

> Registration: `docs/e428_registration.md` (OTS-stamped, committed 7c6be201
> before any query). Run 2026-10-09 11:10 UTC, exit 0. Code:
> `experiments/e428_exposure_audit.py`. Outputs: `results/e428/`
> (`t_structures.json`, `a_chains.json`, `stats.json`, `receipts.json` with a
> SHA-256 per request, `run.log`). 2,258 request receipts (25 are cache re-reads) over 2,233
> unique cached responses: 200 ×1,872, 204 ×340 (empty search result),
> 404 ×21. Audit: `docs/audit_e428.md`, PASS-with-notes; its notes are
> corrected in this document.

## Bottom line

- **T (Terwilliger et al. 2024): primary prediction NOT supported.** Structures
  with a close training-era relative do not have clearly lower published
  r.m.s.d. P1 alone is uninformative (CI about ±0.5 Å); what
  shows that exposure does not explain the low r.m.s.d. is that the 49
  NOVEL structures have a median of 0.943 Å, close to the overall 0.954 Å. This is a negative result for our hypothesis and does not
  challenge their conclusion.
- **A (AlphaFlow ATLAS test split): descriptive only.** 32 of 82 test chains
  (39 %) have a close relative (same UniProt accession or ≥ 95 % identity) in
  PDB entries released by AlphaFold2's cutoff. Three of those 32 are fusion
  constructs; for two of them (6xrx_A, 7p41_D) the close status comes only
  from the fusion partner (see caveat 2), so the share is 30–32 of 82
  (37–39 %). Five test chains have a ≥ 95 % relative
  in AlphaFlow's own training split; in three of them the match is through a
  fusion partner (MBP or T4 lysozyme), leaving two genuine same-protein or
  near-identical-ortholog overlaps. Forty test entries were released on or
  before ESMFold's stated training cutoff (2020-05-01). No performance claim is
  made.

## T — Terwilliger et al. 2024 (n = 102)

Exposure classes (registered rule, first protein entity):

| class | n |
|---|---|
| SAME_ACCESSION | 9 |
| SEQ95 | 15 |
| SEQ30 | 28 |
| NOVEL | 49 |
| QUERY_FAILED | 1 (7DRH) |

Close relative = 24 of the 101 classed structures (24 %). Two structures had no
UniProt reference (7DTR, 7S5L) and were classed from sequence search only (both
NOVEL).

**7DRH.** The RCSB entry endpoint returns 404 (cached). The entry was removed
from the PDB on 2023-09-13 and replaced by 8ITE: RCSB holdings/removed record,
fetched after the run (2026-10-09T11:18Z), cached at
`data/e428/holdings/7DRH_removed.json`, SHA-256 `631f4664…d652`. As registered, it is QUERY_FAILED and excluded from tests. No
substitution of 8ITE was made.

**T-P1 (primary).** Published Cα r.m.s.d. (PRED-RMSDLIST!AR):

| | n | median (Å) |
|---|---|---|
| close relative | 24 | 0.847 |
| other | 77 | 0.956 |

Difference in medians (close − other) −0.109 Å, 95 % bootstrap CI
[−0.412, 0.534] (10,000 resamples within groups, seed 0), two-sided permutation
p = 0.376 (10,000, seed 1). The registered support criterion (CI upper bound
< 0) is **not met**. The direction matches the prediction, but the interval is
wide and includes differences in both directions.

**Secondary (exploratory).**
- S2, mean pLDDT: close 93.19 vs other 93.95; difference −0.76, CI
  [−2.89, 0.91], p = 0.320. No clear difference in confidence.
- S3, OLS of log r.m.s.d. on close (0/1) and mean pLDDT (n = 101):
  coefficient of close −0.254, 95 % case-bootstrap CI [−0.532, 0.011]
  (10,000, seed 2); pLDDT coefficient −0.072 per point. Borderline; a
  close-relative structure has about 22 % lower r.m.s.d. at equal pLDDT, but
  the interval includes zero. Exploratory, one of several secondary analyses,
  not corrected for multiplicity: we do not report it as an effect.
- S4, median r.m.s.d. by class: SAME_ACCESSION 0.757 (n 9), SEQ95 0.873
  (n 15), SEQ30 1.251 (n 28), NOVEL 0.943 (n 49). Not monotone in exposure:
  the NOVEL group is more accurate than the SEQ30 group.
- S1, shares among the 215 structures with a mean pLDDT: **not computed.** The
  script classed only the 102 structures of the r.m.s.d. list. This is a
  deviation from the registration; it was descriptive only and is not needed
  for any conclusion.

**Reading.** In this template-free sample of high-confidence predictions,
having a close relative in the training-era PDB is not associated with a clear
r.m.s.d. advantage, and structures with no relative at all (NOVEL, n 49,
median 0.943 Å) are as accurate as the sample as a whole (0.954 Å). Terwilliger et al.'s own explanation (high confidence of
the selected predictions) is not contradicted by our analysis.

## A — AlphaFlow ATLAS test split (n = 82)

| training-era class (≤ 2018-04-30) | n |
|---|---|
| SAME_ACCESSION | 28 |
| SEQ95 | 4 |
| SEQ30 | 25 |
| NOVEL | 25 |

Close training-era relative: 32 of 82 (39 %). Seven chains had no SIFTS UniProt
mapping (6sms_A, 6okd_C, 6l8s_A, 6d7y_B, 6ypi_A, 7jrq_A, 5znj_A) and were
classed from sequence search only.

Overlap with AlphaFlow's own training split (1,266 train chains; 1,249 with a
SIFTS accession; no date filter):

| relation | test chains |
|---|---|
| same UniProt accession | 3 |
| ≥ 95 % identity hit | 5 |
| ≥ 30 % identity hit | 10 |

The five ≥ 95 % chains, inspected after the run:

| test chain | accessions | train match | reading |
|---|---|---|---|
| 6xds_A | P0AEX9 (MBP), Q9NZC2 (TREM2) | 5hz7_A, 2xz3 | MBP fusion; train-split match through MBP |
| 6xrx_A | P0AEX9 (MBP), Q8T5C5 | 5hz7_A, 2xz3 | MBP fusion; train-split match through MBP |
| 7p41_D | D9IEF7 (T4 lysozyme), Q5VT66 (MARC1) | 4epi | T4 lysozyme fusion; match through T4 lysozyme |
| 7e2s_A | P72583 | 1y6i_A | same protein |
| 7buy_A | P0DTD1 (SARS-CoV-2 Mpro) | 2h2z | SARS-CoV Mpro, ≥ 95 % identical ortholog |

Released on or before 2020-05-01 (ESMFold training cutoff as stated by
AlphaFlow): 40 of 82, re-derived from `release_date` (matches the scoping read).

**Reading.** AlphaFlow's test split is a date split of ATLAS with no homology
filter. Measured against AlphaFold2's training-era PDB, about a third of its
test proteins are the same protein or a ≥ 95 % relative of something that was
available for training. Against AlphaFlow's own MD training split, genuine
same-protein overlap is small (2 of 82). This describes the split; it does not
show that AlphaFlow's reported metrics are inflated, since no per-target scores
were analysed (registered guard).

## Source check (2026-10-09, after the audit, before use in the paper)

- Terwilliger et al. text (Europe PMC PMC10776388): "the median r.m.s.d. is
  1.0 Å" is stated in a paragraph that says "we used all 215 structures".
  The workbook's `NumbersUsedInPaper!A138` = `PRED-RMSDLIST!AR4` =
  `MEDIAN(AR8:AR134)` = 0.954 ("the median rmsd is 1. Å"), i.e. the reported
  value is the median of the 102-row list we used. Verbatim quotes confirmed:
  "perhaps due to the high confidence in prediction in our sample"; "only the
  first chain was included"; "initial AlphaFold predictions (made without
  templates)"; December 2021 to May 2022; free R ≤ 0.30 → 102 structures.
- AlphaFlow (arXiv 2402.04845 PDF): "Using training and validation cutoffs of
  May 1, 2018 and May 1, 2019, we obtain train/val/test splits of 1265/39/82
  ensembles"; test = "all 84 targets whose corresponding PDB entries were
  deposited after May 1, 2019, minus the two targets with sequence length
  greater than 1024"; "These test proteins have minimal structural overlap
  with the training ensembles, providing a stringent test of generalization";
  ESMFold cutoff May 1, 2020. The released `atlas_train.csv` has 1,266 rows
  (paper: 1,265).

## Caveats

1. **Association only (T).** Close-relative structures may differ in size,
   family and crystallisability. The negative result holds for this sample of
   101 structures and for the registered close/other split; it is low-powered
   for small effects (CI width ≈ 0.95 Å).
2. **Fusion constructs (A, post hoc).** The registered rule classes a chain by
   any of its UniProt accessions, so a fusion tag (MBP, T4 lysozyme) can make a
   chain SAME_ACCESSION or SEQ95. Three A chains have two accessions, all three
   fusions. 6xds_A is SAME_ACCESSION through its protein of interest too
   (TREM2, Q9NZC2: 3 training-era entries); 6xrx_A (Q8T5C5) and 7p41_D (MARC1,
   Q5VT66, no hits) are close only through the tag. The 30–32 range is a post
   hoc bound, not a re-run with a different rule. All five ≥ 95 % train-split
   matches of the three fusions go through the tag. No T chain has more than one
   accession.
3. **T query chain.** "First chain" was operationalised as the lowest-id
   protein entity; if Terwilliger et al. used a different chain in a
   heteromeric entry, its class could differ.
4. **Count note.** The registration's disclosure gave 1,265 train rows; the file
   has 1,266 data rows (1,267 lines with header). The scoping read was off by
   one; all analyses use the file.
5. **Search tool.** Classes depend on RCSB's MMseqs2 search (E ≤ 0.1, no
   coverage threshold) at query time (2026-10-09); a short local match can
   create a SEQ30 hit.

## What would falsify these readings

- T: a larger or differently drawn template-free benchmark showing a clear
  r.m.s.d. advantage for close-relative targets would show our negative result
  was a power failure.
- A: if the protein of interest of 6xrx_A had a training-era relative missed
  by the sequence search, the share would be 31/82; if a structure-based (not sequence-based)
  comparison found more overlap with the ATLAS training ensembles, the
  "minimal overlap" description would be weaker than our counts suggest.
