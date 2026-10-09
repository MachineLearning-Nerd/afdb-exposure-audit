# e428 registration — exposure audit of two published "unseen" evaluations

> Status: REGISTERED 2026-10-09, written BEFORE any exposure query was sent and
> before any exposure class or exposure–performance association was computed.
> E-number allocated by main (single-session mode). Requested by the human
> ("do option 1, check a published result"). Descriptive and exploratory
> re-analysis of public data; no AFDB/e422/e427 census data are used, and no
> blinded e422/e427 quantity is touched.

## 1. Question

Do two published evaluations that treat post-cutoff PDB structures as new
contain proteins with close relatives in AlphaFold2's training-era PDB, and
(target T) is the published accuracy associated with that exposure?

- **T — Terwilliger et al. 2024, Nat. Methods 21:110** (doi
  10.1038/s41592-023-02087-4). 102 crystal structures released Dec 2021 – mid
  2022 (R_free ≤ 0.30), AlphaFold predictions made without PDB templates, first
  chain of each entry. Headline: median Cα r.m.s.d. of prediction vs deposited
  model 1.0 Å, which the paper contrasts with DeepMind's 2.3 Å and attributes
  "perhaps" to "the high confidence in prediction in our sample". The paper does
  not discuss homology to AlphaFold's training data.
- **A — AlphaFlow/ESMFlow, Jing et al., ICML 2024** (arXiv 2402.04845). ATLAS
  test split of 82 chains deposited after 2019-05-01, described as having
  "minimal structural overlap with the training ensembles, providing a
  stringent test of generalization"; no homology filter. Fine-tuned from
  AlphaFold2 (PDB ≤ 2018-04-30) and ESMFold weights.

## 2. What was seen before registration (disclosure)

- T: the workbook's sheet/column layout and the median of column AR (0.954 Å
  over the 102 rows), confirmed by a scoping read; the paper's ID list equals
  `PRED-RMSDLIST!AP8:AP109`. No exposure information for any T structure.
- A: split files downloaded (row counts 82/1265/39). A scoping read reported
  that 40 of the 82 test entries have release dates before 2020-05-01; this
  count is re-derived here, not assumed. No exposure information computed.

## 3. Data (cached under `data/e428/`, hashes in results)

- T: `data/e428/terwilliger/AlphaFoldCrystal_M.xlsx`, SHA-256
  `7b3e1473b99e1575019d8c095d292d6a96de702c36580241cf85932e7d69d143`
  (https://phenix-online.org/phenix_data/terwilliger/alphafold_crystallography_2022/,
  downloaded 2026-10-09). Outcome: sheet `PRED-RMSDLIST`, rows 8–109, ID column
  AP, r.m.s.d. column AR. Mean pLDDT: sheet `Summary`, ID column B, column H,
  joined by PDB ID. The workbook has no licence statement: it is not
  redistributed; only derived per-structure exposure classes and the IDs are
  released, with the published r.m.s.d. values cited to the source.
- A: `data/e428/alphaflow/atlas_{train,val,test}.csv` from
  github.com/bjing2016/alphaflow at commit `0408d7c8` (MIT), SHA-256 in
  `data/e428/alphaflow/SHA256SUMS`.
- Lookups (cached request-keyed, SHA-256 of every response): RCSB Data API
  (entry/entity sequences, UniProt references, release dates), RCSB Search API
  (MMseqs2 sequence search and attribute queries), PDBe SIFTS mappings API.

## 4. Method (frozen)

**Query chain and sequence.**
- T: the paper uses "the first chain". Operational rule: the protein polymer
  entity with the lowest entity id in the RCSB entry; query sequence = its
  canonical one-letter sequence; accession(s) = its UniProt reference
  identifiers in RCSB.
- A: chain named in `atlas_test.csv` (`<pdb>_<authchain>`); query sequence =
  `seqres` column; accession = PDBe SIFTS mapping of that author chain
  (UniProt accession(s) mapped to `chain_id`).

**Training-era exposure classes** (identical to the paper's V2 rule, query =
the chain's sequence instead of the canonical UniProt sequence): RCSB entries
with initial release ≤ 2018-04-30 that (i) contain a polymer entity with the
same UniProt accession, (ii) contain a chain at ≥ 95 % sequence identity or
(iii) at ≥ 30 % identity (MMseqs2 via RCSB Search, E-value ≤ 0.1, no coverage
threshold). Class = strongest of SAME_ACCESSION > SEQ95 > SEQ30 > NOVEL. A chain
with no UniProt accession is classed from (ii)–(iii) only and flagged. Any
failed query → QUERY_FAILED (excluded from tests, reported). **Close relative**
= SAME_ACCESSION or SEQ95.

**A only — overlap with AlphaFlow's own training split.** For each test chain:
whether any `atlas_train.csv` entry (by 4-letter PDB ID) is among (a) RCSB
≥ 95 % or (b) ≥ 30 % sequence hits with no date filter, or (c) shares a
UniProt accession (PDBe SIFTS for train chains). And the count of test entries
with `release_date` ≤ 2020-05-01 (ESMFold's training cutoff as stated in the
AlphaFlow paper; quoted in the results).

## 5. Analyses and predictions

**T-P1 (primary).** Among the 102 structures, compare the published Cα
r.m.s.d. (AR) between close-relative and other structures: difference in
medians (close − other), 95 % percentile bootstrap CI (10,000 resamples within
groups, seed 0) and two-sided permutation p (10,000 label permutations, seed
1). Prediction: median r.m.s.d. lower for close-relative structures.
Supported if the CI upper bound < 0. If either group has < 10 structures, report
descriptively and do not test.

**T-S (secondary, exploratory).**
1. Class shares among the 102 (and among the 215 for which a mean pLDDT exists,
   descriptive only).
2. Same comparison as P1 for mean pLDDT (Summary!H), prediction higher for
   close-relative structures.
3. Whether exposure explains r.m.s.d. beyond pLDDT: OLS of log r.m.s.d. on
   close (0/1) and mean pLDDT, coefficient of close with bootstrap CI (10,000
   case resamples, seed 2). Exploratory.
4. Medians by the four classes.

**A-P1 (descriptive).** Numbers and shares of the 82 test chains: by
training-era class; with a ≥ 95 % / ≥ 30 % / same-accession relative in
AlphaFlow's own training split; released ≤ 2020-05-01. No performance analysis
here (per-target scores are not published; recomputing them from released
ensembles would need its own registration).

## 6. Interpretation guards

- Association, not causation: structures with close relatives may also differ
  in size, family, crystallisability and how much they have been studied.
- A training-era relative shows the protein or a homologue was available to
  training, not that it was sampled, nor that it explains accuracy.
- T used no templates, so the template channel does not apply; A-P1 shares are
  descriptive and do not show that AlphaFlow's reported metrics are inflated.
- If T-P1 is not supported, that is reported as a negative result: exposure
  would then not explain the paper's low r.m.s.d. in this sample.

## 7. Compute and outputs

About (102 + 82) × 3 search queries, ~1,500 SIFTS/RCSB lookups, rate-limited
(≥ 250 ms between requests), all cached; CPU only, minutes. Code:
`experiments/e428_exposure_audit.py`. Outputs: `results/e428/` (per-structure
classes, statistics JSON, request receipts) and `docs/e428_results.md`. Audit:
independent recomputation of the T-P1 statistic and of a random 20 % of
classes from the cached responses, `docs/audit_e428.md`.
