# Paper-verification analyses V1–V3 (registration, written BEFORE running V2/V3)

> Status: REGISTERED 2026-10-07 (written before any V2 network request and
> before any V3 computation). Repository HEAD at registration:
> `6041d3bd9ee0f60677e4bbd8f7c13aade56251d9`.
> Requested by the human (2026-10-07) to correct `paper/arxiv/temporal_leakage.tex`.
> Scope: descriptive verification for the paper revision. **No e-number has
> been allocated**; these analyses do not bind, amend, or pre-empt any
> registered e420/e421/e422 outcome. `main` may allocate an e-number and
> re-home these files; provenance is this file plus the scripts beside it.

## Why

A review of the paper (2026-10-07) found that its temporal rule compares the
reference release date with the AFDB `modelCreatedDate`. That date is when the
prediction was run (or relabelled), not the AlphaFold2 training cutoff. AF2 was
trained on PDB structures released up to **2018-04-30** (Jumper et al. 2021;
AlphaFold technical note), while AFDB predictions ran template search against a
PDB copy whose cutoff is not stated for the proteome runs (Tunyasuvunakool et
al. 2021, Methods). The paper also calls post-2022 entries "leakage-free"
without checking whether the same protein or close homologs were in the PDB
before the cutoff (for example, P01112/HRAS and P06241/FYN are in the e421
census).

## Analyses

**V1 — channel-separated recount of the ATLAS audit (no new data).**
Input: `results/e420/label_run_receipt.json` (typed failures carry
per-entry `release=` and `modelCreatedDate=`) and
`results/e420/overflow_temporal_check.json`. For each accession, count:
(T) the *training channel*: every dated reference entry was released on or
before 2018-04-30; (P) the *template channel*: every dated reference entry was
released on or before its AFDB `modelCreatedDate` (the paper's original rule).
Report accession-level and entry-level counts for the admitted pool,
overflow, and union, plus the release-year distribution. The accession whose
AFDB child failed to decode (Q582G4) is reported separately.

**V2 — prior-PDB exposure of the post-snapshot proteins (new public data).**
Universe: the 64 accessions with ≥1 successful e421 census entry, and the 138
accessions in the e422 cumulative batch-1–4 ledger. For each accession, fetch
the UniProt canonical sequence (rest.uniprot.org) and query the RCSB Search
API (search.rcsb.org, public, no credentials) for PDB entries with
`initial_release_date ≤ D`, where D ∈ {2018-04-30 (training cutoff),
the accession's recorded AFDB modelCreatedDate (conservative template bound)}.
Exposure levels per accession and date:
- `SAME_ACCESSION`: an entry whose polymer entity maps to this UniProt accession;
- `SEQ95`: a polymer entity with ≥95% sequence identity (RCSB mmseqs2
  sequence service, e-value ≤ 0.1);
- `SEQ30`: a polymer entity with ≥30% sequence identity (same service);
- `NOVEL`: none of the above.
Exposure class = the strongest level hit. Every request is cached under
`data/paper_verification/` with URL, UTC time, status, and SHA-256. Failed
queries are typed `QUERY_FAILED` and the accession is reported as unresolved,
never imputed as NOVEL.

**V2b — the paper's own rule on the e421 census (added before any V2 run).**
`results/e421/qualifying_entries.json` records AFDB `modelCreatedDate` values of
2025-08-01 (107 entries), 2022-06-01 (49), 2025-07-03 (11), 2025-03-31 (1),
while the census universe admitted releases after 2022-06-02. Fetch each
qualifying entry's `rcsb_accession_info.initial_release_date` (RCSB Data
API, public) and report how many of the 168 entries — and of the 96 successful
entries, with their residue counts — satisfy release > modelCreatedDate.
Entries that fail are exposed to the template channel by the paper's own
definition.

**V3 — exposure-stratified re-analysis (exploratory, descriptive).**
Using V2's training-channel exposure class (D = 2018-04-30):
(a) e421: residue counts, the lDDT<0.60 error rate per pLDDT bin, and Spearman
pLDDT–lDDT, split into `NOVEL ∪ SEQ30` (no ≥95% pre-cutoff relative) vs
`SAME_ACCESSION ∪ SEQ95`;
(b) e422 (cumulative batch-1–4 rows): the **pooled** split-conformal lDDT
interval exactly as registered in `docs/e422_registration.md` §3A.1 (median
centre fit on train; |y − ŷ| scores; ⌈(n+1)(1−α)⌉ order statistic), re-implemented
independently (no import of `experiments/e422_protocol.py`) from the
registration text: canonical sorted accession list over eligible rows;
`np.random.default_rng(seed).permutation(n)`; t = n//2, c = (n−t)//2;
train / calibration / test = perm[:t] / perm[t:t+c] / perm[t+c:]; **R = 200**
seeds 0–199; coverage reported overall (residue-weighted, as registered) and
by test-protein exposure class, plus the protein-averaged coverage. This is a re-implementation for verification of the
paper's citation of P5. It is **not** the registered P5 evaluator, binds
nothing, and does not compute any P6/P7 quantity (no per-pLDDT-bin coverage,
no Mondrian/CQR comparison) so as not to pre-empt those registered analyses.

## Decision rules (fixed now)

- If V1 shows that the training-channel share is materially below the
  paper's 99.4%, the paper must report the channels separately and must not
  call the template-channel rate "training contamination".
- If V2 finds that ≥50% of e421 or e422 accessions are `SAME_ACCESSION` or
  `SEQ95` before 2018-04-30, the paper must not call either census
  "leakage-free"; it may say "free of exact-structure temporal leakage".
- V3 results are reported whatever their direction. A between-stratum
  coverage difference is only described, never claimed, unless it exceeds
  3 combined MCSE, and even then it is labelled exploratory and post hoc.
- If the V3 overall pooled coverage does not reproduce the registered P5 mean
  within 3 MCSE, the paper cites the registered artifact and records the
  discrepancy as an open audit item.

## Addendum A1 (written 2026-10-07 AFTER the first V3 run — post hoc)

The first V3 run showed residue-weighted coverage of 0.905/0.951 for test
proteins with a pre-cutoff same-accession/≥95% relative versus 0.830/0.894
for proteins without, and that proteins without a close relative carry
more rows each (59,101 rows over 51 proteins vs 58,804 over 87). To separate
protein size from exposure, V3 now also reports, per stratum, the
protein-averaged coverage and the rows-per-protein distribution. These are
**not pre-specified**, are labelled as such in the output
(`by_stratum_protein_averaged_NOT_PRESPECIFIED`), and are reported whatever
their direction. One bug fix: a residue with non-finite lDDT made the e421
Spearman coefficient NaN; non-finite values are now excluded (n = 24,456,
matching the committed e421 analysis).

## Addendum A2 (2026-10-07, after A1 — post hoc sensitivity)

A miss-concentration check found that one accession, K7PQ54 (16,947 of
117,905 rows; 19 entries / 38 chains, all batch 4), carries 49.8% of all
test misses at level 0.90. Within the same entries, chains whose ledger rows
start at UniProt position 35–36 give median lDDT 0.99 (d ≈ 0.4 Å) while chains
starting at 70–71 give median lDDT 0.20 (d ≈ 20.7 Å) at median pLDDT 98.4;
PDBe SIFTS maps both chains of 9qj6 to UniProt 35–507. This pattern suggests a
residue-mapping artifact in the e422 label pipeline (unconfirmed). V3 therefore
also reports every e422 summary with K7PQ54 excluded, labelled
`sensitivity_excluding_K7PQ54` (post hoc; does not alter any registered
outcome). The anomaly is reported to the human for routing to the e422 owner.
Mechanism check (post hoc, same addendum): across the 200 splits, K7PQ54 fell
in train/calibration/test 99/44/57 times; mean 0.90-level half-width q was
0.688/0.826/0.708 and protein-averaged coverage of the *other* test proteins
0.891/0.972/0.882 respectively (script inline in the session; to be re-derived
by the checker).

## Addendum A3 (2026-10-08, written BEFORE running V5/V6/V2c)

Context: revision 3 of the paper (referee report recommended measuring the
template channel directly and reporting an exposure effect). All three
analyses below are new and are run only after this addendum is committed.

**V5 — template channel, measured (e422 submitted roster, 138 accessions).**
For each accession, the AFDB model file is the one the e422 census fetched,
read from `results/e422/raw` via the original batch-checkpoint receipts (no
re-download). From each file record the software/version string and every
template in `_ma_template_ref_db_details` with `db_name = PDB`. For each
template PDB ID, fetch the PDBe SIFTS UniProt mapping and the RCSB initial
release date (cached, hash-receipted, same client as V2). Per-accession
template class: SAME_ACC (≥ 1 template maps to the accession), OTHER (templates
present, none maps to it), NONE (no templates listed), UNRESOLVED (a template
lookup failed and no SAME_ACC was found; never folded into OTHER/NONE).
Report: class counts; cross-tab against the V2 training-cutoff class; the
latest template release date overall (sanity check against the documented
2021-02-15 template limit); software versions. Descriptive only; no
thresholds. Labels and pLDDT are not read.

**V2c — homology coverage extension.** Run the unchanged V2 query and class
rule for accessions in the corrected rosters that V2 did not cover (e422
corrected ledger: 3 new accessions; e421 corrected relabel: 22). Report them
separately from the original V2 table.

**V6 — exploratory exposure effect (e421 corrected relabel only).** Input:
`results/e427/e421_relabel/census_results.json` (status `ok`). Unit =
accession (rows pooled across that accession's ok entries). Per accession:
median lDDT; fraction of residues with lDDT < 0.60; fraction with pLDDT ≥ 90
and lDDT < 0.60 ("confident error"); median pLDDT. Groups by training-cutoff
class from V2/V2c: CLOSE (SAME_ACC or SEQ95) vs NOT_CLOSE (SEQ30 or NOVEL);
the four-class breakdown is also reported. Statistic: difference in group
means (CLOSE − NOT_CLOSE) of each per-accession quantity, with a 95%
percentile bootstrap over accessions (10,000 resamples, seed 0). Reported
whatever the direction, labelled EXPLORATORY, with group sizes; no
significance claim. **Not run on any e422 ledger**: pLDDT-conditional
quantities on e422 are registered predictions (P6/P7) that stay blind until
their registered batch-8 binding.

## Addendum A4 (2026-10-08, written BEFORE running V7/F-scope)

**V7 — trigger prevalence (label-free).** Over the PDBe SIFTS UniProt
mapping responses the e422 census fetched for the submitted roster (read from
`results/e422/raw` via checkpoint receipts), count mapped segments, and those
whose start has a null `author_residue_number` (the condition under which the
old fallback switched numbering schemes). For those segments, also report
the label-to-author offset at the first observed residue where determinable
from the reference mmCIF (`label_seq_id` vs `auth_seq_id`), as a
distribution. Report per segment and per accession. No labels, lDDT, or
pLDDT are read.

**F-scope — before/after figure (pLDDT-free).** For every chain segment in
the original (`results/e422/batches`) and corrected (`results/e427/batches`)
batch-4 cumulative ledgers: median lDDT before vs after (scatter; segments
present in only one ledger shown separately), and the residue-level lDDT
histogram before and after. pLDDT is deliberately not plotted: pLDDT-
conditional e422 quantities are registered predictions (P6/P7) that stay
blind until their batch-8 binding.

### Deviation D1 (2026-10-08, after V7 ran — post hoc)
V7 as registered reads no labels. Its output additionally joins the
triggered segments to the F-scope before/after median-lDDT changes
(`link_to_fscope_NOT_IN_A4_V7_TEXT`): all 147 triggered segments with a
nonzero label-to-author offset are among the 151 segments whose median lDDT
changed by > 0.2; the other 4 are zero-offset segments with an internal
numbering gap. This join was requested by the analyst, not pre-specified, and
is reported as post hoc. Extra descriptive counts (null end author numbers,
within-segment offset constancy) are also post hoc.
