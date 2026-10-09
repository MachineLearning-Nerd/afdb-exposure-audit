# Exposure and Residue Correspondence — code, registrations and derived data

This repository accompanies the paper

> Dinesh Jinjala. *Exposure and Residue Correspondence: Two Silent Failure
> Modes When Comparing AlphaFold Database Models with Experimental
> Structures.* arXiv preprint, 2026 (identifier to be added).

The paper measures two hidden failure modes of comparisons between AlphaFold
Database (AFDB) models and experimental structures: (i) exposure of the
reference, or a close relative, to AlphaFold2 through training, templates or
homology, measured on the ATLAS molecular-dynamics corpus and on a
preregistered census of PDB entries released after 11 September 2026; and
(ii) a residue-numbering error (PDBe/SIFTS `label_seq_id` positions applied to
atoms indexed by `auth_seq_id`) that manufactured confident "errors" and made
a preregistered conformal-coverage prediction fail. The original result is
retracted; on corrected labels the registered rerun meets the prediction.

Everything needed to trace each number in the paper is here: the analysis
code (byte-identical to the registered versions), the registrations with
their OpenTimestamps proofs, the derived tables, and the request receipts
whose SHA-256 digests identify the exact third-party bytes used. Raw
third-party structure files are not redistributed; `data/fetch/` re-fetches
them and checks the digests. The small API responses behind the exposure
tables (RCSB searches, UniProt sequences, PDBe/RCSB template lookups) are
shipped as request-keyed caches, so those analyses run offline.

Internal analysis identifiers used throughout (paper, Appendix A):
**e420** ATLAS exposure; **e421** exploratory post-2022 census;
**e422** preregistered census (dispatch/evaluation chain e422–e426);
**e427** residue-numbering repair, replay, relabel and corrected P5 rerun.
**V1–V7** are the paper-verification analyses registered in
`registrations/verification/REGISTRATION.md`.

## Quick check (offline, about one minute)

```sh
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python -I code/paper_numbers.py      # 220 checks of the paper's numbers -> "220/220 checks OK" / "ALL CHECKS OK"
python -I code/reproduce_p5.py       # Table 3 from the ledgers -> "ALL MATCH" with the registered evaluator
```

`paper_numbers.py` prints one line per check with its class, the value in
the paper, the value obtained, and the source:

- **RECOMPUTED (157)**: computed by the script itself from row- or
  request-level inputs in `data/derived/` — the census ledgers, the e420
  receipts (release dates against 2018-04-30 and the AFDB
  `modelCreatedDate`), the e421 census rows, the RCSB/UniProt caches (the V2/
  V2c homology classes are replayed query by query), the AFDB template table
  plus the PDBe/RCSB template-lookup cache (V5 classes), and the V6 group
  means and bootstrap intervals (10,000 protein resamples, seed 0,
  percentile), the post-2022 AFDB model table (model types), and the
  chance identity of residue pairs under a register shift (corrected ledger
  plus cached UniProt sequences), the census flow from enumerated entries to
  scored proteins (checkpoint counts), and the PDBe SIFTS chains of template
  6SWU. P5 means, MCSEs and thresholds are
  recomputed with the script's own implementation of the registered
  protocol. Percentages use computed numerators.
- **RECONCILED (43)**: a stored output compared with that recomputation at
  1e-12 (for example the V1/V2/V5/V6 outputs, the post-2022 `summary.json`,
  the registered evaluator's P5 mean, MCSE and threshold, the two P5
  reimplementations cited in Section 7.5, the per-segment row counts and
  medians of the biotite cross-check against the shipped ledger), or a
  stored summary compared with its own per-record fields (V7, the
  independent label audit, the replay report, the excluded-residue counts,
  the e420 manifest), plus the canonical digests of the 8 ledger
  checkpoints.
- **RE-READ (20)**: numbers that need raw third-party files not shipped
  here (the 9qj6 SIFTS segment and the raw-file realignment for Figure 1,
  the V7 trigger counts, the 36-segment independent label audit, the
  biotite cross-check, its count of SIFTS segments per scored chain segment
  and the 16 + 57 excluded residues in `audits/label_crosscheck_biotite/`,
  and the V5 count of templates confirmed by chain name); each is also covered by a RECONCILED
  check where the shipped data allow one.

Not checked by the script: design parameters stated as definitions (lDDT
radius and thresholds, the 0.60/0.4/90/30-residue/10% cut-offs, 2.5 Å, the
E-value, 200 splits and the 50/25/25 split, the SIFTS identity and coverage
cut-offs 0.90/0.80 of the ATLAS step, the 10,000 bootstrap resamples, the
alignment scores +2/−1/−2), the post-2022 census's candidate count (33,517
X-ray entries after 2022-06-02, RCSB search of 2026-09-05) and its
5,000-atom size filter, both stated in `registrations/e421_registration.md`
rather than in a data file, the 12
planned batches, the ATLAS archive and check dates (2026-08-30, 2026-09-04),
and the registration and amendment dates of the Appendix A timeline (9–10
September), which are attested by the `.ots` proofs rather than by a data
file. The paper's data statement describes
the 220 checks by class (157 / 43 / 20).

The script imports nothing from `code/experiments`, `code/verification` or
`reproduce_p5.py`. It exits 0 only if all 220 checks ran and passed, 1 if a
check differs, and 2 if an input is missing, a section fails, or the number of
checks is not 220. `reproduce_p5.py` is a standalone implementation of the
registered P5 protocol (it imports nothing from `code/experiments/`). Before
computing, it verifies each ledger checkpoint's stored canonical digest and
its `MANIFEST.sha256` entry; it then requires the mean, MCSE and threshold to
equal the registered evaluator's to 1e-12 and the outcome to be identical:
0.894492567684754 / 0.945577705711315 on the corrected labels (thresholds
0.8914 / 0.9449, HOLD) and 0.857965290190139 / 0.915196136779333 on the
original labels (thresholds 0.8760 / 0.9328, FAIL). Exit status 1 on a
mismatch, 2 on a missing or altered input.

**Blinding scope.** The preregistered census (e422 and its corrected ledger)
is blind for its pLDDT-conditional predictions until the batch-8
evaluation. Nothing here computes pLDDT-binned or pLDDT-conditional
coverage, Mondrian or CQR quantities for it. pLDDT is read from the census
ledgers only as the registered P5 eligibility filter and, for the
label-quality checks, as a per-chain median pLDDT (the Table 2 screen for
high-confidence chains with median lDDT below 0.4, and the 9qj6 chain): a
screening statistic, not a coverage quantity, and already reported in the
paper. The exploratory post-2022 (e421) census is not blinded; its
pLDDT-binned error rates are reported in the paper and recomputed here.

## Repository map

```
paper/                       temporal_leakage.tex (final), references.bib, rendered PDF,
                             figures/rev4/*.pdf (final figures) and figures/rev3/*.pdf (earlier revision)   (0.5 MB)
code/
  experiments/               registered pipeline code e420*, e421*, e422*–e427* (byte-identical;
                             e421_conformal.py is included only because the e422 freeze binds it
                             as the source of the registered pLDDT grid)
  tests/                     registered unit tests for e422–e427 (run in a work root, see below)
  verification/              V1–V7 scripts, make_figures_rev3.py and make_figures_rev4.py (byte-identical)
  audit/                     independent-audit scripts whose outputs the paper cites
  reproduce_p5.py            standalone P5 recomputation with input digest checks (new)
  paper_numbers.py           checks every traceable paper number (new; RECOMPUTED/RECONCILED/RE-READ)
  tools/make_workroot.py     assembles the original directory layout in a scratch copy (new)
  tools/extract_p5_outcomes.py  how data/derived/p5_batch04_outcomes.json was produced (new)
  tools/extract_v5_template_table.py  how the AFDB template table was produced (new)
  tools/extract_e421_model_table.py   how the post-2022 AFDB model table was produced (new)
audits/                      independent audit write-ups and scripts (historical, unchanged), and
                             label_crosscheck_biotite/ (new: biotite cross-check of the corrected labels)
registrations/               registrations, freeze/dispatch/binding/receipt records, correction
                             of record, runbook, verification registration, and .ots proofs   (0.4 MB)
  stamped_versions/          earlier exact versions of records whose OTS proof stamps that version
  verification/              V1–V7 registration (current text + the two stamped versions)
data/
  derived/                   derived tables, mirroring the original repository-relative paths (77 MB)
    results/e420/            ATLAS availability census, label-run receipt, overflow temporal check
    results/e421/            exploratory post-2022 census (original, pre-repair outputs)
    results/e422/            preregistered census: batch 1–4 checkpoints (original labels), freezes,
                             binding sidecar, manifest, execution receipt, .ots proofs
    results/e427/            corrected batch 1–4 checkpoints, replay report, post-2022 relabel
                             (with afdb_model_table.json: model type of its 93 AFDB files),
                             execution receipt, audit outputs, .ots proofs
    results/e427_repair/     preliminary (unregistered) repair rows behind the 0.898/0.947 diagnostic
    paper/arxiv/verification/out/    V1–V7 outputs and figure numbers
    paper/arxiv/verification/audit/  9qj6 / K7PQ54 mapping audit (Figure 1 input)
    data/paper_verification/http/ and paper/arxiv/verification/cache/v2c/http/
                             RCSB Search/Data API responses (CC0) and UniProt FASTA responses
                             (CC BY 4.0) used for homology exposure (V2, V2c)
    paper/arxiv/verification/cache/v5/http/
                             PDBe SIFTS mapping (CC BY 4.0) and RCSB entry (CC0) responses for
                             the 413 AFDB template PDB IDs (V5)
    paper/arxiv/verification/v5_afdb_template_table.json
                             templates, software and target of the 138 AFDB model files (V5 input)
    p5_batch04_outcomes.json key-selected P5 outcome blocks of the two batch-4 evaluations
  fetch/                     fetch_raw.py + README: re-fetch raw inputs and verify SHA-256
```

`data/derived/` deliberately mirrors the paths the registered code was frozen
with (`results/e427/batches/...`, `paper/arxiv/verification/out/...`), so the
code runs unchanged. Every copied file is byte-identical to the project
repository; see `CHANGES_FROM_ORIGINAL.md` for the SHA-256 of each file and
its origin, and `MANIFEST.sha256` for the whole tree.

## Environment

CPython 3.12 with the pinned packages in `requirements.txt` (numpy 1.26.4,
scipy 1.17.1, scikit-learn 1.9.0, matplotlib 3.11.0, requests 2.34.2,
pytest 9.1.1, opentimestamps-client 0.7.2). The figures use Liberation Sans if
installed (DejaVu Sans otherwise, which changes glyph metrics but not
content). Run scripts with `python -I`.

## Running the registered code: the work root

Several registered scripts compute paths relative to their own location or
expect `docs/`, `results/` and `paper/arxiv/verification/` side by side. Rather
than edit them, build a scratch copy of the original layout:

```sh
python -I code/tools/make_workroot.py _work      # copies, never links; _work/ is git-ignored
cd _work
W=$PWD                                           # several scripts need an absolute root
```

| Paper item | Script (in `_work/`) | Inputs | Output | Needs |
|---|---|---|---|---|
| Sec. 5 (868, 791, 851, 77, 863, years 1988–2023, median 2008) | `python -I paper/arxiv/verification/v1_channel_recount.py $W rerun/v1.json` | `results/e420/label_run_receipt.json`, `overflow_temporal_check.json` | `out/v1_channel_recount.json` | offline |
| Sec. 4 ATLAS counts (1,938 / 1,735 / 943 / 869) | `experiments/e420_census.py`, `e420_label_run.py`, `e420_overflow_check.py` | ATLAS archive, PDBe, AFDB | `results/e420/source_availability_manifest.json` (`summary`), `label_run_receipt.json`, `overflow_set.json` | network; see limits below |
| Table 1 (homology), Sec. 6.1 (87/138, 128/138, P00698 730 …), Sec. 4 (224 proteins queried; 594 queries with hits, 408 without) | `python -I paper/arxiv/verification/v2_prior_exposure.py $W rerun/v2.json`, then `python -I paper/arxiv/verification/v2c_coverage_extension.py $W rerun/v2c.json` | e421 census, e422 ledger + manifest, cached RCSB and UniProt responses | `out/v2_prior_exposure.json`, `out/v2c_extension.json` | offline (every request served from the shipped caches) |
| Table 1 (templates), Sec. 6.2 (45/88/5; 43 vs 35; Q9F0J8; template dates) | `python -I paper/arxiv/verification/v5_templates.py $W rerun/v5.json` | the 138 AFDB model files in `results/e422/raw` (`fetch_raw.py --set e422`, 38.6 MB), cached PDBe SIFTS and RCSB responses | `out/v5_templates.json` | AFDB model files; lookups offline. Without the model files, `paper_numbers.py` recomputes the classes from `v5_afdb_template_table.json` and the same cache |
| Fig. 3 of the earlier revision (not in the final paper; its numbers are in Table 1) | `python -I paper/arxiv/verification/make_figures_rev3.py $W` | `out/v2…`, `out/v5…` | `paper/arxiv/figures/rev3/f3_exposure.pdf` | offline |
| Sec. 4 (Model versions): post-2022 model types (93 files: 6 ColabFold v1.5.2, 87 AlphaFold Monomer v2.0; no ColabFold model among the 85 proteins of the exposure comparison) | `python -I code/tools/extract_e421_model_table.py data/derived _work/data/e427/e421_raw data/derived/results/e427/e421_relabel/afdb_model_table.json` (from the repo root); counts: `code/paper_numbers.py` | the AFDB files of the post-2022 relabel (`fetch_raw.py --set e421relabel`) | `afdb_model_table.json` | AFDB files for the table; counts offline |
| Sec. 6.3 exposure effect (40 vs 45 proteins; 0.016 [0.003, 0.032] …) | `python -I paper/arxiv/verification/v6_exposure_effect.py $W rerun/v6.json` | `results/e427/e421_relabel/census_results.json`, `out/v2`, `out/v2c` | `out/v6_exposure_effect.json` | offline. V6 reads `out/v2…` and `out/v2c…` from the `out/` directory next to the script, not from the root argument: run the work-root copy so that both point to the same tree |
| Sec. 7.1, Fig. 1 (9qj6: shift 34; 0.198/0.199 → 0.993/0.995; 472/471) | `python -I paper/arxiv/verification/audit/k7pq54_mapping_audit.py --pdb … --sifts … --af … --checkpoints results/e422/batches/e422_batch04_checkpoint.json --output rerun/k7.json`; figure via `make_figures_rev4.py $W rerun/fig` | 3 raw files (`fetch_raw.py --set 9qj6`) | `audit/k7pq54_9qj6_audit.json`; `f1_mechanism.pdf` | 3 files from network; figure offline |
| Table 2, Fig. 2, Sec. 7.2 (488/495 segments, 188 changed, 151 > 0.2, 116 → 0 chains, 13 entries, +6,465/−916 residue positions) | `python -I paper/arxiv/verification/make_figures_rev4.py $W rerun/fig` (F2) and `code/paper_numbers.py` | `results/e422/batches`, `results/e427/batches` | `out/figures_rev3_numbers.json`, `f2_scope.pdf` | offline |
| The two ledgers themselves | original: `experiments/e423_dispatch.py` → `e423_runner.py` (batch 1: `e426_attested_empty.py`); replay + corrected: `experiments/e427_replay_relabel.py` (frozen labeller must reproduce every original checkpoint before the corrected one is written) | raw files (`--set e422`) | `results/e422/batches/*`, `results/e427/batches/*`, `results/e427/replay_report.json` | raw files; replay itself offline |
| Sec. 7 opening (K7PQ54: 14% of rows, half of 0.90-level misses) | `python -I paper/arxiv/verification/v4_k7pq54_diagnostics.py $W rerun/v4.json` | original ledger | `out/v4_k7pq54.json` | offline |
| Sec. 7.2 (independent reimplementation, 36 segments) | `python -I results/e427/audit/independent_label_sample.py` (the release copy is `code/audit/independent_label_sample.py` = `audits/census_correction/independent_label_sample.py`; `make_workroot.py` places it here because it locates the root as `parents[2]` of its own path) | raw files + both ledgers | `results/e427/audit/independent_label_sample.json` (shipped: `data/derived/results/e427/audit/`) | raw files |
| Sec. 7.2 (second label check: 494 of 495 segments; 21xg A; 16 altloc and 57 modified residues excluded) | `audits/label_crosscheck_biotite/*.py` (see its README) | raw files + PDBe updated mmCIF | `audits/label_crosscheck_biotite/*_out.json` | raw files, network, `biotite==1.7.1` |
| Sec. 7.3 (492 segments, 304 triggered, 147 offset, −31…+201) | `python -I paper/arxiv/verification/v7_trigger_prevalence.py $W rerun/v7.json` | raw SIFTS + mmCIF (`--set e422`), both ledgers | `out/v7_trigger_prevalence.json` | raw files |
| Sec. 7.4 post-2022 census (96 → 139 entries, 0.617 → 0.643, 0.58% → 0.29% …) | original: `experiments/e421_census_fixed.py` (**runs the whole census at import — never import it**); relabel: `experiments/e427_e421_relabel.py` | network (re-fetched 2026-10-07) | `results/e421/census_results.json`, `results/e427/e421_relabel/{census_results,summary,receipts}.json` | network |
| Table 3, Sec. 7.5 (registered evaluator) | `experiments/e427_rerun_p5.py` → `e425_execreceipt.run_logged_evaluation` → `e424_evalgate` → `e422_evaldriver` / `e422_protocol` (about 20 h) | `results/e427/batches` | `results/e427/evals/*` (withheld, see below) | offline |
| Table 3 (fast recomputation) | `python -I code/reproduce_p5.py` (from the repo root) | both ledgers | stdout | offline, seconds |
| Sec. 7.5 boundary sensitivity (0.9453 vs 0.9449) | `python -I code/reproduce_p5.py --open-interval` (from the repo root); in `_work/`: `python -I results/e427/audit/phaseB_subagent/sensitivity_p5.py $W` (release copy: `code/audit/e427_phaseB/sensitivity_p5.py` = `audits/census_correction/phaseB_subagent/sensitivity_p5.py`) | ledgers | `results/e427/audit/phaseB_subagent/sensitivity_p5.json` | offline |
| Sec. 7.5 independent reimplementation | in `_work/`: `python -I results/e427/audit/phaseB_subagent/recompute_p5.py $W` (release copy: `code/audit/e427_phaseB/recompute_p5.py`) | corrected ledger | `results/e427/audit/phaseB_subagent/recompute_p5.json` | offline; prints its values, does not compare them (see `audits/README.md`) |
| Sec. 7.5 preliminary diagnostic (122,024 rows; 0.898 / 0.947) | rows: `experiments/e427_offline_relabel.py --out … --rows-out …`; statistic: `code/paper_numbers.py` | raw files | `results/e427_repair/repaired_rows_batches01-04.json.gz` | raw files for rows; statistic offline |
| Appendix A timeline | registrations and their `.ots` proofs | — | — | Bitcoin node or block explorer |
| Unit tests | `python -m pytest -q tests` | — | — | offline (`test_9qj6_real_bytes_repaired` needs `--set 9qj6`; deselect it otherwise) |

Re-running a script writes into `_work/`; compare with the shipped file, for
example `cmp rerun/v1.json paper/arxiv/verification/out/v1_channel_recount.json`.
Outputs that record a run time (`utc`) differ only in that field.

`code/verification/v3_stratified.py` (and `out/v3_stratified.json`) is the
pre-correction exploratory V3 analysis; its e422 part was computed on the
defective labels and is not used by the paper. It is kept for completeness.
Unlike the other V scripts it takes three arguments:
`python -I paper/arxiv/verification/v3_stratified.py $W paper/arxiv/verification/out/v2_prior_exposure.json rerun/v3.json`.

`v1_channel_recount.py` enumerates the admitted ATLAS pool through the typed
failures of the label-run receipt; that is complete only because no admitted
row reached the label ledger (`summary.ledger_rows == 0`). The registered
script is kept byte-identical and does not test this;
`code/paper_numbers.py` does (check "receipt: ledger_rows == 0").

### Verified for this release (2026-10-09)

Offline (network namespace disabled), from a copy of the release tree:
`paper_numbers.py` 208/208 (147 RECOMPUTED, 43 RECONCILED, 18 RE-READ) at revision 4; revision 6 added 12 checks (section `rev6()`: census flow, ATLAS network errors, sample composition, confident-error counts, template chain names, the post-2022 protein no longer scored), and 220/220 pass;
`reproduce_p5.py` verifies the 8 checkpoint digests and manifest entries and
matches the registered evaluator's mean, MCSE and threshold to 1e-12 on both
ledgers. 23 injected faults (including a receipt release date moved across
2018-04-30, a deleted V6 field, a post-2022 entry status changed without
updating `summary.json`, a per-protein V5 class or V7 offset changed without
updating the summary counts, 791 → 790 in `out/v1`, P5 MCSE + 1e-10, a
ledger value + 1e-4, a deleted or altered cache record), and 16 more on the
checks added for the final paper (a ColabFold model relabelled or attributed
to a labelled protein, a model-file digest, the excluded-residue total or one
per-entry count, a biotite row count or deviation, the 21xg detail, a replay
digest, the two P5 reimplementations and the open-interval output at
1e-9/1e-10, a UniProt sequence behind the chance-identity figure, the 9zxa
exclusions, a receipt date, one post-2022 confident error) each make
`paper_numbers.py` exit non-zero. V2 and V2c rerun offline from the shipped
caches and reproduce `out/v2_prior_exposure.json` and `out/v2c_extension.json`
exactly apart from `utc`, receipt order and the receipts' `cached` flag. V5,
given the 138 AFDB model files, reruns offline and reproduces
`out/v5_templates.json` apart from the same fields and the status of the 11
lookups that had returned HTTP 404 (not cached by design; offline they fail
with status −1, giving the same classes). V1, V3
and V4 outputs byte-identical; V6 identical except its `utc` field;
`make_figures_rev3.py` reproduces `figures_rev3_numbers.json` byte-for-byte
and figure rasters identical to `paper/figures/rev3/`, and `make_figures_rev4.py`
(final figures) reproduces the same JSON byte-for-byte and the two PDFs in
`paper/figures/rev4/` pixel-identically (only `/CreationDate` differs); the
paper PDF builds with tectonic without undefined references or overfull boxes; the phase-B
`recompute_p5.py` reproduces all 400 per-split coverages and
`sensitivity_p5.py` its JSON byte-for-byte; 238 of 239 tests pass offline, and
the remaining one (`test_9qj6_real_bytes_repaired`) passes after
`fetch_raw.py --set 9qj6`. With network: the three
Figure 1 source files re-fetched with identical SHA-256 and the 9qj6 mapping
audit reproduced its per-chain results. With the census source bytes and the
PDBe updated mmCIF files, the scripts in `audits/label_crosscheck_biotite/`
reproduced the reviewer's cross-check output exactly (494 of 495 segments).

## What is withheld, and why

- **Full batch-4 evaluation artifacts** (`results/e422/evals/` and
  `results/e427/evals/e422_batch04_evaluation.json`). Besides P5 they contain
  sections (pLDDT bins, six-bin floor table, decided grid) that feed the
  registered predictions P6/P7, which stay blind until the census's batch-8
  evaluation. `data/derived/p5_batch04_outcomes.json` holds only the P5-90/P5-95
  outcome blocks, `bound_at_batch`, the stored canonical digest and the file
  SHA-256, extracted by key selection (`code/tools/extract_p5_outcomes.py`).
  The `.ots` proofs of both full artifacts are included; their stamped digests
  equal the `file_sha256` values in that extract. The execution receipts and
  binding-state sidecars contain no binned quantities and are included. The
  full artifacts will be released after the census closes (batch 12, late
  November 2026). The authors have not computed pLDDT-conditional coverage on
  the released ledgers.
- **Raw third-party files** (about 405 MB of e422 source bytes, including the
  38.6 MB of AFDB model files V5 reads, plus the e420 and e421 source bytes):
  re-fetch with `data/fetch/fetch_raw.py`. The template information V5 takes
  from the AFDB files is shipped as `v5_afdb_template_table.json`.
- **Not part of this paper**: the ATLAS label-run pool freeze
  (`results/e420/pool_freeze.json`, 36 MB, OpenTimestamps-stamped) and the rule
  files `docs/g1_successor_*` that `e420_label_run.py` reads belong to a
  separate analysis line; `e420_label_run.py` is included for provenance of
  its receipt but cannot be re-run from this release. These, the
  pre-repair e421 label ledger and internal audit write-ups are available on
  request.
- **Not reproducible from this release**: the step that reduced the
  post-2022 census's candidates (X-ray entries at 2.5 Å or better released
  after 2022-06-02; 33,517 on 2026-09-05) to the 168 entries of
  `data/derived/results/e421/qualifying_entries.json`. It included a filter on
  entry size (fewer than 5,000 atoms; `registrations/e421_registration.md`),
  and its code is not in this release; only the resulting list is. That census
  also scored author chain A of each entry only (`CHAIN = "A"` in
  `code/experiments/e427_e421_relabel.py`; `chain = "A"` in
  `e421_census_fixed.py`).

## Provenance notes

- **Two kinds of replay.** The preregistered census was replayed from its
  original stored, hash-checked source bytes: with the original labeller the
  replay reproduced every stored batch exactly
  (`results/e427/replay_report.json`), and the corrected labeller was then run
  on the same bytes. The post-2022 relabel (Section 7.4) used source files
  **fetched again in October 2026**, because the September files of the
  exploratory e421 census had not been kept; it is exploratory and its
  receipts identify the October bytes.
- **Registered code is unmodified.** All files under `code/experiments`,
  `code/tests`, `code/verification` and `code/audit`, and the historical
  scripts under `audits/census_correction` and `audits/paper_verification`,
  are byte-identical to the project repository (`CHANGES_FROM_ORIGINAL.md`);
  improved checks live in new files (`code/paper_numbers.py`,
  `code/reproduce_p5.py`, `audits/label_crosscheck_biotite/`). The six modules bound by
  the e422 gate-(a) freeze (`e422_protocol`, `e422_manifest`, `e422_status`,
  `e422_runner`, `e422_evaldriver`, `e421_conformal`) match the SHA-256
  digests recorded in `registrations/e422_freeze_gate_a_20260909.md`. Comments in some
  registered files and registration texts refer to internal project records
  (a decisions log, status files, cron prompts) and to working-session
  identifiers of the form `w1:pX`; those records are not part of this release
  and the texts are left byte-identical so that their digests and timestamps
  remain verifiable.
- **Derived files keep their original bytes**, including absolute paths of
  the machine that produced them in two files
  (`results/e420/label_run_receipt.json` `summary`, and the `inputs` /
  `checkpoint` fields of `paper/arxiv/verification/audit/k7pq54_9qj6_audit.json`),
  because their SHA-256 digests are recorded elsewhere.
- **Not blind.** As the paper states, a preliminary, unregistered P5
  diagnostic (0.898 / 0.947) was known before the registered rerun; it is
  reproducible from `results/e427_repair/`.
- `e421_census_fixed.py` executes its census (network requests) at module
  import. Run it only deliberately, never import it.
- `e422_dispatch.py` is a refuse-to-run tombstone; registered dispatch runs
  through `e423_dispatch.py` (see `registrations/e423_dispatch_repair_record_20260910.md`).

## OpenTimestamps proofs

Each `X.ots` proves that a file with the SHA-256 of `X` existed when the
proof's Bitcoin transaction was mined.

```sh
pip install opentimestamps-client
ots info  registrations/e422_registration.md.ots          # shows the stamped SHA-256 and attestations
ots verify registrations/e422_registration.md.ots         # needs a local Bitcoin Core node
ots upgrade some_file.ots                                 # completes proofs still pending at the calendars
```

Without a local node, `ots verify` checks that the file hash matches and then
reports that it cannot reach a Bitcoin node; the attested block height printed
by `ots info` can be checked against any block explorer, or the file and proof
can be dropped into <https://opentimestamps.org>. Proofs already anchored in
Bitcoin at release: e422 registration and freeze (block 966254), e423 and e424
records (966262, 966273), e425 record (966287), e421 post-run checkpoint
(965898), e422 batch 1–4 checkpoints (967512, 967779, 968772, 969800), the
original batch-4 evaluation, receipt and binding sidecar (969934). The proofs
made on 7–9 October 2026, and some older ones, were completed with
`ots upgrade` on 2026-10-09; the upgraded proofs are shipped (their own
SHA-256 changed, the stamped digests did not; see `CHANGES_FROM_ORIGINAL.md`):
e427 registration, corrected batch 2–4 checkpoints and replay report (first
attestation in block 970355), post-2022 relabel outputs (970359), verification
registration A3/A4 and the 8 October addendum of the correction of record
(970417), corrected batch-4 evaluation, receipt and binding sidecar, and the
9 October addendum (970581), correction of record and `AUDIT_E427.md`
(970582). Every shipped `.ots` now carries at least one Bitcoin attestation.

Where a record was amended after it was stamped, the current text is shipped
together with the exact stamped version in `registrations/stamped_versions/`
(or `registrations/verification/REGISTRATION_A3_*.md` and `_A4_*.md`, and
`data/derived/results/e422/manifest_checkpoint.stamped_440566b4.json`), and the
proof sits next to the version it stamps. `registrations/verification/REGISTRATION.md`
(current) adds deviation D1, written after V7 ran, to the stamped A4 version.

## Web sources cited in the paper

Archived 2026-10-08 (not redistributed; re-fetch and compare):

| Source | URL | SHA-256 |
|---|---|---|
| AFDB FAQ text (templates released before 2021-02-15) | https://alphafold.ebi.ac.uk/chunk-7VK47MF2.js | `63fbdfbee9b98701104042203dcb31afa8ee76b24c5c0d6d5f0b37b578c14684` |
| PDBe AFDB release notes (v6 metadata; coordinates from v4) | https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes | `0ba8ab8c500a15b055e735b3bf9925cc4468736ccb5b6c3c8b72bf9f9a0dda63` |

## Licences

- **Code** (`code/`, `data/fetch/*.py`): MIT (`LICENSE`).
- **Paper text, figures and derived data** (`paper/`, `registrations/`,
  `data/derived/`): CC BY 4.0 (`LICENSE-CC-BY-4.0`), except where third-party
  terms below apply.
- **Third-party data keep their own terms.** PDB data (RCSB/wwPDB, including
  the cached RCSB API responses): CC0 1.0. AlphaFold DB (including the
  template table extracted from AFDB model files): CC BY 4.0. PDBe/SIFTS
  (including the cached SIFTS mapping responses) and UniProt (including the
  cached UniProt FASTA responses; The UniProt Consortium): CC BY 4.0. ATLAS: CC BY-NC 4.0 (as recorded in
  `registrations/e420_registration.md`; see the ATLAS website,
  https://www.dsimb.inserm.fr/ATLAS); `data/derived/results/e420/` lists ATLAS
  entries and chains and is therefore distributed under CC BY-NC 4.0 terms for
  its ATLAS-derived content.

## Citation

See `CITATION.cff`. The arXiv identifier will be added on submission.
