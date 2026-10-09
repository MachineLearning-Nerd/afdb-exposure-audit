# Auditing AlphaFold Database Comparisons — code, registrations and derived data

This repository accompanies the paper

> Dinesh Jinjala. *Auditing AlphaFold Database Comparisons: Exposure Channels
> and a Residue-Numbering Error.* arXiv preprint, 2026 (identifier to be added).

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
third-party files are not redistributed; `data/fetch/` re-fetches them and
checks the digests.

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
python -I code/paper_numbers.py      # recomputes 127 numbers quoted in the paper -> "127/127 checks OK"
python -I code/reproduce_p5.py       # Table 4 from the ledgers -> "ALL MATCH" with the registered evaluator
```

`paper_numbers.py` prints, for each number, the value in the paper, the value
recomputed from `data/derived/`, and the source file. `reproduce_p5.py` is a
standalone implementation of the registered P5 protocol (it imports nothing
from `code/experiments/`) and reproduces the registered evaluator's means
exactly: 0.894492567684754 / 0.945577705711315 on the corrected labels
(thresholds 0.8914 / 0.9449, HOLD) and 0.857965290190139 / 0.915196136779333
on the original labels (thresholds 0.8760 / 0.9328, FAIL).

## Repository map

```
paper/                       temporal_leakage.tex, references.bib, figures/rev3/*.pdf, rendered PDF   (0.4 MB)
code/
  experiments/               registered pipeline code e420*, e421*, e422*–e427* (byte-identical;
                             e421_conformal.py is included only because the e422 freeze binds it
                             as the source of the registered pLDDT grid)
  tests/                     registered unit tests for e422–e427 (run in a work root, see below)
  verification/              V1–V7 scripts and make_figures_rev3.py (byte-identical)
  audit/                     independent-audit scripts whose outputs the paper cites
  reproduce_p5.py            standalone P5 recomputation (new)
  paper_numbers.py           recomputes every traceable paper number (new)
  tools/make_workroot.py     assembles the original directory layout in a scratch copy (new)
  tools/extract_p5_outcomes.py  how data/derived/p5_batch04_outcomes.json was produced (new)
registrations/               registrations, freeze/dispatch/binding/receipt records, correction
                             of record, runbook, verification registration, and .ots proofs   (0.4 MB)
  stamped_versions/          earlier exact versions of records whose OTS proof stamps that version
  verification/              V1–V7 registration (current text + the two stamped versions)
data/
  derived/                   derived tables, mirroring the original repository-relative paths (66 MB)
    results/e420/            ATLAS availability census, label-run receipt, overflow temporal check
    results/e421/            exploratory post-2022 census (original, pre-repair outputs)
    results/e422/            preregistered census: batch 1–4 checkpoints (original labels), freezes,
                             binding sidecar, manifest, execution receipt, .ots proofs
    results/e427/            corrected batch 1–4 checkpoints, replay report, post-2022 relabel,
                             execution receipt, audit outputs, .ots proofs
    results/e427_repair/     preliminary (unregistered) repair rows behind the 0.898/0.947 diagnostic
    paper/arxiv/verification/out/    V1–V7 outputs and figure numbers
    paper/arxiv/verification/audit/  9qj6 / K7PQ54 mapping audit (Figure 1 input)
    data/paper_verification/http/ and paper/arxiv/verification/cache/v2c/http/
                             RCSB Search/Data API responses used for homology exposure (CC0)
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
| Table 1; Sec. 5 (868, 791, 851, 77, 863, years 1988–2023, median 2008) | `python -I paper/arxiv/verification/v1_channel_recount.py $W rerun/v1.json` | `results/e420/label_run_receipt.json`, `overflow_temporal_check.json` | `out/v1_channel_recount.json` | offline |
| Sec. 4 ATLAS counts (1,938 / 1,735 / 943 / 869) | `experiments/e420_census.py`, `e420_label_run.py`, `e420_overflow_check.py` | ATLAS archive, PDBe, AFDB | `results/e420/source_availability_manifest.json` (`summary`), `label_run_receipt.json`, `overflow_set.json` | network; see limits below |
| Table 2 (homology), Sec. 6 (87/138, 128/138, P00698 730 …), App. B (224 accessions; 594 × 200, 408 × 204) | `python -I paper/arxiv/verification/v2_prior_exposure.py $W rerun/v2.json`, then `v2c_coverage_extension.py` | e421 census, e422 ledger + manifest, cached RCSB responses | `out/v2_prior_exposure.json`, `out/v2c_extension.json` | UniProt sequences fetched live (199 requests); all RCSB queries served from the shipped cache |
| Table 2 (templates), Sec. 6 (45/88/5; 43 vs 35; Q9F0J8; template dates) | `python -I paper/arxiv/verification/v5_templates.py $W rerun/v5.json` | AFDB model files in `results/e422/raw` (re-fetch, set `e422`), PDBe SIFTS, RCSB | `out/v5_templates.json` | raw files + network |
| Fig. 3 | `python -I paper/arxiv/verification/make_figures_rev3.py $W` | `out/v2…`, `out/v5…` | `paper/arxiv/figures/rev3/f3_exposure.pdf` | offline |
| Sec. 6 exposure effect (40 vs 45 proteins; 0.016 [0.003, 0.032] …) | `python -I paper/arxiv/verification/v6_exposure_effect.py $W rerun/v6.json` | `results/e427/e421_relabel/census_results.json`, `out/v2`, `out/v2c` | `out/v6_exposure_effect.json` | offline |
| Sec. 7.1, Fig. 1 (9qj6: shift 34; 0.198/0.199 → 0.993/0.995; 472/471) | `python -I paper/arxiv/verification/audit/k7pq54_mapping_audit.py --pdb … --sifts … --af … --checkpoints results/e422/batches/e422_batch04_checkpoint.json --output rerun/k7.json`; figure via `make_figures_rev3.py` | 3 raw files (`fetch_raw.py --set 9qj6`) | `audit/k7pq54_9qj6_audit.json`; `f1_mechanism.pdf` | 3 files from network; figure offline |
| Table 3, Fig. 2, Sec. 7.2 (488/495 segments, 188 changed, 151 > 0.2, 116 → 0 chains, 13 entries, +6,465/−916 residue positions) | `make_figures_rev3.py` (F2) and `code/paper_numbers.py` | `results/e422/batches`, `results/e427/batches` | `out/figures_rev3_numbers.json`, `f2_scope.pdf` | offline |
| The two ledgers themselves | original: `experiments/e423_dispatch.py` → `e423_runner.py` (batch 1: `e426_attested_empty.py`); replay + corrected: `experiments/e427_replay_relabel.py` (frozen labeller must reproduce every original checkpoint before the corrected one is written) | raw files (`--set e422`) | `results/e422/batches/*`, `results/e427/batches/*`, `results/e427/replay_report.json` | raw files; replay itself offline |
| Sec. 7.2 (K7PQ54: 14% of rows, half of 0.90-level misses) | `python -I paper/arxiv/verification/v4_k7pq54_diagnostics.py $W rerun/v4.json` | original ledger | `out/v4_k7pq54.json` | offline |
| Sec. 7.2 (independent audit, 36 segments) | `python -I results/e427/audit/independent_label_sample.py` | raw files + both ledgers | `results/e427/audit/independent_label_sample.json` | raw files |
| Sec. 7.3 (492 segments, 304 triggered, 147 offset, −31…+201) | `python -I paper/arxiv/verification/v7_trigger_prevalence.py $W rerun/v7.json` | raw SIFTS + mmCIF (`--set e422`), both ledgers | `out/v7_trigger_prevalence.json` | raw files |
| Sec. 7.4 post-2022 census (96 → 139 entries, 0.617 → 0.643, 0.58% → 0.29% …) | original: `experiments/e421_census_fixed.py` (**runs the whole census at import — never import it**); relabel: `experiments/e427_e421_relabel.py` | network (re-fetched 2026-10-07) | `results/e421/census_results.json`, `results/e427/e421_relabel/{census_results,summary,receipts}.json` | network |
| Table 4, Sec. 7.5 (registered evaluator) | `experiments/e427_rerun_p5.py` → `e425_execreceipt.run_logged_evaluation` → `e424_evalgate` → `e422_evaldriver` / `e422_protocol` (about 20 h) | `results/e427/batches` | `results/e427/evals/*` (withheld, see below) | offline |
| Table 4 (fast recomputation) | `python -I code/reproduce_p5.py` (from the repo root) | both ledgers | stdout | offline, seconds |
| Sec. 7.5 boundary sensitivity (0.9453 vs 0.9449) | `python -I code/reproduce_p5.py --open-interval`; `results/e427/audit/phaseB_subagent/sensitivity_p5.py $W` | ledgers | `sensitivity_p5.json` | offline |
| Sec. 7.5 independent reimplementation | `python -I results/e427/audit/phaseB_subagent/recompute_p5.py $W` | corrected ledger | `recompute_p5.json` | offline |
| Sec. 7.5 preliminary diagnostic (122,024 rows; 0.898 / 0.947) | rows: `experiments/e427_offline_relabel.py --out … --rows-out …`; statistic: `code/paper_numbers.py` | raw files | `results/e427_repair/repaired_rows_batches01-04.json.gz` | raw files for rows; statistic offline |
| Appendix A timeline | registrations and their `.ots` proofs | — | — | Bitcoin node or block explorer |
| Unit tests | `python -m pytest -q tests` | — | — | offline (`test_9qj6_real_bytes_repaired` needs `--set 9qj6`; deselect it otherwise) |

Re-running a script writes into `_work/`; compare with the shipped file, for
example `cmp rerun/v1.json paper/arxiv/verification/out/v1_channel_recount.json`.
Outputs that record a run time (`utc`) differ only in that field.

`code/verification/v3_stratified.py` (and `out/v3_stratified.json`) is the
pre-correction exploratory V3 analysis; its e422 part was computed on the
defective labels and is not used by the paper. It is kept for completeness.

### Verified for this release (2026-10-09)

Offline, from a fresh work root: `paper_numbers.py` 127/127; `reproduce_p5.py`
matches the registered evaluator to all printed digits on both ledgers; V1, V3
and V4 outputs byte-identical; V6 identical except its `utc` field;
`make_figures_rev3.py` reproduces `figures_rev3_numbers.json` byte-for-byte
and figure rasters identical to `paper/figures/rev3/`; the phase-B
`recompute_p5.py` reproduces all 400 per-split coverages and
`sensitivity_p5.py` its JSON byte-for-byte; 238 of 239 tests pass offline, and
the remaining one (`test_9qj6_real_bytes_repaired`) passes after
`fetch_raw.py --set 9qj6`. With network: the three
Figure 1 source files re-fetched with identical SHA-256, the 9qj6 mapping
audit reproduced its per-chain results, and V2 reproduced its 199-accession
table exactly with every RCSB query served from the shipped cache.

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
- **Raw third-party files** (about 405 MB of e422 source bytes, plus the e420
  and e421 source bytes): re-fetch with `data/fetch/fetch_raw.py`.
- **Not part of this paper**: the ATLAS label-run pool freeze
  (`results/e420/pool_freeze.json`, 36 MB, OpenTimestamps-stamped) and the rule
  files `docs/g1_successor_*` that `e420_label_run.py` reads belong to a
  separate analysis line; `e420_label_run.py` is included for provenance of
  its receipt but cannot be re-run from this release. These, the
  pre-repair e421 label ledger and internal audit write-ups are available on
  request.

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
  `code/tests`, `code/verification` and `code/audit` are byte-identical to the
  project repository (`CHANGES_FROM_ORIGINAL.md`). The six modules bound by
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
original batch-4 evaluation, receipt and binding sidecar (969934). Proofs made
on 7–9 October 2026 (e427 registration, correction of record, corrected
ledgers and evaluation, relabel outputs, verification registration) were still
pending calendar confirmation when this release was assembled; `ots upgrade`
completes them.

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
  the cached RCSB API responses): CC0 1.0. AlphaFold DB: CC BY 4.0.
  PDBe/SIFTS and UniProt: CC BY 4.0. ATLAS: CC BY-NC 4.0 (as recorded in
  `registrations/e420_registration.md`; see the ATLAS website,
  https://www.dsimb.inserm.fr/ATLAS); `data/derived/results/e420/` lists ATLAS
  entries and chains and is therefore distributed under CC BY-NC 4.0 terms for
  its ATLAS-derived content.

## Citation

See `CITATION.cff`. The arXiv identifier will be added on submission.
