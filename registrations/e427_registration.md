# e427 registration — e422/e421 label-mapping repair, correction of record, and corrected-ledger continuation

> **Status: REGISTERED 2026-10-07** (written and committed before the formal
> relabel, the corrected-ledger P5 rerun, the dispatch rewire, the pause
> lift, and the e421 relabel). Executed by session w1:pQ under the human's
> directives of 2026-10-07 ("pause e422 and correct code"; "do finish all
> 4"); main (Research1, w1:p1) informed. E-number e427 is taken by that
> directive; main may re-label it without changing content.
> **The frozen e422 registration text (git-blob 6e1b9ce1…) is not edited.**
> This document amends how e422 labels are constructed and which ledger the
> remaining e422 evaluations read; every e422 analysis rule, seed, grid,
> threshold, prediction, and schedule is carried unchanged.

## 1. Defect (confirmed)

The e422 label path — `parse_mappings` in the frozen `experiments/e422_runner.py`
and its live copy in `experiments/e423_runner.py` — builds each SIFTS residue
map from `author_residue_number or residue_number`. PDBe's `residue_number`
is in mmCIF `label_seq_id` space; the reference CA atoms are keyed by
`auth_seq_id`. Whenever a mapped segment starts on an unobserved residue the
author number is null and the map silently switches numbering spaces,
shifting every residue of the segment. The e422 label path has no per-residue
amino-acid identity check, so shifted pairs were labelled. Confirmed
independently (paper/arxiv/verification/audit/AUDIT_REVISION2.md F1): 9qj6
chains A/B, ledger lDDT 0.198/0.199 vs 0.993/0.995 with label-space mapping.
The e421 census (`experiments/e421_census_fixed.py`) uses the same fallback but
DOES drop amino-acid mismatches, so there the defect mainly excluded residues
and entries (e.g. `insufficient_aligned`) rather than mislabelling them.

## 2. Repair (code of record)

`experiments/e427_label_repair.py` (tests: `tests/test_e427_label_repair.py`):
1. residue maps keyed on `label_seq_id` and bound to `struct_asym_id`; no
   author-number fallback; a segment lacking either is the typed
   `E422_MAPPING_PARSE_FAILED`;
2. equal-span SIFTS segments map by offset; unequal-span segments (internal
   insertions/deletions) map by global alignment (match +2, mismatch −1,
   gap −2) of `_entity_poly_seq` to the UniProt sequence carried by the AF
   model, keeping identical aligned pairs only;
3. per-residue identity check: mismatched pairs are excluded and counted; a
   segment with > 10% mismatched pairs is `E427_SEQUENCE_IDENTITY_FAILED`.
All other label rules (Kabsch, lDDT-Cα, B-factor robust z, pLDDT from the AF
B column, MIN_ALIGNED = 30) are the e421-P §4 rules, unchanged.

## 3. Correction of record (e422 batch 4, P5)

The anchored batch-4 outcome P5-90 FAIL (0.8580) / P5-95 FAIL (0.9152) was
computed on labels produced by the defective path. It is retained byte-for-byte
as the **original-ledger outcome** and is classified as a **label artifact,
not a coverage finding**. The P5 outcome of record becomes the
corrected-ledger outcome of §4.3, reported beside the original.

## 4. Procedure

**4.1 Replay relabel of batches 1–4 (no network).** `experiments/e427_replay_relabel.py`
re-executes the registered `e423_runner.run_batch` for each batch with a
replay transport that serves every request from `results/e422/raw` by the
receipted URL and raw-byte SHA-256 (hash-verified), with the original
`completed_utc` as the clock and a fresh manifest chained batch to batch.
Self-checks (all must hold, else STOP and report): (a) with the frozen
labeller the replay reproduces each original checkpoint's `rows`,
`exclusions`, `manifest_checkpoint_sha256`, and `freeze` counts exactly;
(b) the replayed manifest after each batch equals the recorded
`manifest_checkpoint_sha256`. Then the replay with the e427 labeller writes
the corrected checkpoints to `results/e427/batches/e422_batch0k_checkpoint.json`
(same schema; `receipts` = the original fetch receipts; an added
`e427_relabel` block records the original checkpoint sha256, labeller, and
replay time; `sha256` recomputed by the registered canonical rule). Batch 1
(attested empty, e426) is carried byte-identical.

**4.2 Ledger layout from here on.** `results/e422/batches/` is frozen
(original ledger). The corrected cumulative ledger is `results/e427/batches/`;
batches 5–12 are dispatched there with the e427 labeller; the manifest
(`results/e422/manifest_checkpoint.json`) and raw cache (`results/e422/raw`)
continue unchanged. Evaluations read `results/e427/batches/` and write
`results/e427/evals/`, binding sidecar `results/e427/batches/e422_binding_state.json`.
The dispatch CLI refuses `--out-dir results/e422/batches` and the evaluation
CLI refuses `results/e422/*` ledger/result paths (typed
`E427_FROZEN_LEDGER_DIR`).

**4.3 Corrected-ledger P5 rerun.** The registered evaluation command
(`experiments/e425_execreceipt.py --batch 4`, unchanged evaluator, seeds
0–199, R = 200, registered thresholds) on `results/e427/batches/`. Its
P5-90/P5-95 outcomes are the P5 outcomes of record. P6/P7a/P7b remain
suppressed until batch 8 exactly as registered and bind on the corrected
ledger.

**4.4 Disclosure.** Before this registration, an unregistered diagnostic on
offline-relabelled rows (`results/e427_repair/`) gave pooled coverage 0.8980
(threshold 0.8913) and 0.9474 (threshold 0.9448). The rerun is therefore not
blind. No analysis choice is open: the evaluator, seeds, grid, thresholds,
and decision rule are the frozen e422 ones, and the repair is fixed by the
data format (label vs author numbering), not by outcomes. The 10% identity
threshold was set before any relabel ran; the observed maximum is 9
mismatches in a segment of hundreds of residues, so it binds nowhere.

**4.5 e421 corrected relabel.** The e421 census kept no raw bytes, so the 168
qualifying entries are re-fetched from the same public lanes (PDBe mappings,
PDBe entry files, AFDB model files) with bytes persisted under
`data/e427/e421_raw/` and per-request receipts, then labelled with the e427
labeller (chain-A and accession restrictions as in e421). Outputs:
`results/e427/e421_relabel/` (corrected census, per-entry old/new status, and
the e421 descriptive quantities recomputed: lDDT<0.60 rate by pLDDT bin,
Spearman, four-cell counts). e421 stays EXPLORATORY; files under
`results/e421/` are not modified. Source bytes may have been revised since
September; this is disclosed, not corrected for.

**4.6 Audit and pause lift.** LeakagePaperChecker (w1:pR) audits from a
fresh context: repair code and tests, replay self-checks, corrected
checkpoints, rewired CLIs, the §4.3 rerun, and §4.5. The pause
(`results/e422/PAUSE.json`) is lifted only after an audit verdict of PASS or
PASS-with-notes with no blocking item; then the dispatch schedule is re-armed
with the corrected `--out-dir`. Missed batch windows are dispatched late under
the registered missed-batch clause.

## 5. Stop conditions
Any self-check failure in 4.1, a test failure, an evaluator typed failure, or
an audit FAIL stops the procedure at that step; the pause stays in force and
the failure is recorded.
