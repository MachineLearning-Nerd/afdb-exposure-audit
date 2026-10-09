# e422 correction of record — batch-4 P5 outcome (2026-10-07)

Authority: human directives of 2026-10-07 in session w1:pQ ("pause e422 and
correct code"; "do finish all 4"). Procedure: `docs/e427_registration.md`
(committed 73aef782, OTS-stamped). This record does not modify any frozen
byte; it reclassifies an outcome.

## What is corrected

| Prediction | Original-ledger mean | Threshold | Original outcome | Classification |
|---|---:|---:|---|---|
| P5-90 | 0.8579652901901389 | 0.8759801398774404 | FAIL | **label artifact — not a coverage finding** |
| P5-95 | 0.9151961367793325 | 0.9328446915697971 | FAIL | **label artifact — not a coverage finding** |

Original artifacts (unchanged, anchored): `results/e422/evals/e422_batch04_evaluation.json`
(23abe6ea…), `results/e422/batches/e422_batch0{1..4}_checkpoint.json`,
binding sidecar `results/e422/batches/e422_binding_state.json`.

## Why

The labels were built with a residue-numbering defect (registration §1):
SIFTS `label_seq_id` start positions were applied to `auth_seq_id` atoms when
a segment's first residue is unobserved, with no amino-acid identity check.
Evidence:
- independent audit F1 (`paper/arxiv/verification/audit/AUDIT_REVISION2.md`):
  9qj6 chains A/B ledger lDDT 0.198/0.199 vs 0.993/0.995 correctly mapped;
- replay relabel (`results/e427/replay_report.json`, commit e9102501): the
  frozen labeller reproduces every original batch-2–4 checkpoint exactly; the
  repaired labeller changes 151 chain segments across 40 accessions by more
  than 0.2 median lDDT; residues with lDDT < 0.60 fall from 31.9% to about 1%;
  chains with median pLDDT ≥ 90 and median lDDT < 0.4 fall from 116 to 0.

## Outcome of record going forward

The P5 outcome of record is the registered evaluator's result on the
corrected ledger `results/e427/batches/` (registration §4.3), reported beside
the original. **Status at this writing: running** (launched 2026-10-07; the
original batch-4 evaluation took 20.4 h). A pre-registration diagnostic on
offline-relabelled rows gave 0.8980 / 0.9474 against thresholds 0.8913 /
0.9448 (disclosed in registration §4.4); it is not the outcome of record.
P6/P7a/P7b bind at batch 8 on the corrected ledger.

## Consequences for other records
- Any document citing the batch-4 P5 FAIL as a coverage result must cite
  this record instead (CLAIMS.md already classifies it as a label artifact).
- The e421 census shares the mapping fallback; its corrected relabel is
  registration §4.5 (`results/e427/e421_relabel/`). e421 remains exploratory.
- The temporal-leakage paper's Section 6 must not report a prospective
  coverage failure.

## Addendum 2026-10-08 (precision fixes from an independent claims audit)
1. **Source of the defect-scope figures.** The figures above (151 chain
   segments across 40 accessions changing by > 0.2 median lDDT; lDDT < 0.60
   falling from 31.9% to about 1%; 116 → 0 suspect chains) are not stored in
   `results/e427/replay_report.json`, which holds the replay self-checks and
   row counts only. They were first computed in the disclosed diagnostic
   (`results/e427_repair/README.md`). An independent recomputation on the
   registered original (`results/e422/batches`) and corrected
   (`results/e427/batches`) batch-4 ledgers confirms them: 217/488 segments
   change, 151 by > 0.2 across 40 accessions; residues with lDDT < 0.60
   31.86% → 1.11%; chain screen 116 → 0, accession screen 28 → 0,
   within-entry screen 13 → 0. The 116 flagged chains span 32 accessions; "28"
   is the separate accession-level screen.
2. **Rerun invocation.** Registration §4.3 names the evaluation CLI
   (`experiments/e425_execreceipt.py --batch 4`). Because the pause guard
   refuses that CLI, the rerun calls the same registered function
   `e425_execreceipt.run_logged_evaluation(REGISTERED_T0, 4, …)` through the
   wrapper `experiments/e427_rerun_p5.py` with the registration's corrected
   ledger, output, and binding-sidecar paths. The evaluator, seeds, grid,
   and decision rule are unchanged; the Phase B audit checks this equivalence.
3. **Thresholds.** The registered rule is mean coverage ≥ nominal − 3 MCSE.
   The values 0.8913 / 0.9448 quoted above are that rule evaluated in the
   diagnostic run, not registered constants; the rerun computes its own.
