# Independent audit reports

Write-ups and scripts from the independent checks cited in the paper. The
files in `census_correction/` and `paper_verification/` are reproduced
unchanged from the working repository (they are the record of what was run,
so they are not edited, not even to fix paths); internal session names that
appear in them identify which reviewer process ran each check.
`label_crosscheck_biotite/` was added for this release.

## census_correction/
Audit of the residue-numbering repair, the replay of the preregistered census,
the post-2022 relabel and the registered coverage rerun.

- `BRIEF.md` — what the auditor was asked to check (items 1–8).
- `AUDIT_E427.md` (+ `.ots`) — Phase A (items 1–7, 7 Oct 2026) and Phase B
  (item 8, the registered rerun, 9 Oct 2026). Verdict: PASS-with-notes.
- `independent_label_sample.*`, `independent_e421_audit.*`,
  `replay_rerun_report.json`, `independent_p5_audit.*` — the auditor's own
  scripts and outputs (written without importing the study's code).
- `PHASE_B_independent_subagent.md` (+ `.ots`) and `phaseB_subagent/` — a
  second, separate audit of the rerun (digest checks and an independent
  coverage recomputation). Verdict: PASS-with-notes.

## paper_verification/
Checks run on earlier drafts of the paper: `AUDIT_REVISION2.md` (the audit that
found the numbering error), recounts of the defect scope and of the exposure
tables, and the scripts behind them.

## label_crosscheck_biotite/
A recomputation of all 495 corrected chain segments with biotite 1.7.1 and
residue-level SIFTS from the PDBe updated mmCIF: 494 of 495 identical (every
lDDT within 5e-5); the remaining one, 21xg chain A, is an alignment tie next
to two crosslinked residues (18 values move by up to 0.054, median by 0.003,
share below 0.60 unchanged). Also counts the residues the labeller excludes
(57 modified residues other than MSE, 16 altloc-not-A). See its README.

## What each script checks and how to run it

"Layout" says where a script must sit, because the historical scripts compute
their root from their own location: `ROOT = HERE.parents[2]` means the script
must be at `<root>/results/e427/audit/<script>` inside a work root built by
`python -I code/tools/make_workroot.py _work` (which places
`independent_label_sample.py` and the two `phaseB_subagent` P5 scripts there;
copy the others the same way). "Compares automatically" says whether the
script itself tests its result against the registered or stored value; where
it does not, compare its printed output by hand (or use
`code/paper_numbers.py` / `code/reproduce_p5.py`, which do).

| Script | What it checks | Compares automatically | Layout | Needs (not in the release) |
|---|---|---|---|---|
| `census_correction/independent_label_sample.py` (= `code/audit/independent_label_sample.py`) | re-derives 36 corrected chain segments (26 targeted, 10 random unchanged) from raw bytes with its own parser and lDDT; checkpoint digests and replay provenance | yes: per-row values and checkpoint digests (`all_rows_match`, `checkpoint_integrity_pass`); its numeric tolerance (1e-4 on 4-decimal values) is looser than needed (5e-5), it samples single-segment chains only and trusts the stored self-check flags | `results/e427/audit/` | census raw files (`fetch_raw.py --set e422`) |
| `census_correction/independent_e421_audit.py` | re-derives the post-2022 relabel (status per entry, label rows) and recomputes `summary.json`, including the **pLDDT-binned lDDT < 0.60 error rates of the post-2022 (e421) census** — allowed, because that exploratory census is not blinded; it never touches the preregistered census | yes: against `results/e427/e421_relabel/summary.json` and `census_results.json` | `results/e427/audit/` | October 2026 source bytes (`fetch_raw.py --set e421relabel`) |
| `census_correction/independent_p5_audit.py` | P5 recomputation on the corrected ledger; digests of checkpoints, evaluation, receipt and binding sidecar | yes: against the P5 fields of the evaluation artifact | `results/e427/audit/` | the withheld batch-4 evaluation (`results/e427/evals/e422_batch04_evaluation.json`); cannot run from this release |
| `census_correction/phaseB_subagent/recompute_p5.py` (= `code/audit/e427_phaseB/recompute_p5.py`) | P5 recomputation, all 200 per-split coverages, k-convention diagnostic | no (prints mean, MCSE, threshold) | `python -I results/e427/audit/phaseB_subagent/recompute_p5.py $W` (writes next to itself) | nothing |
| `census_correction/phaseB_subagent/sensitivity_p5.py` (= `code/audit/e427_phaseB/sensitivity_p5.py`) | P5 on the original ledger; open-interval boundary sensitivity on the corrected ledger | no | `python -I …/sensitivity_p5.py $W` (writes to `$W/results/e427/audit/phaseB_subagent/`) | nothing |
| `census_correction/phaseB_subagent/verify_digests.py` | canonical digests of checkpoints, evaluation, receipt, binding sidecar; extracts allowed fields | yes (digest `match` flags) | `python -I verify_digests.py $W` | the withheld evaluation artifact; **see the warning below** |
| `paper_verification/p5_independent_reproduction.py` | P5 on the original ledger from the registration text, plus coverage by exposure stratum (not by pLDDT) | no: prints its values next to the registered evaluator's | explicit arguments (`--checkpoints`, `--v2`, `--registered-evaluation`, `--output`) | the withheld original batch-4 evaluation (`--registered-evaluation`) |
| `paper_verification/e422_signature_screen.py` | screens the original ledger for chains with median pLDDT ≥ 90 and median lDDT < 0.4 (the Table 3 screen; per-chain median pLDDT, not coverage) | no | explicit arguments (`--checkpoints`, `--output`) | nothing |
| `paper_verification/k7pq54_mapping_audit.py` (= `code/audit/k7pq54_mapping_audit.py`) | 9qj6 chains A/B mapped from SIFTS to coordinates by an independent parser; register shift; repaired lDDT | partly: reports whether its rows match the checkpoint rows | explicit arguments (`--pdb`, `--sifts`, `--af`, `--checkpoints`, `--output`) | 3 raw files (`fetch_raw.py --set 9qj6`) |
| `label_crosscheck_biotite/*.py` | see that directory's README | `full_ledger_biotite.py` prints the 494/495 count | argument-driven (run from the repository root) | census raw files, PDBe updated mmCIF (network), `biotite==1.7.1` |

### Interval boundary convention

The registered evaluator tests coverage as `c − q ≤ y ≤ c + q`, with the
bounds computed in floating point; the registration text ("ŷ ± q") does not
pin the form. Testing `|y − c| ≤ q` instead gives a different answer for
residues lying exactly on the boundary: about 500 residues of the **original**
ledger, raising the original-ledger means by 7e-5 (0.90) and 8e-5 (0.95).
This is why `AUDIT_REVISION2.md` (which used `|y − c| ≤ q`) reports
0.858036 / 0.915276 where the registered evaluator gives 0.857965 / 0.915196.
No residue of the **corrected** ledger is affected, so the corrected-ledger
values (0.894493 / 0.945578) are the same under either form.
`code/reproduce_p5.py` and `code/paper_numbers.py` use the registered form.

### Warning: blinding of the preregistered census

None of these files contains the pLDDT-binned results of the preregistered
census, which remain blind until its batch-8 evaluation. `verify_digests.py`
prints, besides the P5 fields, the binding *status strings* of the non-P5
predictions recorded in an evaluation artifact. Do not run it on an
evaluation or binding artifact from batch 8 or later before those
predictions have been evaluated and released. The post-2022 (e421) census
is exploratory and not blinded; its pLDDT-binned error rates (computed by
`independent_e421_audit.py`, `code/paper_numbers.py` and the e427 relabel)
are reported in the paper.
