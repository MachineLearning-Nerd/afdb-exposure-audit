# e427 Phase B (brief item 8): independent audit of the registered corrected-ledger P5 rerun

**Who ran this:** a subagent of the `main` session ran this audit from a fresh context. It is not the designated checker. The designated checker pane (LeakagePaperChecker, codex) could not start its command host, so this subagent did the work instead. The subagent started from the recorded artifacts and did not rely on the builder's narrative or the commit message. Read-only outside `results/e427/audit/phaseB_subagent/` and this file. No commit. `AUDIT_E427.md` and `BRIEF.md` were not edited.

**Audit date:** 2026-10-09. **HEAD:** f79778d2 (branch paper-rigor).

**Blinding:** the scripts extract only these from the evaluation artifact, receipt, and sidecar:
- digests;
- ledger and provenance fields;
- P5 fields and pooled marginal aggregates;
- P6/P7 binding *status strings*.

The scripts did not print or compute any per-bin, floor-table, Mondrian, CQR, joint-grid, decided-grid, or P6/P7 quantities. To verify the evaluation digest, the full object is hashed but never displayed.

**Not blind:** this rerun is not blind. Registration §4.4 discloses an earlier unregistered diagnostic on offline-relabelled rows: 0.8980 against threshold 0.8913, and 0.9474 against 0.9448. Report the result as a registered, deterministic rerun with that disclosure, not as a blind replication.

## (a) Digests, ledger references, row count — PASS

Canonical rule, re-implemented in `phaseB_subagent/verify_digests.py` without importing project code: `json.dumps(obj minus "sha256", sort_keys=True, separators=(",",":"))`, UTF-8, then SHA-256. This rule was read from `experiments/e422_manifest.py::canonical_json_bytes`.

| Artifact | Recorded `sha256` | Recomputed | Match |
|---|---|---|---|
| `results/e427/evals/e422_batch04_evaluation.json` | f579c541…fd2c0c | f579c541…fd2c0c | yes |
| `results/e427/evals/e422_batch04_execution_receipt.json` | 1cf9a3e1…f3a6fd | 1cf9a3e1…f3a6fd | yes |
| `results/e427/batches/e422_binding_state.json` | c3ee0f6d…bc41f | c3ee0f6d…bc41f | yes |
| corrected checkpoint b01 | 4f54057d… | 4f54057d… | yes (0 rows) |
| corrected checkpoint b02 | a18a2274… | a18a2274… | yes (31,868 rows) |
| corrected checkpoint b03 | d0677b47… | d0677b47… | yes (45,903 rows) |
| corrected checkpoint b04 | ad6166b5… | ad6166b5… | yes (45,683 rows) |

- **Ledger references.** The evaluation's `ledger` block lists the digests a18a2274, d0677b47, and ad6166b5 for batches 2–4, with n_rows 31,868, 45,903, and 45,683. These are the **corrected** checkpoints in `results/e427/batches`. The original frozen checkpoints in `results/e422/batches` are 312d9ab0, c356084f, and d860ab55, with 29,054, 45,405, and 43,446 rows. Batch 1 is identical in both ledgers (4f54057d, empty).
- **Row counts.** `n_rows_cumulative` = inner `n_rows_input` = **123,454**, which equals the sum of the corrected checkpoint rows. The inner evaluation reports 141 accessions and an all-zero exclusion tally. The SHA-256 of its canonical accession list (7234a502…) equals the one my own loader derives from the corrected checkpoints.
- **Other evaluation fields:**
  - registration citation: blob 6e1b9ce1, content dc7a3627, commit 9dd2380b. `git hash-object docs/e422_registration.md` = 6e1b9ce1, so the frozen text is unchanged;
  - seeds [0, 199]; levels [0.9, 0.95]; scope `BATCH4_P5_ONLY_REGISTERED_SUPPRESSION`; status `EVALUATED`;
  - `floor_decision` = null and `evaluation_on_decided_grid` = null, as expected at batch 4;
  - `prior_binding_path` = `results/e427/batches/e422_binding_state.json`.
- **Receipt.** Started 2026-10-08T07:07:27Z and completed 2026-10-09T03:33:17Z. Wall time is 73,550 s (20.43 h), within the 24 h registered budget, and `deviation_note` is null. `evaluated_utc` matches the artifact.
- **File identity.** The working-tree evaluation file is byte-identical to the blob committed in f79778d2 (file sha256 a003c5c3…).
- **OTS.** `.venv/bin/ots verify` on all three `.ots` sidecars found no hash mismatch. All calendars report "Pending confirmation in Bitcoin blockchain", so the stamps are not yet Bitcoin-confirmed.
- **Timeline:**
  1. registration commit 73aef782 (2026-10-07 15:17Z);
  2. corrected checkpoints committed and stamped in e9102501 (15:29Z);
  3. rerun attempt 2 started 2026-10-08 07:07Z.

  Attempt 1 began about the same time as e9102501. Its log says the memory reaper killed it after about 10.5 h with no outputs written. `results/e427/evals/` contains no extra files.

## (b) Registered evaluator invoked unchanged — PASS-with-notes

- **No evaluator changes after the e427 registration.** `git diff 73aef782 HEAD` over the evaluator stack shows changes only in `experiments/e425_execreceipt.py`, and only in `main()`. The stack is `e422_protocol.py`, `e422_evaldriver.py`, `e424_evalgate.py`, `e422_manifest.py`, `e420_census.py`, and `e425_execreceipt.py`.
- **Not since the original batch-4 evaluation either.** The same holds against 9d8e9f28, the original e422 batch-4 evaluation commit. `run_logged_evaluation`, the gate, the driver, and the protocol are byte-unchanged.
- **What changed in `main()`:**
  - the CLI defaults now point to `results/e427`;
  - a pause guard was added;
  - an `E427_FROZEN_LEDGER_DIR` refusal was added.

  The wrapper does not call `main()`.
- **Changes since the e422 freeze.** Since 9dd2380b, `e422_protocol.py` changed in one commit, ca8e3ea3, two minutes after the freeze. That commit replaced only the `REGISTRATION_CITATION` string; this was the gate-(a) bind. The driver and manifest code are unchanged since the freeze. The gate commits from 2026-09-10 predate the original batch-4 run.
- **No uncommitted changes.** `git diff --quiet e9102501 HEAD` and `git diff --quiet HEAD` are both clean for the stack and the wrapper. The code on disk during the run was therefore the committed code.
- **The wrapper.** `experiments/e427_rerun_p5.py` (committed in e9102501, before the run) does three things:
  1. imports `e425_execreceipt`;
  2. calls `run_logged_evaluation(REGISTERED_T0, 4, "results/e427/batches", "results/e427/evals", prior_binding_path="results/e427/batches/e422_binding_state.json")` with seeds = None, which the driver defaults to `range(200)`;
  3. prints the two P5 bindings.

  It only redirects paths. It bypasses `main()`, and therefore the pause guard; registration §4.3 authorizes this one evaluation during the pause.
- **Note:** the registration names the CLI `experiments/e425_execreceipt.py --batch 4`. The actual call went to the same function through a wrapper. This invocation-route deviation is disclosed in the commit message. It does not change the computation, because the gate, the driver, and the T0 pin are all still exercised.

## (c) Independent recomputation of P5-90 / P5-95 — PASS (exact match)

**Implementation.** `phaseB_subagent/recompute_p5.py` imports no `experiments/` or `moluq` module. It is written from `docs/e422_registration.md` §3, §3A.1, and §7:
- it loads the corrected checkpoints 1–4 directly and checks their digests itself;
- eligibility is all four channels finite and pLDDT in [0, 100.01);
- the pool is `sorted(set(accessions))`;
- the split uses `default_rng(seed).permutation(n)` with t = n//2, c = (n−t)//2, so 70/35/36 accessions for n = 141;
- the center is the median of the train rows;
- the score is |y − center|;
- q is the ⌈(n_cal+1)(1−α_c)⌉-th order statistic;
- the interval is center ± q;
- coverage is the mean over test rows;
- MCSE = sd(ddof=1)/√200 and threshold = level − 3·MCSE.

I read the evaluator source, without importing it, to pin two details the registration leaves implicit:
- the interval is closed (`lo <= y <= hi`, with lo = center − q and hi = center + q);
- α_c is the float `1.0 − level`.

For every seed and level, the float-arithmetic k equals the exact-rational k. Runtime is 1.3 s.

| | Artifact (evaluator) | Independent recompute | Δ |
|---|---|---|---|
| P5-90 mean | 0.894492567684754 | 0.894492567684754 | 0 (bitwise) |
| P5-90 sd | — | 0.04074192962552704 | |
| P5-90 MCSE | 0.0028808894716835266 | 0.0028808894716835266 | 0 |
| P5-90 threshold | 0.8913573315849495 | 0.8913573315849495 | 0 |
| P5-90 outcome | HOLD | **HOLD** (margin +0.00314 = 1.09 MCSE) | |
| P5-95 mean | 0.9455777057113146 | 0.9455777057113146 | 0 (bitwise) |
| P5-95 sd | — | 0.024048384118961957 | |
| P5-95 MCSE | 0.0017004775487096876 | 0.0017004775487096876 | 0 |
| P5-95 threshold | 0.9448985673538709 | 0.9448985673538709 | 0 |
| P5-95 outcome | HOLD | **HOLD** (margin +0.00068 = 0.40 MCSE) | |

- **Per-seed comparison:** not possible. The registered artifact stores only aggregates (mean, MCSE, r_eff) and no per-seed pooled coverage. My per-seed values are in `phaseB_subagent/recompute_p5.json`. Because the 200-seed mean, sd, and MCSE match exactly, the per-seed values must agree at least to floating-point precision. A direct per-seed check was not possible.
- **Implementation cross-check.** I ran the same code on the frozen **original** ledger (`results/e422/batches`: 117,905 rows, 138 accessions). It gives 0.8579652901901389 (FAIL; threshold 0.8759801398774404) and 0.9151961367793325 (FAIL; threshold 0.9328446915697971). Both are bitwise equal to the anchored original-ledger P5 in `results/e422/evals/e422_batch04_evaluation.json`.
- **Sensitivity (not registered, descriptive only).** With an open interval (`lo < y < hi`) on the corrected ledger, both levels still HOLD:
  - P5-90: 0.89387 vs threshold 0.89131;
  - P5-95: 0.94534 vs threshold 0.94487.

  The registered closed-interval convention is therefore not what decides the 0.95 outcome, though that margin stays below 1 MCSE (`phaseB_subagent/sensitivity_p5.json`).
- **Single-seed range** (descriptive): 0.787–0.969 at 0.90 and 0.875–0.985 at 0.95.

## (d) Binding scope — PASS

- **Sidecar.** `results/e427/batches/e422_binding_state.json` (schema `e422-binding-state-v1`, batch_index 4) has binding keys **exactly `p5_90` and `p5_95`**, both `bound_at_batch: 4`, with no other keys. The sidecar was a new file in f79778d2, so this was a fresh batch-4 chain start, as the gate requires for batch 4.
- **Evaluation artifact.** The binding view shows p6, p7a, and p7b with outcome string `DESCRIPTIVE_AT_BATCH4_FIRST_EVALUABLE_BATCH8`. The inner predictions contain only `p5`. The per-level aggregates carry no `per_bin` key, and the inner evaluation has no `joint_grid`. This is consistent with `_scope_batch4` suppression. The aggregate *key names* include `r_bin` and `mondrian_deferred_bins_tally`; their values were not read.

## Negative / residual findings (non-blocking)

1. **Not blind (§4.4).** The P5-95 HOLD margin is 0.40 MCSE. A different but equally defensible label repair could plausibly flip it. Registration §4.4 fixes the repair by data format, not by outcome, and Phase A item 3 independently confirmed the corrected labels. However, the outcome of record rests on a margin smaller than its own Monte-Carlo error. It should be reported with that caveat, alongside the original-ledger FAIL.
2. **Invocation route.** The rerun called the registered function through a wrapper instead of the registered CLI, bypassing the pause guard by design. The computation is unchanged; see (b).
3. **No per-seed record.** The artifact stores no per-seed P5 values, so a per-seed audit can only go through recomputation. This is not a defect, because the registered schema does not include them.
4. **OTS stamps** on the three Phase B artifacts are still pending Bitcoin confirmation.
5. **Corrected population.** The corrected ledger adds 5,549 rows and 3 accessions (141 vs 138) relative to the original. Phase A item 7 already notes this change in the evaluated population.

## Verdict

**PASS-with-notes.** No blocking item. All three Phase B digests and all four corrected-checkpoint digests verify under an independently re-implemented canonical rule. The evaluation reads the corrected `results/e427/batches` ledger (123,454 rows, 141 accessions). The evaluator stack is unchanged since the registration and since the original batch-4 run, and the wrapper only redirects paths. An independent, non-importing recomputation reproduces P5-90 = 0.894492567684754 (HOLD) and P5-95 = 0.9455777057113146 (HOLD) bitwise. Only P5 is bound, and P6/P7a/P7b remain descriptive and unbound. Notes: the rerun is not blind (§4.4), and the P5-95 margin (0.40 MCSE) is smaller than its MCSE.

## Files written

- `results/e427/audit/PHASE_B_independent_subagent.md` (this report)
- `results/e427/audit/phaseB_subagent/verify_digests.py`, `digests.json`
- `results/e427/audit/phaseB_subagent/recompute_p5.py`, `recompute_p5.json` (per-seed pooled marginal coverage, seeds 0–199)
- `results/e427/audit/phaseB_subagent/sensitivity_p5.py`, `sensitivity_p5.json`
