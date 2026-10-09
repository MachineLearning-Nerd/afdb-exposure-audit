# Independent audit of e427

**Phase A verdict (items 1–7): PASS-with-notes. Phase B verdict (item 8): PASS-with-notes. Overall verdict: PASS-with-notes.** No blocking audit finding remains. Three OpenTimestamps proofs are present; verification reports pending Bitcoin confirmations and no digest mismatch. **The e422 pause remains in force; this audit did not change it.**

Audit dates: 2026-10-07 (Phase A), 2026-10-09 (Phase B). Work was read-only outside this directory. No commit was made; `paper/` and e422 artifacts were not touched. The pause file remained present at the end of these checks.

## Phase A findings

### 1. Repair code and tests — PASS-with-notes

`experiments/e427_label_repair.py` maps SIFTS `residue_number` as `label_seq_id` and binds atoms by `struct_asym_id`; it has no author-number fallback. Equal-span segments use the SIFTS offset. Unequal spans use global sequence alignment of `_entity_poly_seq` to the AFDB model sequence, retaining identical pairs. It rejects a segment above the registered 10% mismatch fraction and excludes mismatching residue pairs below that threshold.

The tests include the old 9qj6 register shift, corrected real 9qj6 bytes, a synthetic indel, point mutations, missing `struct_asym_id`, altloc/sequence fixtures, CLI pause/path refusals, and replay transport hash/order behavior. Results:

- `tests/test_e427_label_repair.py`: **15 passed**.
- `tests/test_e42*.py`: **553 passed, 2 skipped** in 123.68 seconds.

Edge-case review: label numbering avoids author insertion-code shifts; MSE maps to methionine through the unchanged e421 residue map; atom parsers accept altloc IDs `.`, `?`, blank, and `A`. `?` remains accepted under the inherited rule, but no `?` Cα occurred in the 36-segment sample (which contained `.` and `A`). The raw-mapping scan found eight multi-segment accession/entry/chain groups, with no overlapping label or UniProt intervals. All eight are explicitly excluded as `E422_ABSENT_FROM_SNAPSHOT`, so those current records do not exercise multi-segment labeling. In the independent sample, there were no insertion-coded Cα rows, MSE rows, repeated label positions, or entity-sequence variant positions; altloc A occurred in 26 sampled segments. The sampled multi-asym case (A0A4Q0WMG2) maps author chains A/B to struct-asym IDs C/A and reproduces exactly.

One auditability note: the repair function returns `n_identity_mismatch`, but `e423_runner.run_batch` retains only `labeled["rows"]` on successful segments. Thus the per-segment mismatch counts used by the repair are not persisted in the corrected checkpoint. The filtering itself passed the independent sample; persisting these counts would make later review easier. For the D1MPT3 sample, all four segments used sequence alignment and matched the corrected rows; Q9HWK6 had one identity mismatch per chain, also excluded consistently.

### 2. Replay self-check — PASS-with-notes

The recorded `results/e427/replay_report.json` has verdict `OK`. For batches 2–4, every frozen-labeller self-check is true: rows, exclusions, census counts, mapped-entry counts, manifest additions, verified accessions, manifest checkpoint digests, and receipt consumption. The corrected pass has the same manifest digest as the original and leaves zero receipts unconsumed.

I reran `experiments/e427_replay_relabel.py` from `/tmp`, using the raw e422 receipt cache, with outputs under `/tmp/e427_audit_replay_checkpoints_20261007/` and the rerun report under this audit directory. Every corrected batch row list is identical to the recorded corrected checkpoint. Batch 1 is byte-identical to the original. For batches 2–4, full checkpoint bytes and SHA-256 values differ because `e427_relabel.replay_utc` is regenerated and the canonical checkpoint digest consequently changes. After removing only that replay timestamp and the dependent `sha256`, every field is identical. The replay report and comparison are in [replay_rerun_report.json](replay_rerun_report.json) and [independent_label_sample.json](independent_label_sample.json).

The replay transport keeps receipt attempt order within each `(method, URL)` queue and verifies cached raw-byte hashes. It does not assert a single global order across different URLs; all recorded receipts were consumed and the frozen replay reproduced the recorded rows and manifest. This is a traceability note, not a data mismatch.

### 3. Corrected checkpoints and independent labels — PASS

Each corrected checkpoint has a valid canonical digest, the original receipts, the original manifest digest, and an `e427_relabel` provenance block whose original-checkpoint hash matches. Batch 1 is byte-identical to the original checkpoint. Corrected row counts are:

| Batch | Original rows | Corrected rows | Change |
|---:|---:|---:|---:|
| 2 | 29,054 | 31,868 | +2,814 |
| 3 | 45,405 | 45,903 | +498 |
| 4 | 43,446 | 45,683 | +2,237 |
| **Total** | **117,905** | **123,454** | **+5,549** |

Using a separate standard-library mmCIF/SIFTS parser and independent NumPy Kabsch/lDDT calculation, I relabeled **36 raw segments**: 26 targeted segments and 10 seed-427 sampled unchanged segments. The targeted sample includes K7PQ54 (9qj6, 9qj7, 9qj8, 9qjd, 9qjf), all four Q5JF22/9sq1 chains, P07342/9y8y, A0A4Q0WMG2/9wvj, all four Q9HWK6/31rm chains, and both chains from D1MPT3 entries 31al and 31jt. All selected mapping, reference, and AF model receipts were hash-checked. **All 36 row sets and values match the corrected checkpoint to its four-decimal precision; maximum numeric difference was 0.**

The corrected K7PQ54 rows number 17,879 across 38 entry-chain groups, versus 16,947 original rows. For 9qj6 specifically, the independent map yields 472/471 aligned Cα rows for A/B and median lDDT 0.9929/0.9950. The former low-score rows are not retained under the corrected mapping.

For 9zxa, the cached mapping response is now parseable and contains eight Q92802 mappings. No successful reference-CIF receipt or reference bytes exist. The corrected checkpoint records eight `E422_REF_CIF_FETCH_FAILED` failures with `status=None`; it produces no fabricated residue rows.

Across 366 hash-verified SIFTS responses, the scan covered 829 accession-entry-chain groups. Eight had multiple disjoint mapping segments; none had overlapping label/UniProt intervals. The current e427 manifest excludes those eight accessions, so their multi-segment labeling behavior was not exercised by the corrected rows.

### 4. Dispatch/evaluation rewire and pause guard — PASS

The dispatch entry point passes `parse_mappings_labelspace`, `label_triple_labelspace`, and the e427 labeller ID. `run_batch` keeps its previous defaults when those optional callables are omitted, which is how the frozen replay self-check runs. The evaluation defaults point to `results/e427/batches`, `results/e427/evals`, and the e427 binding sidecar; paths under `results/e422/` are refused.

`results/e422/PAUSE.json` was present. With the live guard, dispatch and evaluation both returned `E427_PROGRAM_PAUSED` before side effects. In an isolated Python process only, I redirected the guard’s in-memory path to a nonexistent `/tmp` path and confirmed that dispatch and evaluation then refused the frozen `results/e422/batches` path with `E427_FROZEN_LEDGER_DIR`. No ledger path was written.

### 5. e421 relabel — PASS-with-notes

All 466 successful HTTP response receipts (of 469 total) match the SHA-256 and byte length of their content-addressed raw files. The three HTTP 404s are recorded AFDB v4 URLs; their response-body hashes and lengths are in the receipt log, but those error bodies were not persisted and therefore could not be re-hashed. The 404 bodies were not calculation inputs. All 410 files in `data/e427/e421_raw/` are receipted and hash-named, and no successful URL has multiple body hashes.

I independently recomputed the old and corrected summaries from their row-level inputs. Every reported field matches: status counts and transitions, successful entries/accessions, rows, per-cell row and unique-residue counts, pLDDT bins, error rates, Spearman correlation, and finite sample sizes. Corrected values are 168 entries, 139 `ok` entries / 85 accessions, 37,531 rows, cells A/U/OC/O/B = 5,047/31,765/354/363/2, and Spearman 0.6430108209634237 (n=37,531). The original recomputation matches 96 `ok` entries / 64 accessions, 24,457 rows, and Spearman 0.6170601466142406 (n=24,456). Thus the e421 descriptive population changes by +43 successful entries and +13,074 rows; it remains exploratory and B remains sparse (2 rows).

The `kabsch`, `lddt`, and `bfactor_z` function ASTs in `e427_e421_relabel.py` are identical to the frozen copies in `e421_census_fixed.py`. The replacement `cell_of` agrees with e421's inline cell rule over the tested threshold truth table (45/45 combinations). Independent raw-byte spot checks cover every changed-status category: 10af/10ah and alignment-path 10qf (`insufficient_aligned -> ok`) reproduce 168, 355, and 500 rows exactly; 10cy reproduces `mapping_acc_absent`; 11ns/11ny/11pu reproduce identity failures at 151/158, 148/156, and 154/170 mismatches; and 10pa/10ps reproduce `ok -> sequence_identity_failed` at 58/316 and 58/555 mismatches. The five unchanged `insufficient_aligned` status group was also sampled (10dv, zero aligned pairs). All recomputed status transitions match the summary.

One implementation-description note: the e427 script header calls its altloc handling a deviation because the e421 parser “kept” other altlocs. The frozen e421 parser actually filters to the same accepted IDs (`.`, `?`, blank, `A`) used by the e427 reference parser. e427 additionally suppresses repeated UniProt positions with `seen`; a scan of all 149 fetched, mapped reference entries found no duplicate accepted-altloc Cα label positions, so this difference did not affect these outputs. Two entries (10ee and 10ef) have four chain-A mapping segments each; their label and UniProt spans are disjoint, and no selected accession has multiple struct-asym IDs. The source-drift limitation remains material: September e421 source bytes were not retained, so this audit validates the new receipted run and its calculations, not byte identity with September sources.


### 6. Registration, correction record, and cron prompts — PASS-with-notes

The registration and its `.ots` sidecar are present. `ots verify -f docs/e427_registration.md docs/e427_registration.md.ots` found no document/hash mismatch, but returned exit 1 because calendar attestations are still waiting for six Bitcoin confirmations. The timestamp is not yet fully confirmed.

The correction-of-record document preserves the original P5 result as a label artifact and discloses the earlier offline diagnostic. The old and e427 cron schedules have the same fire times and purposes. Dispatch/evaluation output paths now use `results/e427/`; raw bytes and the manifest remain under `results/e422/`. The weekly prompt explicitly checks PAUSE, and both dispatch/evaluation CLIs enforce the guard before side effects. Batch-8 and batch-12 prompts rely on those CLI guards rather than repeating a prompt-level PAUSE check. The readiness probe remains informational and has no write or dispatch step.

### 7. Adversarial review — PASS-with-notes

The corrected ledger is intentionally a different label set from the original: it adds 5,549 rows across batches 2–4. That changes residue weighting and the evaluated population of labeled residues; it is not a byte-preserving correction of the original outcome. The registration is explicit that the corrected-ledger P5 result becomes the outcome of record while the original remains alongside it.

The rerun is not blind. Registration §4.4 discloses an earlier offline diagnostic of 0.8980/0.9474 against thresholds 0.8913/0.9448. The repair rule, identity threshold, seeds, evaluator, and thresholds are fixed in the registration. Report Phase B as a registered deterministic rerun on the corrected ledger, with the prior diagnostic disclosed; do not describe it as a blind replication.

The e422 replay avoids source drift by using the original receipted raw bytes. The e421 relabel re-fetches public source files, and the registration discloses that the September bytes were not retained. Phase B and item 5 therefore answer different provenance questions and should not be merged into one claim.


## Phase B findings

### 8. Registered P5 rerun — PASS-with-notes

The final evaluation and execution receipt are byte-identical to the versions in commit `f79778d2`. The receipt links to batch 4, the registered T0 and role, and the evaluation timestamp; it records completion at **2026-10-09 03:33:17 UTC (09:03:17 IST)** after 73,550 seconds (20.43 h), within the registered 24 h budget. The final log emits both P5 binding records. The separate attempt-1 log says that attempt was killed after about 10.5 h before writing outputs, then restarted as a systemd service; it does not replace or contaminate the completed run.

`experiments/e427_rerun_p5.py` directly calls `e425.run_logged_evaluation(REGISTERED_T0, 4, "results/e427/batches", "results/e427/evals", prior_binding_path="results/e427/batches/e422_binding_state.json")`, as disclosed in the correction-of-record addendum because the CLI pause guard refuses the CLI. The wrapper then runs the e424 digest gate and the registered e422 driver. It does not override `seeds`, so the frozen driver uses seeds 0–199. The `run_logged_evaluation` function AST matches the registered baseline commit `73aef782`; the frozen `e422_evaldriver.py`, `e422_protocol.py`, and `e422_manifest.py` are byte-identical to that baseline.

All internal canonical SHA-256 digests recompute exactly. Evaluation `f579c541…`, execution receipt `1cf9a3e1…`, and binding sidecar `c3ee0f6d…` match their stored digests. The evaluation ledger equals the corrected checkpoint ledger, including each canonical digest and row count:

| Batch | Rows | Checkpoint SHA-256 prefix |
|---:|---:|---|
| 1 | 0 | `4f54057d…` |
| 2 | 31,868 | `a18a2274…` |
| 3 | 45,903 | `d0677b47…` |
| 4 | 45,683 | `ad6166b5…` |
| **Total** | **123,454** | |

All three `.ots` stamps were checked with `.venv/bin/ots verify`. Each proof recognizes the file but reports pending Bitcoin confirmation from its calendars; there is no digest-mismatch report. Thus file integrity is verified, while public timestamp finality is still pending.

I independently recomputed only pooled P5 from the corrected checkpoint rows, without importing any e42x or moluq module. Eligibility produced 123,454 rows across 141 accessions; each of seeds 0–199 used the registered accession-disjoint 70/35/36 train/calibration/test split. For each level I used the training-residue median, the registered absolute-residual conformal order statistic, inclusive test-interval coverage, and the sample MCSE across all 200 repeats. Both outcomes exactly match the evaluation artifact and binding state:

| Prediction | Artifact mean / independent mean | MCSE | Registered threshold (nominal − 3 MCSE) | Margin above threshold | Outcome |
|---|---:|---:|---:|---:|---|
| P5-90 | 0.894492567684754 / 0.894492567684754 | 0.002880889471684 | 0.891357331584950 | 0.003135236100 | HOLD |
| P5-95 | 0.945577705711315 / 0.945577705711315 | 0.001700477548709 | 0.944898567353871 | 0.000679138357 | HOLD |

P5-95 is a narrow registered HOLD: the margin is about 0.40 MCSE. The run is deterministic and reproducible, but not blind, as the earlier diagnostic is disclosed in Phase A. The batch-4 binding sidecar contains only the two P5 keys; P6/P7 have no bound entries.

## Audit artifacts

- [independent_label_sample.py](independent_label_sample.py)
- [independent_label_sample.json](independent_label_sample.json)
- [replay_rerun_report.json](replay_rerun_report.json)
- [independent_e421_audit.py](independent_e421_audit.py)
- [independent_e421_audit.json](independent_e421_audit.json)
- [independent_p5_audit.py](independent_p5_audit.py)
- [independent_p5_audit.json](independent_p5_audit.json)

The independent label and P5 recomputations use no e427, e42x, or moluq module imports. All audit report and script/JSON outputs are under `results/e427/audit/`; the one replay destination was the `/tmp` directory required by the brief.

## Overall verdict

**PASS-with-notes; no blocking audit item.** The corrected-ledger P5-90 and P5-95 outcomes reproduce exactly, both are registered HOLDs, all linked ledger and artifact digests verify, and the batch-4 sidecar has no P6/P7 bindings. The P5-95 margin is narrow and the OTS proofs still await Bitcoin confirmations. The e422 pause remains present; lifting it is outside this audit task.
