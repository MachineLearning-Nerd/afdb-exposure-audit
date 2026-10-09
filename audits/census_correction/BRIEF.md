# LeakagePaperChecker — independent audit brief: e427 (label repair, corrected ledger, rewire, e421 relabel)

Independent, fresh-context audit required by `docs/e427_registration.md` §4.6
before the e422 pause can be lifted. Work read-only on the repository except
for your outputs under `results/e427/audit/` (scripts, JSON, report). Do not
commit, do not modify experiments/, docs/, tests/, results/e422/, or other
results/e427 files, do not touch paper/. Use `python -I` from a cwd outside
the repo where possible, and do not import the e427 modules for the
independent recomputations (you may run them and the test suite).

## Phase A (now)
1. **Repair code** — `experiments/e427_label_repair.py`: confirm label-space
   mapping (label_seq_id + struct_asym_id, no author fallback), the indel
   alignment path, and the identity check behave as the registration §2
   says; look for remaining numbering hazards (insertion codes, microheterogeneity,
   multi-segment chains, chains with several struct_asym ids, MSE, altlocs).
   Run `tests/test_e427_label_repair.py` and the full `tests/test_e42*.py`.
2. **Replay self-check** — `experiments/e427_replay_relabel.py`,
   `results/e427/replay_report.json`: verify the replay serves original
   receipts in order with hash checks and that the frozen-labeller pass is a
   genuine reproduction (rows, exclusions, manifest sha). Re-run it into a temp
   dir (`--out-dir /tmp/...`) and confirm byte-identical corrected checkpoints.
3. **Corrected checkpoints** — `results/e427/batches/e422_batch0{1..4}_checkpoint.json`:
   digests valid under the registered canonical rule; batch 1 identical to the
   original; `receipts` equal the originals; provenance block present.
   Independently relabel a sample of ≥ 30 segments from `results/e422/raw`
   with YOUR OWN label-space mapper (e.g. the one you wrote for 9qj6), including
   K7PQ54, Q5JF22, P07342, A0A4Q0WMG2, Q9HWK6, D1MPT3 (indel case), and ≥ 10
   random unchanged segments; compare rows. Check 9zxa's handling (newly
   parseable entry whose reference CIF was never fetched).
4. **Rewire** — `experiments/e423_dispatch.py`, `experiments/e423_runner.py`,
   `experiments/e425_execreceipt.py`: dispatch passes the e427 labeller and
   refuses `results/e422/batches`; evaluation defaults to `results/e427` and
   refuses `results/e422/*`; pause guard refuses before side effects; default
   `run_batch` behaviour unchanged.
5. **e421 relabel** — `experiments/e427_e421_relabel.py`,
   `results/e427/e421_relabel/{census_results,summary,receipts}.json`,
   raw bytes `data/e427/e421_raw/`: receipts hash-match bytes; e421 label
   rules preserved (verbatim functions); recompute the summary numbers
   independently; spot-check entries whose status changed.
6. **Records** — `docs/e427_registration.md` (+ .ots),
   `docs/e422_correction_of_record_20261007.md`, `docs/e427_cron_prompts.md`
   (every path now results/e427/batches; PAUSE check; nothing else changed vs
   `docs/e422_cron_prompts_20261007.md` — diff them).
7. **Adversarial pass** — strongest objections to treating the corrected
   ledger as the e422 ledger of record (post hoc concerns, disclosed
   diagnostic, selection effects from newly admitted segments, source drift).

## Phase B (when `results/e427/evals/` exists and `results/e427/rerun_p5.log` shows completion)
8. **Registered P5 rerun** — verify the evaluation artifact, execution receipt,
   and binding sidecar digests; confirm it used the unchanged registered
   evaluator on `results/e427/batches`; independently recompute P5-90/P5-95
   (pooled interval, seeds 0–199) from the corrected checkpoints; confirm no
   P6/P7 quantity is bound.

## Output
`results/e427/audit/AUDIT_E427.md`: verdict line for Phase A now (PASS /
PASS-with-notes / FAIL, with any blocking item), then append Phase B when it
runs. Facts and recomputed numbers side by side; negative findings matter.
