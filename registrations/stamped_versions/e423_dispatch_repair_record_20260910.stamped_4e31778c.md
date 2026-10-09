# e423 dispatch-machinery repair — record + delta registration (2026-09-10, w1:pG)

**Allocation.** e423 allocated 2026-09-10 by w1:pG under main's standing
delegation (DECISIONS 2026-09-09T20:31+05:30).  Trigger: the gate-(e)
dispatch-preparation shakedown (below) exposed three defects in the
frozen e422 dispatch implementation.  Per the gate-(a) freeze discipline
(docs/e422_freeze_gate_a_20260909.md: any change to a bound object is a
new E-number), the repair is this new E-number.

**Inheritance (what e423 does NOT change).**  The e422 FROZEN
registration text — git-blob 6e1b9ce128ac3a1be4678463bdc777c25b50ca8d,
content sha256 dc7a362773fbf3a03e759aa7d178125963a14d88d9bbd61b66389a0
4ba8d7b3b, freeze commit 9dd2380b, OTS-anchored (CONFIRMED 2026-09-10:
BitcoinBlockHeaderAttestation(966254), shared txid 85d3598a…, receipt
docs/e422_anchor_receipt_gate_a_20260909.json) — is UNCHANGED
and remains the registration of record for the study.  T0 =
2026-09-11T00:00:00Z unchanged.  The four other bound e422 modules
(protocol, manifest, status, evaldriver) are UNTOUCHED — their frozen
content SHAs stand.  The defect was confined to e422_runner.py (bound
at content sha256 ac7729a5…) and its thin CLI, plus the shared
e420_census.py scheduler the runner calls into.

## The three defects (found live, never silent in effect — the shakedown
## caught all three BEFORE any registered batch window opened)

1. **Invalid census request_options (defect in e422_runner.py
   build_census_query).**  `"scoring_strategy": "none"` is not in the
   RCSB search v2 enum; the API returned HTTP 400 on EVERY page.
   Live-verified repair: the field is dropped (scoring-only, irrelevant
   to a metadata predicate); HTTP 200 confirmed 2026-09-10.
2. **Silent-drop on non-2xx (defect in the scheduler fall-through,
   e420_census.py RequestScheduler.get and the e422 POST subclass).**
   Any status outside the 429/5xx retry class (400/404/…) was marked
   ok=True and its ERROR BODY parsed as data — the shakedown produced a
   well-formed n_rows=0 checkpoint from an HTTP 400 error document,
   exit 0.  The worst possible failure mode for a prospective
   registration: an empty study indistinguishable from a quiet week.
   Repair: non-2xx statuses are typed failures
   (E420_HTTP_UNEXPECTED_STATUS, receipted, body head in the failure
   detail); the e423 runner additionally validates the census document
   strictly (dict; integer total_count >= 0; result_set list of
   identifier dicts) with the new typed code
   E423_CENSUS_RESPONSE_MALFORMED — a legitimate empty window
   (total_count 0) is only ever accepted on a VALIDATED 200.
3. **Parser/fixture echo (defect in e422_runner.py parse_mappings and
   its test fixture).**  The parser assumed mapping items carry a
   nested `"chain": {"chain_id": …}`; the live PDBe v2 payload is FLAT
   (`chain_id` at top level, with `struct_asym_id`) — exactly as the
   field-proven e420 resolve_chain and e421_pipeline parsers read it.
   The nested shape existed only in the parser's own test fixture
   (KeyError 'chain' on the first live receipt; and the KeyError
   propagated UNTYPED mid-batch).  Repair: flat-shape parser; every
   structural surprise in a mapping item raises the typed
   E422_MAPPING_PARSE_FAILED and is tallied per-entry, never crashes a
   batch.  Fixtures now mirror receipted live bytes.

**Shared-infra 404 semantics preserved.**  The e420 census pipeline's
exclusion semantics (AFDB prediction absent → excluded with
availability recorded; segment/primary unavailable) previously reached
their 404 branches because the buggy scheduler marked 404 ok.  Under
the repaired scheduler the 404 classification is hoisted ABOVE the
hard-failure gate in fetch_primary / fetch_afdb / _process_entry
(tested by the unchanged test_afdb_prediction_absent_excludes_without_failure).

**Supersession.**  All dispatch execution — the registered batch
dispatch, gate-(e) — runs through experiments/e423_runner.py +
e423_dispatch.py from now on.  The e422 batch-checkpoint artifact
schema, filenames (e422_batchNN_checkpoint.json / _freeze.json),
allowlist lane roles, and study identity are unchanged (the e422
evaldriver consumes the checkpoints as-is); checkpoints record runner
provenance "e423 (dispatch-machinery repair; DECISIONS 2026-09-10)".

## Review + tests

- Detached review (2026-09-10): TWO independent detached read-only
  auditors over the e423 diff (commit bdc1d0a9) against this record —
  auditor A (adversarial code lens: dynamic execution of the repaired
  paths, SHA recomputation, shakedown-artifact re-verification) and
  auditor B (governance / registry-consistency lens: freeze-discipline
  honor, bound-object integrity, ledger verification, live re-runs).
  Verdicts + findings + the disposition are below.
- Test state AT bdc1d0a9 (precise, correcting the commit message):
  test_e422_* suites 117 (runner 21, dispatch 5, protocol 43, manifest
  26, evaldriver 17, status 5); test_e423_* suites 38 (runner 33,
  dispatch 5); shared-infra test_e420_census.py 20 — the
  repair-relevant total is 175 passing.  The bdc1d0a9 commit message
  said "172 e42x" — an arithmetic slip (off by three); reviewer A's
  recount of "155" also missed test_e420_census.py.  All e420/e422/e423
  test files at bdc1d0a9 collect 471 (469 passed, 2 skipped
  pre-existing e420 pool skips).  Full repo suite at bdc1d0a9: 1382
  passed / 2 skipped (the +4 delta vs the earlier working-tree count
  came from unrelated pre-existing edits to tests/test_ratios.py and
  tests/test_e418_post_auth_*.py — not e423 content).  The one e422
  test that pinned the old silent-ok 404 path was AMENDED with a dated
  comment — the frozen e422 module bytes are unchanged; tests are not
  frozen artifacts.

## Live shakedown receipts (dev-only; machinery validation, NOT evidence)

Non-registered historical window (2026-08-28T00:00:00Z,
2026-09-04T00:00:00Z], --max-entries 25, scratch dirs
results/e423_dev_shakedown/ (dev-only paths; exploratory discipline —
no headline numbers, no claim promotion, ledger never cites these rows
as results):

- Round 1 (frozen e422 runner, commit aa6eb3e8 lineage): HTTP 400 every
  page → silent empty checkpoint n_rows=0 — defects 1+2 reproduced.
- Round 2 (partial repair, dispatch still importing e422_runner): same
  silent empty — caught the stale import before it could matter.
- Round 3 (e423 complete, commit bdc1d0a9): census 200 → 25 entries →
  22 with mappings (2 MAPPING_UNAVAILABLE typed, 1 NO_PROTEIN_MAPPING)
  → 12 accessions manifest-admitted → 5,852 rows labeled, 7 accessions
  verified batch-time → typed exclusions: 21 × E422_ABSENT_FROM_SNAPSHOT
  (fan-out of the one E422_CHILD_URL_NOT_V6 accession — registered
  fail-closed), 1 × E422_INSUFFICIENT_ALIGNED.  75 receipted fetches;
  checkpoint digest verifies.  Checkpoint:
  results/e423_dev_shakedown/batches/e422_batch01_checkpoint.json —
  internal canonical payload sha256 0e580c7fac2cbba9… (the checkpoint's
  own "sha256" field), file sha256
  cf24affb27630131eb2e95228168adc58d18e8285d66f1a85aacda36bb4606d9
  (863,925 bytes); freeze receipt file sha256
  826bfc60cec20aedc98c1caeb911312a2e932b7445a8fbf018340908443217b3;
  manifest checkpoint file sha256
  09cf02fdecc28bbb04ceb0d4a52c759054f71d1dfce2b8b826fa4477c6149dbf.
  results/ is gitignored; these DEV-ONLY artifacts are deliberately
  UNTRACKED (never to be mistaken for registered artifacts) — the
  digests above are their record.  Raw fetch bytes (29 MB) also
  untracked; every receipt (URL, status, sha256 of bytes) is inside
  the checkpoint's receipts array.

## Dispatch authority (gate (e))

Batch 1 opens 2026-09-11T00:00:00Z over (T0, T0+7d].  Per main's
dispatch authorization (DECISIONS 2026-09-09), batch dispatch decisions
are delegated to w1:pG under the standing governance rules.  The
shakedown validated: census query validity, pagination, typed-failure
paths, manifest fail-closed admission, freeze-receipt-before-labels
ordering, four-channel row schema, checkpoint digest discipline.

## Detached review — verdicts + findings (2026-09-10)

**Auditor A (adversarial code lens) — READY_WITH_FINDINGS: 0 BLOCKER /
2 MAJOR / 5 MINOR / 4 NOTE.**  Dynamically confirmed: the
scoring_strategy omission is complete; typed non-2xx handling is
complete with no receipt loss; e422_runner.py bytes still match the
frozen content SHA ac7729a5…; e420_census.py is bound by nothing (the
shared-infra change is freeze-legal); the frozen evaldriver consumes
e423 checkpoints unchanged; 404 exclusion semantics are preserved; all
shakedown artifact digests verify (5,852 rows; 75 receipts: 73×200 +
2×404 typed); the batch-1 window is NOT invalidated.

- A-MAJOR-1: the superseded e422_dispatch.py remained a live silent-fail
  path at HEAD — one command from reproducing the empty-checkpoint
  defect; a refuse-to-run guard is freeze-legal (not a bound object).
- A-MAJOR-2: parse_mappings left five untyped crash classes OUTSIDE the
  per-item try — the loop head, accession-value checks, and the res_map
  arithmetic (non-dict accession values, list/None UniProt block, null
  unp_end, non-numeric unp_start, non-list mappings) could crash a
  batch untyped.
- A-MINOR-3: metadata-lane untyped paths (UnicodeDecodeError on
  undecodable bodies; scalar/null/ill-typed metadata items; null
  latestVersion reaching the frozen module's int(None); a non-object
  verification payload reaching verify_against_manifest).
- A-MINOR-4: total_count accepted via int() coercion (int("7"),
  int(7.9), int(True) all passed); a page stream ending short of
  total_count on an empty page silently truncated the population.
- A-MINOR-5: test_typed_census_error_exits_nonzero was vacuous
  (simulated the guard inline, asserted nothing about it); the
  __main__ guard caught only CensusError.
- A-MINOR-6: record/commit test counts (see the correction above).
- A-MINOR-7: the DECISIONS e423 ALLOCATED entry was missing.
- A-NOTEs (4): verification-positive observations and documentation
  suggestions — accepted; the documentation items are covered by the
  notes below (e420 re-run semantics; tombstone ledger note).

**Auditor B (governance / registry-consistency lens) —
READY_WITH_FINDINGS: 0 BLOCKER / 2 MAJOR / 4 MINOR / 8 NOTE.**  All
SHAs recomputed live from git objects; all counts re-run live.

- B-F1 (MAJOR): the entire e423 event was unledgered — no DECISIONS
  entry existed while the amended test, both artifact writers, and this
  record cited "DECISIONS 2026-09-10" (dangling pointers).  Discharge:
  dated DECISIONS entries BEFORE gate (e)/T0.
- B-F2 (MAJOR): frozen-text/object conflict — the frozen registration
  names e422_runner.py as the committed implementation (§2 amendment;
  §6 bound set) while dispatch execution moved to e423 modules; the
  freeze discipline resolves exactly this by a ledger note, not by
  touching the frozen text.  Discharge: DECISIONS supersession note
  (with B-F13's delegation-scope interpretation).
- B-F3 (MINOR): this record was uncommitted and still mutating during
  the audit; its review gate read "wf_<pending>" with nothing appended.
  Discharge: this commit; review section filled above.
- B-F4 (MINOR): the bdc1d0a9 commit message said "172 e42x"; actual
  175 under the repair-relevant reading; "Full e42x" was loose (all
  e420/e422/e423 files collect 471).  Discharge: precise per-suite
  statement above.
- B-F5 (MINOR): the record's "1382 passed / 2 skipped" holds at HEAD
  bdc1d0a9, not in the then-current working tree (unrelated working-
  tree test edits explained above).  Discharge: both numbers now
  stated with their baselines.
- B-F6 (MINOR): "--max-entries … (recorded in the checkpoint)" was not
  literally true — no checkpoint field existed.  Discharge: explicit
  "max_entries" checkpoint field (None for a registered dispatch) +
  regression test.
- B-F7 through B-F12 (NOTEs, all PASS): freeze discipline honored
  (all five bound modules + frozen text byte-identical at bdc1d0a9 and
  HEAD); tests are not bound objects and the amendment is dated and
  attributed; e420_census.py is bound by no registration anywhere; the
  dev-only shakedown is fully gitignored-contained with no ledger
  citation; every shakedown receipt and number verified; the
  e422→e423 delta is exactly the three claimed repairs with no
  registered rule moved (MIN_ALIGNED, page size, page cap, resolution
  bound, window arithmetic, freeze-before-labels, filenames, schemas,
  lane roles, DEFAULT_T0 all byte-identical).
- B-F13 (NOTE): the delegation reads "for THIS registration only:
  E-number allocation …" — allocating a successor E-number under it is
  an interpretation that belongs in the DECISIONS entry.  Accepted.
- B-F14 (NOTE): nothing mechanically prevented gate (e) from running
  the frozen, known-broken e422_dispatch path (the frozen POST override
  bypasses the repaired shared scheduler).  Discharge: same as
  A-MAJOR-1 (tombstone).

## Review disposition (executed 2026-09-10, this commit)

| finding | action |
|---|---|
| A-MAJOR-1 / B-F14 | experiments/e422_dispatch.py rewritten as a refuse-to-run TOMBSTONE: main() unconditionally prints TYPED_FAILURE E423_DISPATCH_SUPERSEDED and exits 2 (even --dry-run, whose printed query is the invalid one).  Not a freeze-bound object; the freeze artifact's informational SHA row (12211ca7…) is superseded by this ledger note + DECISIONS. |
| A-MAJOR-2 | parse_mappings typed envelope extended over the WHOLE walk (isinstance dict check on the UniProt block; per-accession try covering loop head, accession-value/mappings-list type checks, residue-bound int strictness incl. bool rejection, and the res_map arithmetic); detail now carries the exception message; 7 new typed regressions (6 direct + 1 run_batch-tallied). |
| A-MINOR-3 | Stage-3 admission: decode errors typed; pre-validation (non-empty array of JSON objects) before the frozen module sees the payload; residual catch (TypeError/AttributeError/KeyError/IndexError/ValueError) → TYPED MANIFEST_SCHEMA_INVALID (covers the frozen int(None) latestVersion gap).  Stage-5 verification: same pre-validation + residual catch → TYPED METADATA_FETCH_FAILED.  4 new regressions. |
| A-MINOR-4 | total_count strict bool/int isinstance (no coercion); post-pagination cross-check: len(entries) < total_count → TYPED CENSUS_RESPONSE_MALFORMED "truncated pagination".  3 new regressions (bool, float, truncation). |
| A-MINOR-5 | Vacuous test replaced by two real runpy __main__-path tests (stubbed requests; HTTP 400 end-to-end → TYPED_FAILURE E422_CENSUS_QUERY_FAILED exit 2; residual RuntimeError → TYPED_FAILURE E423_UNEXPECTED_RuntimeError exit 2); __main__ guard widened with the typed-by-class clause. |
| A-MINOR-6 / B-F4 | Test counts corrected with per-suite precision (section above); commit-message slip recorded here rather than rewritten (bdc1d0a9 is history). |
| A-MINOR-7 / B-F1 | Dated DECISIONS 2026-09-10 entries appended (e423 ALLOCATED + review disposition; gate-(d) CONFIRMED amendment) — the amended test's and the artifact writers' "DECISIONS 2026-09-10" pointers now resolve. |
| B-F2 / B-F13 | DECISIONS supersession note records the frozen-text-vs-e423 conflict resolution + the delegation-scope interpretation.  Frozen text bytes UNTOUCHED. |
| B-F3 | This record committed with reviews + disposition; "wf_<pending>" placeholder replaced with the actual auditor identification. |
| B-F5 | Baseline-precise full-suite numbers stated (1382/2 at bdc1d0a9; post-disposition below). |
| B-F6 | Explicit "max_entries" checkpoint field + regression test; both dispatch docstrings' claim now literally true. |
| B-F9 | Accepted + documented: e420_census.py is bound by no registration; the scheduler repair makes any future re-run of e420/e421 fetch semantics strictly safer (404 exclusion semantics preserved and tested; all other non-2xx now typed instead of silently ok). |
| A/B PASS NOTES | No action — recorded as verification-positive. |

Post-disposition test state (this commit): test_e422_* suites 115
(dispatch 5→3 under the dated tombstone amendment; all other e422 files
unchanged); test_e423_* suites 54 (runner 48, dispatch 6); shared-infra
test_e420_census.py 20 — repair-relevant total 189.  All e420/e422/e423
files: 483 passed, 2 skipped (pre-existing e420 pool skips).  FULL REPO
SUITE: 1400 passed, 2 skipped, 0 failed.

## e423 record anchoring

This record is OTS-stamped (4 BTC calendars) after its commit;
STAMPED_PENDING until block inclusion, then upgraded + amended per the
e421/e422 pattern (DECISIONS 2026-09-10).
