# e422 prospective dispatch runbook (§6 cadence, 2026-09-10)

Operational procedure for the 12-week registered cadence after the
gate-(a) freeze (T0 = 2026-09-11T00:00:00Z).  This runbook is
INFORMATIONAL — the binding rules live in the FROZEN registration
(git-blob 6e1b9ce1…, content sha256 dc7a3627…) and the bound modules;
any conflict is resolved by the frozen text, never by this file.
Dispatch machinery is e423 (docs/e423_dispatch_repair_record_20260910.md);
changes to any bound object are a NEW E-number.

## Registered schedule (all UTC)

ACCELERATION AMENDMENT (DECISIONS 2026-09-10, user-approved): registered
batch 1 is a 48-HOUR window, (2026-09-11T00:00:00Z, 2026-09-13T00:00:00Z]
— dispatch closes 5 days earlier than the frozen weekly cadence; batches
2–12 remain 7-day and shift by batch 1's early close (the schedule stays
gap-free: start(k) = close(k-1)).  The frozen registration text is never
edited; this table mirrors the AMENDED arithmetic implemented in
e423_runner.batch_window_utc (registered-T0-only — dev/historical T0s
keep the pure 7-day windows).  Nominal day labels are driver ROLES keyed
to batch index: the "day-28/56/84" evaluations land at actual days
23/51/79 post-T0.

| batch | window (start, close] | role |
|---|---|---|
| 1 | 2026-09-11T00:00:00Z → 2026-09-13T00:00:00Z (48 h) | dispatch |
| 2 | 2026-09-13 → 2026-09-20 | dispatch |
| 3 | 2026-09-20 → 2026-09-27 | dispatch |
| 4 | 2026-09-27 → 2026-10-04 | dispatch + **batch-4 evaluation (nominal day-28; actual day 23)**: P5-90/P5-95 bind; per-bin tables + §3B joint grid + P6/P7 SUPPRESSED |
| 5 | 2026-10-04 → 2026-10-11 | dispatch |
| 6 | 2026-10-11 → 2026-10-18 | dispatch |
| 7 | 2026-10-18 → 2026-10-25 | dispatch |
| 8 | 2026-10-25 → 2026-11-01 | dispatch + **batch-8 evaluation (nominal day-56; actual day 51)**: floor decision receipt BEFORE any coverage table; P6/P7a/P7b bind on the decided grid |
| 9 | 2026-11-01 → 2026-11-08 | dispatch |
| 10 | 2026-11-08 → 2026-11-15 | dispatch |
| 11 | 2026-11-15 → 2026-11-22 | dispatch |
| 12 | 2026-11-22 → 2026-11-29 | dispatch + **close (nominal day-84; actual day 79)**: unbound predictions → DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE |

## Per-batch dispatch (every batch k, AFTER its window closes)

Run FROM THE REPO ROOT (the registered relative paths below resolve
from there).  MECHANICALLY ENFORCED since 2026-09-10 (detached review
wf_d8bff144): `--t0` is PINNED — any value other than the frozen
gate-(a) bind refuses with typed `E423_T0_MISMATCH` before window
arithmetic, network, or filesystem side effects (any t0 yields
self-consistently closed windows, so the pin is load-bearing).

```
.venv/bin/python experiments/e423_dispatch.py \
    --t0 2026-09-11T00:00:00Z --batch k \
    --raw-dir results/e422/raw --out-dir results/e422/batches \
    --manifest results/e422/manifest_checkpoint.json
```

- The census query covers exactly the registered window; the manifest
  accumulates incrementally (first appearance per accession, §2);
  the freeze receipt precedes all label construction (§7, runner-enforced).
- `--max-entries` is for SHAKEDOWNS ONLY (dev-only scratch dirs, never
  results/e422/).  A registered dispatch never bounds the census.
- No pre-T0 or pre-close registered dispatch: the query over an
  unclosed window is not the registered population.  MECHANICALLY
  ENFORCED since 2026-09-10 (DECISIONS): the e423 runner refuses with
  typed E423_WINDOW_NOT_CLOSED before any network or filesystem side
  effect; `now == close` is the allowed boundary (the (start, close]
  interval is complete).  Dry-run stays informational and prints the
  window so an operator can see the close instant before dispatch day.
- One dispatch per batch, MECHANICALLY ENFORCED since 2026-09-10
  (DECISIONS): an existing e422_batchNN_checkpoint.json refuses a
  second run with typed E423_CHECKPOINT_EXISTS before any side effect —
  a checkpoint is an anchored artifact and is never silently
  overwritten.  A documented re-dispatch is a dated DECISIONS entry
  plus MANUAL removal of the prior checkpoint, then the run.
- Typed failure (exit 2, `TYPED_FAILURE <code>` on stderr): LEDGER IT,
  never retry-by-force, never widen the allowlist silently.  A typed
  failure is a receipt, not an embarrassment.  Codes possible on
  dispatch day: `E423_WINDOW_NOT_CLOSED` (a MISSED batch dispatch is
  legal — the guard only refuses UNCLOSED windows; dispatch whenever
  the window is found closed and ledger the gap), `E423_CHECKPOINT_EXISTS`
  (one dispatch per batch — re-dispatch is a dated DECISIONS entry +
  MANUAL removal), `E423_T0_MISMATCH` (only if `--t0` was edited —
  run the command verbatim), `E423_CENSUS_DUPLICATE_ENTRIES` /
  `E422_CENSUS_*` (malformed census — ledger, do not dedupe silently),
  `E422_MANIFEST_*` (chain integrity — ledger, investigate before any
  re-run).

## Per-batch anchor + ledger (immediately after each dispatch)

1. OTS-stamp the checkpoint (+ manifest checkpoint if it changed):
   `.venv/bin/ots stamp results/e422/batches/e422_batchNN_checkpoint.json`
   (`results/e422/…` paths force-added to git; raw/ fetch bytes stay
   untracked — receipt SHAs live inside the checkpoint).
   The force-add is LOAD-BEARING: `.gitignore` line 15 excludes
   `results/` wholesale, and git cannot re-include children of an
   excluded directory, so the `!results/e4*/` negations are dead —
   every registered artifact needs explicit `-f`.  That is the safe
   configuration: plain `git add` REFUSES loudly on an ignored path,
   so an untracked checkpoint can never slip through silently; do not
   "fix" the dead negations (that would also un-ignore the dev-only
   results/e423_dev_shakedown/).
2. Commit: checkpoint + .ots + manifest checkpoint (explicit paths).
3. Ledger append (DECISIONS + docs/agent_status/lab.md): batch index,
   window, n_rows, typed-exclusion tally, checkpoint payload sha256
   (the checkpoint's own "sha256" field), commit id.
4. Upgrade the .ots after Bitcoin confirmation (`ots upgrade`) and
   amend by dated DECISIONS entry (e421 ce43caa2 pattern).

## Scheduled evaluations (batches 4, 8, 12 only)

Run FROM THE REPO ROOT, via the registered e425 CLI (T0-pinned,
double-evaluation-guarded, typed envelope; DECISIONS 2026-09-10):

```
.venv/bin/python experiments/e425_execreceipt.py --batch 4
.venv/bin/python experiments/e425_execreceipt.py --batch 8
.venv/bin/python experiments/e425_execreceipt.py --batch 12
```

(registered defaults: `--t0 2026-09-11T00:00:00Z`,
`--out-dir results/e422/batches`, `--result-dir results/e422/evals`,
`--prior-binding-path results/e422/batches/e422_binding_state.json` —
pass nothing else.)

The e425 execution-receipt layer (DECISIONS 2026-09-10) wraps the e424
gate (docs/e424_binding_gate_record_20260910.md) which wraps the FROZEN
driver: the binding sidecar is digest-verified before the evaluation
(the sidecar path is REQUIRED at every batch — at batches >= 8 a
missing FILE is a chain break, never a fresh study; at batch 4 the
supplied path is what makes the driver CREATE the sidecar), re-verified
after, and the evaluation's wall-clock duration is receipted against
the registered §6 budgets (24/24/72 h) — an exceedance is an
execution-deviation report, never a truncation.  Use this wrapper, not
the raw driver or the bare gate, for all registered evaluations.
Anchor + commit the execution receipt with the evaluation artifact.

Evaluation-day typed codes and remedies (all exit 2,
`TYPED_FAILURE <code>` on stderr; all LEDGER-IT, never retry-by-force):

- `E424_T0_MISMATCH`: the command was edited — run it verbatim.
- `E424_EVALUATION_EXISTS`: an evaluation artifact, floor-decision
  receipt, or execution receipt for this batch already exists — an
  evaluation is never silently re-run; re-run = dated DECISIONS entry
  + MANUAL removal of the prior outputs, then the run.
- `E424_BINDING_STATE_MISSING`: a chain break — the sidecar file is
  gone (batches >= 8) or no `--prior-binding-path` was supplied.  A
  missing file is restored from the last COMMITTED + anchored sidecar
  (see below), never papered over by letting the study start fresh.
- `E424_BINDING_STATE_DIGEST_MISMATCH` / `_SCHEMA` /
  `_CHAIN_POSITION` / `_WRITE_UNVERIFIED`: the sidecar is tampered,
  stale (wrong chain position), carries unregistered keys, or failed
  post-write verification — ledger it and STOP; do not proceed to the
  driver.
- `E422_CHECKPOINT_*` / `E422_UNSCHEDULED_EVALUATION`: frozen-driver
  ledger integrity / batch typo — ledger, fix, re-run (a refused run
  wrote nothing, so a re-run after a fix is safe and legal).

- seeds default = range(200) (registered R=200); do not override.
- batch 4: artifact carries P5-90/P5-95 + descriptive floor table ONLY
  (registered suppression; driver-enforced).
- batch 8: the floor-decision receipt is written BEFORE any coverage
  evaluation exists (decided-before-unblinding fence, §4); on a
  merged/restricted decision the conditional evaluation runs on the
  decided grid via the registered use_grid path.
- batch 12: pass the batch-8 binding state (chained through
  prior_binding_path sidecars); unbound predictions close as
  DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE (a registered outcome).
- Binding-sidecar semantics (driver-enforced, smoke-verified
  2026-09-10): the binding state is read + written AT the
  prior_binding_path — the batch-4 invocation CREATES it
  (e422-binding-state-v1); batches 8 and 12 read + rewrite the same
  path (the chain).  Never point two scheduled evaluations at
  different sidecar paths.
- Harness wiring was smoke-validated end-to-end against a synthetic
  12-batch ledger under results/e423_dev_shakedown/eval_smoke/
  (DEV-ONLY machinery validation — typed unscheduled/missing-batch
  failures, batch-4 suppression, batch-8 floor-decision fence +
  merged-grid conditional, batch-12 close; its numbers are dev
  numbers, never results).
- Multi-batch ACCUMULATION was dispatch-validated on a real second
  fetch chain (dev batch-2, historical window (2026-09-04,
  2026-09-11], --max-entries 25, 2026-09-10): manifest grew 12→13 with
  0 rewrites (first-appearance, §2); checkpoint→freeze and
  checkpoint→manifest linkage digests hold; freeze strictly precedes
  labels; the 1..2 chain loads under eval-side digest enforcement with
  a zero-row batch-2 (a validated quiet window chains as a legitimate
  registered outcome).
- Evaluation DETERMINISM was verified on the same dev ledger
  (2026-09-10): two independent 4→8→12 chains with pinned
  (ledger, seeds, now_iso, sidecar-input) produced BYTE-IDENTICAL
  artifacts — sidecar, all three evaluation artifacts, and both
  floor-decision receipts (six files, pairwise sha256 equal).
  (Dispatch shakedowns must use strictly-past
  windows — the runner's pre-close guard refuses anything else.)
  RE-VERIFIED 2026-09-10 through the FULL hardened registered stack
  (e425 → e424 with the detached-review guards → frozen driver;
  dev-only results/e423_dev_shakedown/e425_receipt_smoke/
  determinism_stack.py): 6/6 byte-identical.  The exact reproducibility
  contract: an evaluation's bytes are a function of (the COMMITTED
  cumulative checkpoints 1..k, seeds = range(200), the evaluated_utc
  recorded in the COMMITTED evaluation artifact, and the
  prior_binding_path STRING).  The frozen driver embeds the
  prior_binding_path string in the artifact, so run the registered
  command verbatim FROM THE REPO ROOT (relative sidecar path) — a
  reproduction from a different working directory lands different
  bytes with identical content.
- Anchor + commit the evaluation artifact and the floor-decision
  receipt exactly as in the per-batch procedure.  ALSO anchor + commit
  the binding sidecar `results/e422/batches/e422_binding_state.json`
  with EVERY scheduled evaluation (`git add -f` + `.ots` stamp + the
  commit): the sidecar IS the chain state, and the E424_BINDING_STATE_MISSING
  remedy ("restore the last committed sidecar") only works if a
  committed copy exists at every link.

## Integrity invariants (check every step)

- Checkpoint digests re-verify (`load_cumulative_ledger` enforces this
  at every evaluation — tamper = typed E422_CHECKPOINT_DIGEST_MISMATCH).
- The AFDB v6 manifest is immutable once written: admission only ADDS
  records, never rewrites (first-appearance rule); the temporal check
  is fail-closed (PRE_SNAPSHOT_RELEASE / ABSENT_FROM_SNAPSHOT).
- results/e423_dev_shakedown/ and any future shakedown dirs are
  DEV-ONLY: gitignored, never cited as results, never the dispatch
  target.

## Standing rules (unchanged)

Only explicit git paths; Co-Authored-By trailers; never force-push;
never touch lead.md / handoff.md; typed failures only, no silent
drops; exploratory ledger is development-only — NO headline numbers
from it; amendments via dated DECISIONS entries.

## DATED AMENDMENT 2026-09-18: batch 1 dispatched ATTESTED EMPTY via e426

Batch 1 was dispatched 2026-09-18 (commit 440566b4) as a QUIET WINDOW —
census_entries=0, n_rows=0, checkpoint sha256 4f54057d…, manifest
6508e0b2… — through experiments/e426_attested_empty.py (user-approved
same day).  During the RCSB serving incident both registered dispatch
attempts (2026-09-13) typed-failed E422_CENSUS_QUERY_FAILED on HTTP 204;
after index recovery the registered census filter ATTESTED the window
empty on a healthy lane (historical week == 140; 09-16 release event ==
156 filtered; batch-1 window 204).  e426 re-attests those three probes
live and fail-closed on EVERY run, then hands the attested-empty census
to the FROZEN e423 run_batch via its transport injection point — the
dev-verified quiet-window path.  Frozen bytes unchanged; the receipt
results/e422/e426_batch1_empty_attestation.json (1aaacdc7…) carries the
probe evidence and the synthesized census body verbatim.

Standing rule for FUTURE quiet windows: e426 is scoped to batch 1 ONLY.
Any later batch whose window attests empty under a healthy lane gets
its own dated E-numbered amendment — never a silent zero-row write,
never a widened query.

## DATED AMENDMENT 2026-10-07: dispatch schedule PRE-ARMED (canonical cron texts committed)

The weekly self-chaining design (each fire CronCreates only the NEXT
Sunday's one-shot) was replaced by a fully PRE-ARMED schedule at user
direction, after the 2026-10-07 context compaction dropped the in-memory
cron store — CronList returned EMPTY (all armed jobs gone: the
self-chaining design's single point of failure, realized).

One-shots now exist for EVERY remaining fire (schedule of record +
verbatim texts in docs/e422_cron_prompts_20261007.md; job IDs
session-local):

- Sat 2026-10-10 09:07 IST — batch-5 readiness probe (informational only)
- Sun 07:13 IST — 10-11, 10-18, 10-25, 11-08, 11-15, 11-22 (batches 5,6,7,9,10,11)
- Sun 10:07 IST — 11-01 (batch-8 dispatch + evaluation), 11-29 (batch-12 + close)

Weekly-chain text changes (canonical in that file): step (6) "RESCHEDULE
next Sunday" became "SCHEDULE SELF-HEAL: CronList; recreate any missing
one-shot from the canonical texts; never duplicate"; step (1) gained an
eval-failure exception — if k in {8,12} is closed+undispatched AFTER its
eval one-shot's fire time +30 min, the weekly chain dispatches it as a
legal late dispatch (missed-batch clause) but runs NO evaluation machinery;
step (5) gained a DECISIONS-commit caveat (never sweep a concurrent
main-phase edit of docs/DECISIONS.md into an e422 commit — commit lab.md
only and let DECISIONS.md ride).

A rebuild gap never loses registered time: windows are calendar-defined
and late dispatch is registry-neutral (missed-batch clause, §per-batch).
Reversal: CronDelete the weekly one-shots; restore single-link chaining
with the pre-2026-10-07 step (6) (preserved in this date's DECISIONS
entry and session transcript).  Entry written 2026-10-07T17:15:23+05:30.
