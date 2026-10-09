# e425 execution receipts — record (2026-09-10, w1:pG)

**Allocation.** e425 allocated 2026-09-10 by w1:pG under main's standing
delegation, from the frozen-text conformance sweep (heartbeat-9
registered improvement work; the first end-to-end read of the FROZEN
registration against the operational machinery).

## Finding (conformance sweep)

§6 promises: "Any evaluation exceeding its budget is reported as an
execution deviation in the evaluation artifact, never silently
truncated (no repeat is dropped for budget reasons)."  The FROZEN
implementation set carries NO timing machinery — verified: zero
timing/duration/budget references in experiments/e422_evaldriver.py
and experiments/e422_protocol.py — and the evaluation artifact is
self-digested with a frozen schema, so a timing field cannot be added
without touching frozen bytes (forbidden by the gate-(a) freeze
discipline; any change is a new E-number).

Related registered facts verified in the same sweep:
- Parallelism is PERMITTED, not required ("CQR fits MAY run in
  parallel … registered 4-process lane"); the frozen protocol's
  `run_evaluation` is serial.
- Estimated serial durations: batch-4 ≈ 15 CPU-h (within the ≤24 h
  budget even serially); batch-8 ≈ 56 CPU-h and batch-12 ≈ 127 CPU-h
  exceed serially — both are compliant IF reported (budget exceedance
  is an execution deviation, never a violation).
- Everything else swept conforms: T0/batch windows, manifest mechanism,
  batch-4 suppression (`_scope_batch4`), binding rules (pinned
  c2ec2b43), §7 decision-rule constants (e.g. `MARGIN_M = 0.01`),
  batch-12 AT_CLOSE, §8 disjoint-streams definition (e421-P HOLD), §9
  fences.

## Decision

**(1) Interpretation (pre-dispatch DECISIONS amendment).** The frozen
artifact schema cannot carry a timing field; the §6 deviation-report
promise is discharged by the evaluation artifact SET — the evaluation
artifact plus a per-run execution receipt, committed and anchored
together.  The receipt's own note field records this interpretation.

**(2) e425 implementation.** experiments/e425_execreceipt.py:
`run_logged_evaluation` = wall-clock timing + receipt + e424 binding
gate + FROZEN driver.  Registered budgets {4: 24 h, 8: 24 h, 12: 72 h}
(asserted against the registration by test); receipt carries
started/completed (real clock — a per-run measurement, deliberately
outside the determinism property), wall_seconds, within_budget, and a
typed deviation note on exceedance (reporting the expected serial
cause; "never truncated; no repeat is dropped").  Gate refusals
precede the driver AND the receipt — a refused run produced no
execution.  The runbook's scheduled-evaluation command now routes
through the e425 wrapper; the receipt is anchored + committed with the
evaluation artifact.

**(3) Serial-first, registered decision point.** Batch-4 (2026-10-09)
runs serially — within budget by the registered arithmetic.  The
serial-vs-registered-4-process-lane decision for batches 8/12 is
deferred to a REGISTERED DECISION POINT: no later than the batch-7
close (2026-10-30), informed by batch-4's MEASURED serial duration
(e425 receipt) against the §6 estimates.  Both paths are compliant
(parallel permitted by the §6 execution detail; serial exceedance
reportable); building a parallel executor before the measurement
exists would add divergence risk for time that is not yet known to be
needed.  A parallel executor, if chosen, is a new E-number with
byte-identity verification against the serial path (the §6 execution
detail's own criterion; the determinism harness provides the method).

## Tests

tests/test_e425_execreceipt.py — 4 passing (2026-09-10): receipt
written beside the artifact with digest verifying, budgets matching the
registration, deviation note on forced exceedance (monkeypatched
budget), and NO receipt when the e424 gate refuses (driver never
reached, sentinel).

## Anchoring

This record is OTS-stamped after its commit; STAMPED_PENDING until
block inclusion, then upgraded + amended per the e421/e422/e423/e424
pattern.
