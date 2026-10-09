# e424 binding-state integrity gate — record (2026-09-10, w1:pG)

**Allocation.** e424 allocated 2026-09-10 by w1:pG under main's standing
delegation (DECISIONS 2026-09-09T20:31+05:30; successor-E-number
interpretation per DECISIONS 2026-09-10).  Trigger: the heartbeat-5
binding state machine audit (registered improvement work while gate (e)
waits on the 2026-09-18 batch-1 close).

## Finding

The registered evaluation chain digest-enforces every link EXCEPT one.
The FROZEN e422_evaldriver (bound content sha256 9307dd1a…) reads the
binding sidecar (`prior_binding_path`) with plain `json.loads`
(evaldriver lines ~309–311) and NEVER verifies the sidecar's own
"sha256" field, and tolerates a MISSING sidecar at batches 8/12 as an
empty prior (`prior = {}`).  Every other link is enforced: cumulative
checkpoints (E422_CHECKPOINT_DIGEST_MISSING/MISMATCH), the floor-decision
receipt (self-digest, asserted in tests), the evaluation artifacts
(self-digest).  Consequence: a tampered, stale, or wrong-path sidecar
would silently feed forged prior bindings into the day-56/day-84
artifacts — the silent-acceptance class e423 repaired on the census
path (its defect 2), here on the binding path.  The P6/P7a/P7b
conclusions flow through exactly this sidecar chain.

## Decision

Repair as a NEW E-number per the gate-(a) freeze discipline (any change
to a bound object is a new E-number): the frozen driver is UNTOUCHED;
experiments/e424_evalgate.py wraps it.

- batch < 8 (day-28): sidecar optional (the driver creates it); if a
  file already exists at the path it is digest-verified BEFORE the
  evaluation — unexpected pre-existing state is refused, never silently
  overwritten.
- batch >= 8 (day-56/day-84): sidecar REQUIRED — missing → typed
  E424_BINDING_STATE_MISSING (a chain break must not look like a fresh
  study); present → digest-verified (E424_BINDING_STATE_DIGEST_MISMATCH
  on failure); malformed/wrong-schema → E424_BINDING_STATE_SCHEMA.
- post-write: after the driver (re)writes the sidecar it is re-read and
  re-verified (E424_BINDING_STATE_WRITE_UNVERIFIED on failure).
- refusals precede every side effect (the driver is not called; tests
  assert this with a sentinel).
- accepted TOCTOU window between gate verification and the driver's own
  read: single-operator local machinery; the gate catches
  corruption/staleness, it is not an adversarial control (documented in
  the module docstring).

## Registration deltas

- The runbook's scheduled-evaluation invocation now routes through
  `e424_evalgate.run_gated_evaluation` (same signature/return as the
  frozen entry point).  All registered evaluations use the gate.
- No frozen object changes: the five bound e422 modules and the frozen
  registration text are byte-identical before/after (verified at
  commit).  e423 dispatch machinery untouched (its pre-close guard is
  complementary — dispatch-side and evaluation-side gates).

## Tests

tests/test_e424_evalgate.py — 11 passing (2026-09-10):
- verification function: absent+optional → None; absent+required →
  MISSING; valid → doc; tampered binding → DIGEST_MISMATCH; malformed
  JSON / scalar / wrong-schema (with recomputed valid digest!) →
  SCHEMA.
- gated end-to-end vs the REAL frozen driver on synthetic ledgers:
  batch-4 creates + post-write verifies; batch-4 refuses a pre-existing
  tampered sidecar (sentinel proves the driver is never reached);
  batch-8 requires the sidecar (absent → refusal, driver unreached);
  batch-8 refuses tampered prior state; healthy chain — AND the
  previously unpinned pin: P6/P7a/P7b bind at batch 8 on the decided
  six-bin grid (bound_at_batch == 8 whenever an outcome exists) with
  P5 immutability carried from batch 4.

## Anchoring

This record is OTS-stamped after its commit; STAMPED_PENDING until
block inclusion, then upgraded + amended per the e421/e422/e423
pattern.
