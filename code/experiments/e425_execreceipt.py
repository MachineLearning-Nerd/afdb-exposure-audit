"""e425 execution receipts — registered wall-clock deviation reporting.

Allocated 2026-09-10 by w1:pG under main's standing delegation
(DECISIONS 2026-09-10T05:10+05:30), from the frozen-text conformance
sweep: §6 promises "any evaluation exceeding its budget is reported as
an execution deviation in the evaluation artifact, never silently
truncated", but the FROZEN implementation set carries no timing
whatsoever (verified: zero timing references in e422_evaldriver.py and
e422_protocol.py), and the artifact schema is self-digested — a timing
field cannot be added without touching frozen bytes (forbidden; any
change is a new E-number).  This module is the reporting vehicle.

Registered wall-clock budgets (§6): batch-4 <= 24 h, batch-8 <= 24 h,
batch-12 <= 72 h.  Parallelism is PERMITTED (registered 4-process
CQR lane), not required — a serial batch-8 (~56 CPU-h) or batch-12
(~127 CPU-h) run exceeds its budget and is REPORTED here, which the
registration treats as a compliant execution deviation, never a
violation (no repeat is dropped for budget reasons).

Layering: run_logged_evaluation = timing + receipt + e424 gate
(binding-state digest verification) + FROZEN driver.  The receipt is a
PER-RUN MEASUREMENT with real-clock timestamps — deliberately excluded
from the determinism property (the evaluation ARTIFACTS remain
bit-reproducible; the receipt measures this particular execution).
The registered batch-4 serial run (~15 CPU-h estimated) is within
budget; the measured value from its receipt informs the batch-8/12
serial-vs-parallel decision (registered decision point: no later than
the batch-7 close, 2026-10-30; DECISIONS 2026-09-10).
"""
from __future__ import annotations

import datetime
import os
import sys

import e420_census as census
import e422_evaldriver as driver
import e422_manifest as manifest_mod
import e424_evalgate as gate

RECEIPT_SCHEMA = "e422-execution-receipt-v1"
BUDGET_HOURS = {4: 24.0, 8: 24.0, 12: 72.0}


def run_logged_evaluation(
    t0_utc: str,
    batch_index: int,
    out_dir: str,
    result_dir: str,
    *,
    seeds=None,
    prior_binding_path: str | None = None,
    now_iso: str | None = None,
) -> dict:
    """Timed wrapper around the gated evaluation; writes the execution
    receipt next to the evaluation artifact.  Same signature and return
    as run_gated_evaluation / run_scheduled_evaluation.  Gate refusals
    precede the driver AND the receipt (a refused run produced no
    execution — there is nothing to report)."""
    started = datetime.datetime.now(datetime.timezone.utc)
    artifact = gate.run_gated_evaluation(
        t0_utc, batch_index, out_dir, result_dir, seeds=seeds,
        prior_binding_path=prior_binding_path, now_iso=now_iso)
    # budget lookup only AFTER the gated call: the gate (then the frozen
    # driver) types every refusal first, so an unscheduled batch is the
    # driver's typed E422_UNSCHEDULED_EVALUATION — never a receipt-layer
    # KeyError.  Post-gate, batch is in SCHEDULE == BUDGET_HOURS keys.
    budget_h = BUDGET_HOURS[batch_index]
    completed = datetime.datetime.now(datetime.timezone.utc)
    wall_seconds = (completed - started).total_seconds()
    within = wall_seconds <= budget_h * 3600.0
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "e_number": "e422",
        "t0_utc": t0_utc,
        "batch_index": batch_index,
        "role": driver.SCHEDULE[batch_index],
        "evaluated_utc": artifact["evaluated_utc"],
        "started_utc": started.isoformat(),
        "completed_utc": completed.isoformat(),
        "wall_seconds": wall_seconds,
        "registered_budget_hours": budget_h,
        "within_budget": within,
        "deviation_note": None if within else (
            f"registered wall-clock budget exceeded "
            f"({wall_seconds / 3600.0:.2f} h > {budget_h:g} h); the "
            "expected cause is serial execution without the registered "
            "4-process CQR lane (section 6 execution detail) — reported "
            "as an execution deviation, never truncated; no repeat is "
            "dropped for budget reasons"),
        "note": ("execution-deviation reporting vehicle per section 6; "
                 "the frozen evaluation artifact carries no timing field, "
                 "so the deviation is reported in the evaluation artifact "
                 "SET (evaluation artifact + this receipt, committed and "
                 "anchored together) per DECISIONS 2026-09-10"),
    }
    receipt["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in receipt.items() if k != "sha256"})
    census._atomic_write_json(
        os.path.join(result_dir,
                     f"e422_batch{batch_index:02d}_execution_receipt.json"),
        receipt)
    return artifact


REGISTERED_T0 = gate.REGISTERED_T0  # single source: the gate pins it


def main(argv: list[str] | None = None) -> int:
    """Registered scheduled-evaluation entry point (runbook §"Scheduled
    evaluations").  Typed-failure envelope: a GateError / DriverError /
    CensusError prints ``TYPED_FAILURE <code>`` on stderr and returns 2 —
    never a bare traceback, never a silent drop."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t0", default=REGISTERED_T0,
                        help="T0 instant (default: the frozen gate-(a) bind)")
    parser.add_argument("--batch", type=int, required=True,
                        help="scheduled batch (4, 8, or 12)")
    # e427 §4.2: defaults point at the CORRECTED ledger; the original
    # results/e422 ledger/evals are frozen (label defect) and refused.
    parser.add_argument("--out-dir", default="results/e427/batches",
                        help="cumulative-ledger directory (dispatch checkpoints)")
    parser.add_argument("--result-dir", default="results/e427/evals",
                        help="evaluation-artifact output directory")
    parser.add_argument("--prior-binding-path",
                        default="results/e427/batches/e422_binding_state.json",
                        help="binding-sidecar path (the chain state)")
    args = parser.parse_args(argv)
    # e427 program pause (human directive 2026-10-07): refuse before any
    # side effect.
    from e427_pause import CODE as PAUSE_CODE, pause_reason
    reason = pause_reason()
    if reason is not None:
        print(f"TYPED_FAILURE {PAUSE_CODE}: e422 program paused — {reason}; "
              "see results/e422/PAUSE.json", file=sys.stderr)
        return 2
    frozen = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "results", "e422")
    for path in (args.out_dir, args.result_dir, args.prior_binding_path):
        ap = os.path.abspath(path)
        if ap == frozen or ap.startswith(frozen + os.sep):
            print(f"TYPED_FAILURE E427_FROZEN_LEDGER_DIR: {path!r} is inside the "
                  "frozen original e422 ledger; evaluate results/e427 "
                  "(docs/e427_registration.md section 4.2)", file=sys.stderr)
            return 2
    try:
        run_logged_evaluation(
            args.t0, args.batch, args.out_dir, args.result_dir,
            prior_binding_path=args.prior_binding_path)
    except (gate.GateError, driver.DriverError, census.CensusError) as exc:
        # typed failure: code on stderr, nonzero exit — never silent
        print(f"TYPED_FAILURE {exc.code}: {exc.detail}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
