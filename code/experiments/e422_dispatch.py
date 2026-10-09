"""SUPERSEDED — refuse-to-run tombstone (e423 dispatch-machinery repair,
DECISIONS 2026-09-10).

This CLI still wired the FROZEN e422_runner.py, whose census
request_options carry the invalid "scoring_strategy" field (RCSB v2
rejects it with HTTP 400 on every page) and whose POST lane marks
non-2xx responses ok — the shakedown-reproduced silent
empty-checkpoint defect.  Running it would manufacture a well-formed
EMPTY batch checkpoint: the worst possible failure mode for a
prospective registration.  main() therefore REFUSES unconditionally
(typed E423_DISPATCH_SUPERSEDED, exit 2); even --dry-run is refused,
because the query it prints is the invalid one.

Registered dispatch runs through experiments/e423_dispatch.py
(docs/e423_dispatch_repair_record_20260910.md).  The e422 FROZEN
registration text, T0, artifact schemas, and study identity are
UNCHANGED.  This file is NOT a freeze-bound object (informational
wiring in docs/e422_freeze_gate_a_20260909.md); its superseded
informational SHA row is recorded in the e423 repair record and
DECISIONS 2026-09-10.
"""
from __future__ import annotations

import sys

SUPERSEDED_CODE = "E423_DISPATCH_SUPERSEDED"


def main(argv: list[str] | None = None) -> int:
    """Refuse unconditionally — see the module docstring."""
    print(
        f"TYPED_FAILURE {SUPERSEDED_CODE}: experiments/e422_dispatch.py is "
        "SUPERSEDED (e423 dispatch-machinery repair; DECISIONS 2026-09-10). "
        "It wires the FROZEN e422_runner.py, whose census query carries the "
        "invalid scoring_strategy field (HTTP 400 on every page) and whose "
        "POST lane marks non-2xx responses ok (silent empty checkpoint). "
        "Dispatch through experiments/e423_dispatch.py instead.",
        file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
