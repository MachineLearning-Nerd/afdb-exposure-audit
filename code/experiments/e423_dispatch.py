"""e423 gate-(e) dispatch CLI (dispatch-machinery repair): execute
registered batch k over the live lanes.

Thin wiring around experiments/e423_runner.run_batch — no registered rule
lives here.  Bounded dispatch: the registered census bound (CENSUS_MAX_PAGES
x CENSUS_PAGE_ROWS) applies inside the runner; --max-entries further bounds
a dispatch for shakedown runs (recorded in the checkpoint).  Every fetch is
receipted and its raw bytes persisted (the scheduler's discipline); every
failure is typed and exits nonzero — never silent.

Usage:
  .venv/bin/python experiments/e423_dispatch.py --t0 2026-09-11T00:00:00Z \
      --batch 1 --raw-dir results/e422/raw --out-dir results/e422/batches \
      [--manifest results/e422/manifest_checkpoint.json] [--max-entries 25]

Offline sanity mode (no network): --dry-run prints the batch window and the
exact census query and exits.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
from e423_runner import (  # noqa: E402
    CENSUS_PAGE_ROWS,
    build_census_query,
    batch_window_utc,
    run_batch,
)

import e427_label_repair as label_repair  # noqa: E402

DEFAULT_T0 = "2026-09-11T00:00:00Z"  # bound in the frozen registration (gate a)
# e427 (docs/e427_registration.md §4.2): every dispatch labels with the
# label-space repair; the id is recorded in each checkpoint.
E427_LABELLER_ID = "e427_label_repair (label_seq_id mapping + identity check)"

USER_AGENT = ("moluq-e423/1.0 (prospective registration; contact: repo owner)")


class PostTransport(census.RequestsTransport):
    """GET + receipted-JSON POST transport (the one allowlisted POST is the
    RCSB search API census query; the scheduler receipts every attempt)."""

    def post_json(self, url: str, payload: bytes, timeout_s: float):
        import requests
        r = requests.post(
            url, data=payload,
            headers={"User-Agent": self._user_agent,
                     "Content-Type": "application/json"},
            timeout=timeout_s)
        return census.HttpResponse(r.status_code, dict(r.headers), r.content)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t0", default=DEFAULT_T0,
                        help="T0 instant (default: the frozen gate-(a) bind)")
    parser.add_argument("--batch", type=int, required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest", default=None,
                        help="snapshot-manifest checkpoint path (admission)")
    parser.add_argument("--max-entries", type=int, default=None,
                        help="bound the dispatch to the first N census entries")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the window + census query; no network")
    args = parser.parse_args(argv)

    if args.t0 != DEFAULT_T0:
        # e423 T0 pin (detached review 2026-09-10, wf_d8bff144): t0 is
        # otherwise procedural-only — any t0 yields self-consistently
        # closed windows.  The registered entry point pins the frozen
        # gate-(a) bind; a different t0 is a typed refusal BEFORE any
        # window arithmetic, network, or filesystem side effect.
        print(f"TYPED_FAILURE E423_T0_MISMATCH: --t0 {args.t0!r} is not "
              f"the frozen registration bind {DEFAULT_T0!r} — the "
              "registered dispatch path pins the T0 identity "
              "(frozen text dc7a3627)", file=sys.stderr)
        return 2

    # e427 program pause (human directive 2026-10-07): refuse after the
    # side-effect-free T0 pin and before any fetch or write; --dry-run
    # stays available (it fetches and writes nothing).
    if not args.dry_run:
        from e427_pause import CODE as PAUSE_CODE, pause_reason
        reason = pause_reason()
        if reason is not None:
            print(f"TYPED_FAILURE {PAUSE_CODE}: e422 program paused — {reason}; "
                  "see results/e422/PAUSE.json", file=sys.stderr)
            return 2

    # e427 §4.2: results/e422/batches is the frozen ORIGINAL ledger (label
    # defect); the corrected cumulative ledger is results/e427/batches.
    frozen_dir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "results", "e422", "batches")
    if not args.dry_run and os.path.abspath(args.out_dir) == frozen_dir:
        print("TYPED_FAILURE E427_FROZEN_LEDGER_DIR: --out-dir "
              f"{args.out_dir!r} is the frozen original e422 ledger; dispatch "
              "to results/e427/batches (docs/e427_registration.md section 4.2)",
              file=sys.stderr)
        return 2

    start_iso, close_iso = batch_window_utc(args.t0, args.batch)
    query = build_census_query(start_iso, close_iso)
    print(f"batch {args.batch}: window ({start_iso}, {close_iso}]")
    print(f"census query: {json.dumps(query, sort_keys=True)}")
    if args.dry_run:
        print("dry-run: no fetch performed")
        return 0

    try:
        checkpoint = run_batch(
            args.t0, args.batch,
            transport=PostTransport(user_agent=USER_AGENT),
            raw_dir=args.raw_dir,
            out_dir=args.out_dir,
            manifest_path=args.manifest,
            max_entries=args.max_entries,
            parse_mappings_fn=label_repair.parse_mappings_labelspace,
            label_triple_fn=label_repair.label_triple_labelspace,
            labeller_id=E427_LABELLER_ID,
        )
    except (census.CensusError, manifest_mod.ManifestError) as exc:
        # typed failure: code on stderr, nonzero exit — never silent
        # (ManifestError widened in 2026-09-10: a tampered manifest
        # checkpoint refuses typed in the runner's pre-network phase)
        print(f"TYPED_FAILURE {exc.code}: {exc.detail}", file=sys.stderr)
        return 2
    print(f"checkpoint: n_rows={checkpoint.get('n_rows')} "
          f"census_entries={checkpoint.get('census_entries')} "
          f"sha256={checkpoint.get('sha256', '')[:16]}…")
    exclusions = checkpoint.get("exclusions", {})
    print("typed exclusions: "
          + json.dumps({k: len(v) if isinstance(v, list) else v
                        for k, v in exclusions.items()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001 - typed-by-class, never silent
        # e423 hardening (detached review MINOR-5, 2026-09-10): a residual
        # non-typed exception still exits 2 with its class on stderr —
        # loud, receipted, never a bare traceback
        print(f"TYPED_FAILURE E423_UNEXPECTED_{type(exc).__name__}: {exc}",
              file=sys.stderr)
        raise SystemExit(2)
