"""e426 — ATTESTED-EMPTY-WINDOW dispatch (registered batch 1).

ALLOCATED 2026-09-18 (w1:pG under the standing delegation; DECISIONS
2026-09-18, "batch-1 deep-dive + attested-empty disposition"): the
registered batch-1 window (2026-09-11, 2026-09-13] is ATTESTED EMPTY
on a demonstrably healthy lane.  During the RCSB serving incident
(search index frozen 2026-09-10..2026-09-16+) the two registered
dispatch attempts (cron 2026-09-13 09:41 IST; OS failsafe 10:30 IST)
typed-failed E422_CENSUS_QUERY_FAILED on HTTP 204 — the deliberate
204-as-fail-closed design, which cannot distinguish "empty window"
from "broken lane".  After index recovery the registered-filter probe
suite observed (all via build_census_query, DECISIONS 2026-09-18):

    historical week (2026-09-03, 2026-09-10]   == 140  (stable lane)
    day (2026-09-10, 2026-09-11]               == 204 (empty)
    BATCH-1 (2026-09-11, 2026-09-13]           == 204 (empty)
    day (2026-09-13, 2026-09-14]               == 204 (empty)
    day (2026-09-14, 2026-09-15]               == 204 (empty)
    09-16 release event (2026-09-15, 2026-09-16] == 156  (345 raw)

i.e. the archive resumed publishing (PDBj: "345 new PDB entries have
been released on 2026-09-16"; 156 pass the registered census filter)
while batch-1's window is empty in the archive itself.

This module converts that state into a FIRST-CLASS registered
artifact: it (1) re-runs the attestation probes LIVE, fail-closed;
(2) only when all three pass, feeds the attested-empty census document
to the FROZEN e423 run_batch through its transport injection point,
which then produces the batch-1 checkpoint exactly as it would have
had the search API returned {"total_count": 0} — the dev-verified
quiet-window path (n_rows=0 chains cleanly, DECISIONS 2026-09-10).
NO frozen module bytes change; the attestation receipt records the
live probe evidence and the synthesized census body verbatim.

Typed failure codes (exit 2, TYPED_FAILURE <code> on stderr — receipt,
never retry-by-force):
    E426_CHECKPOINT_EXISTS        batch-1 checkpoint already anchored
    E426_LANE_HEALTH_FAILED       historical week != 140 (lane changed)
    E426_NO_POST_FREEZE_RELEASES  no post-freeze content (lane frozen)
    E426_WINDOW_NOT_EMPTY         batch-1 window has entries — the
                                  registered runner path applies, not
                                  this module
    E426_TRANSPORT_GET_UNEXPECTED the frozen code attempted a GET on a
                                  zero-entry batch — the attestation
                                  and the run disagree; loud refusal
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
from typing import Any

import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census  # noqa: E402
import e423_runner as runner  # noqa: E402

REGISTERED_T0 = runner.REGISTERED_T0
BATCH_INDEX = 1
SEARCH_URL = runner.SEARCH_URL
UA = {"User-Agent": "moluq-e426-attestation/1.0"}

# Registered paths (docs/e422_batch_runbook.md).
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.path.join(REPO, "results", "e422", "raw")
OUT_DIR = os.path.join(REPO, "results", "e422", "batches")
MANIFEST_PATH = os.path.join(REPO, "results", "e422",
                             "manifest_checkpoint.json")
ATTEST_PATH = os.path.join(REPO, "results", "e422",
                           "e426_batch1_empty_attestation.json")

# A1: the historical week's count under the registered census filter is
# the stable rehearsal figure (140, unchanged 2026-09-10..2026-09-18).
# An exact match is required: drift means the lane changed and the
# attestation is void (investigate; never silently pass).
HISTORICAL_WINDOW = ("2026-09-03T00:00:00Z", "2026-09-10T00:00:00Z")
HISTORICAL_EXPECTED = 140
# A2: post-freeze evidence — the 2026-09-16 release event, stamped
# exactly 2026-09-16T00:00:00Z (345 raw; 156 pass the registered
# filter).  Existence (>= 1) is required; the exact count is recorded.
EVENT_WINDOW = ("2026-09-15T00:00:00Z", "2026-09-16T00:00:00Z")

# The synthesized census body handed to the frozen runner: byte-identical
# to the canonical empty result document, recorded in the receipt.
ATTESTED_EMPTY_BODY = json.dumps(
    {"result_set": [], "total_count": 0}, sort_keys=True).encode("utf-8")


class E426Error(Exception):
    """Typed failure: exit 2, TYPED_FAILURE <code> on stderr."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


def _post_query(query: dict) -> tuple[int, bytes]:
    """The one live network call of the attestation layer (indirected for
    tests: monkeypatch this to inject probe responses)."""
    resp = requests.post(SEARCH_URL, json=query, headers=UA, timeout=30)
    return int(resp.status_code), resp.content


def _probe(start_iso: str, close_iso: str) -> dict:
    """One live registered-filter census probe; returns a receipt dict."""
    query = runner.build_census_query(start_iso, close_iso)
    probed_utc = census._utc_now()
    try:
        status, content = _post_query(query)
    except requests.RequestException as exc:
        raise E426Error(
            "E426_LANE_HEALTH_FAILED",
            f"probe ({start_iso}, {close_iso}] transport error: "
            f"{type(exc).__name__}: {exc}") from exc
    rec = {
        "window": {"start": start_iso, "close": close_iso},
        "http_status": int(status),
        "total_count": None,
        "query_sha256": census._sha256_hex(
            json.dumps(query, sort_keys=True).encode("utf-8")),
        "probed_utc": probed_utc,
    }
    if status == 200:
        try:
            doc = json.loads(content.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise E426Error(
                "E426_LANE_HEALTH_FAILED",
                f"probe ({start_iso}, {close_iso}]: undecodable body") from exc
        rec["total_count"] = doc.get("total_count")
    return rec


def attest_empty_batch1() -> list[dict]:
    """Run the three live probes fail-closed; returns receipt dicts.

    Raises E426Error (typed) unless batch-1 is attested empty on a
    healthy, recovered lane."""
    a1 = _probe(*HISTORICAL_WINDOW)
    if a1["total_count"] != HISTORICAL_EXPECTED:
        raise E426Error(
            "E426_LANE_HEALTH_FAILED",
            f"historical week ({HISTORICAL_WINDOW[0]}, "
            f"{HISTORICAL_WINDOW[1]}] returned "
            f"{a1['total_count']!r}, expected {HISTORICAL_EXPECTED} "
            "(the stable registered-filter count); lane changed — "
            "attestation void")
    a2 = _probe(*EVENT_WINDOW)
    if not (isinstance(a2["total_count"], int) and a2["total_count"] >= 1):
        raise E426Error(
            "E426_NO_POST_FREEZE_RELEASES",
            f"09-16 event window ({EVENT_WINDOW[0]}, {EVENT_WINDOW[1]}] "
            f"returned {a2['total_count']!r} — no post-freeze content; "
            "the lane is still frozen and emptiness is NOT attestable")
    start_iso, close_iso = runner.batch_window_utc(REGISTERED_T0, BATCH_INDEX)
    a3 = _probe(start_iso, close_iso)
    if not (a3["http_status"] == 204
            or a3["total_count"] == 0):
        raise E426Error(
            "E426_WINDOW_NOT_EMPTY",
            f"batch-{BATCH_INDEX} window ({start_iso}, {close_iso}] "
            f"returned http={a3['http_status']} "
            f"total_count={a3['total_count']!r} — the window has a "
            "population; dispatch via the registered e423 path, not "
            "this module")
    return [a1, a2, a3]


@dataclasses.dataclass(frozen=True)
class _Response:
    status: int
    headers: dict
    content: bytes


class AttestedEmptyTransport:
    """The transport handed to the frozen run_batch.

    post_json returns the attested-empty census document (HTTP 200) for
    EVERY POST — for a zero-entry batch the frozen census loop issues
    exactly one page-1 POST and stops; any additional call means the
    attestation and the run disagree, and the body shape itself will
    fail the frozen validation loudly (total_count 0 with a non-empty
    result_set cannot occur; a second page request cannot legitimately
    follow a total_count of 0).  get_json is never reachable on a
    zero-entry batch (no mappings, no metadata fetches) and refuses.
    """

    def __init__(self, body: bytes, attestation_sha256: str) -> None:
        self._body = body
        self._attestation_sha256 = attestation_sha256
        self.post_calls = 0

    def post_json(self, url: str, payload: bytes,
                  timeout_s: float) -> _Response:
        self.post_calls += 1
        return _Response(
            status=200,
            headers={"X-E426-Attested-Empty": self._attestation_sha256,
                     "Content-Type": "application/json"},
            content=self._body,
        )

    def get(self, url: str, timeout_s: float) -> _Response:
        raise E426Error(
            "E426_TRANSPORT_GET_UNEXPECTED",
            f"GET {url} on an attested-empty batch — the frozen run "
            "attempted a per-entry fetch, which contradicts the "
            "attestation; refusing loudly")


def _write_attestation_receipt(probes: list[dict],
                               run_pid: int) -> str:
    receipt = {
        "schema": "e426-empty-attestation-v1",
        "e_number": "e426",
        "registers": "batch 1 attested empty (DECISIONS 2026-09-18)",
        "t0_utc": REGISTERED_T0,
        "batch_index": BATCH_INDEX,
        "probes": probes,
        "synthesized_census_body_utf8": ATTESTED_EMPTY_BODY.decode("utf-8"),
        "synthesized_census_sha256": census._sha256_hex(ATTESTED_EMPTY_BODY),
        "frozen_runner": "e423 (run_batch; transport injection point)",
        "frozen_bytes_changed": False,
        "pid": run_pid,
    }
    receipt["sha256"] = census._sha256_hex(
        json.dumps(receipt, sort_keys=True).encode("utf-8"))
    census._atomic_write_json(ATTEST_PATH, receipt)
    return receipt["sha256"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="e426: attested-empty dispatch of registered batch 1")
    ap.add_argument("--execute", action="store_true",
                    help="perform the live attestation and dispatch "
                         "(default: print the plan and exit)")
    args = ap.parse_args(argv)

    checkpoint_path = os.path.join(
        OUT_DIR, f"e422_batch{BATCH_INDEX:02d}_checkpoint.json")
    if os.path.exists(checkpoint_path):
        raise census.CensusError(
            "E426_CHECKPOINT_EXISTS",
            f"{checkpoint_path} already exists; refusing to overwrite an "
            "anchored dispatch artifact")
    start_iso, close_iso = runner.batch_window_utc(REGISTERED_T0, BATCH_INDEX)
    if not args.execute:
        print(f"e426 plan: attest ({start_iso}, {close_iso}] empty on a "
              f"healthy lane (A1 historical == {HISTORICAL_EXPECTED}, "
              f"A2 09-16 event >= 1), then dispatch batch {BATCH_INDEX} "
              "via the frozen e423 run_batch with the attested-empty "
              "census.  Re-run with --execute.")
        return 0

    probes = attest_empty_batch1()
    attest_sha = _write_attestation_receipt(probes, os.getpid())
    print(f"e426: attestation PASSED and receipted -> {ATTEST_PATH} "
          f"(sha256 {attest_sha})")
    for rec in probes:
        print(f"  probe ({rec['window']['start']}, "
              f"{rec['window']['close']}]: http={rec['http_status']} "
              f"total_count={rec['total_count']}")

    transport = AttestedEmptyTransport(ATTESTED_EMPTY_BODY, attest_sha)
    checkpoint = runner.run_batch(
        REGISTERED_T0,
        BATCH_INDEX,
        transport=transport,
        raw_dir=RAW_DIR,
        out_dir=OUT_DIR,
        manifest_path=MANIFEST_PATH,
    )
    print(f"e426: frozen run_batch produced checkpoint "
          f"{checkpoint_path} (census_entries="
          f"{checkpoint.get('census_entries')}, n_rows="
          f"{checkpoint.get('n_rows')}, sha256={checkpoint.get('sha256')})")
    if transport.post_calls != 1:
        raise census.CensusError(
            "E426_TRANSPORT_GET_UNEXPECTED",
            f"expected exactly one census POST for a zero-entry batch, "
            f"saw {transport.post_calls}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except census.CensusError as exc:
        print(f"TYPED_FAILURE {exc.code}: {exc.detail}", file=sys.stderr)
        sys.exit(2)
    except E426Error as exc:
        print(f"TYPED_FAILURE {exc.code}: {exc.detail}", file=sys.stderr)
        sys.exit(2)
