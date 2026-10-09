"""e427 §4.1 replay relabel of e422 batches 1-4 (no network).

Re-executes the registered e423_runner.run_batch for batches 2-4 with a
replay transport that serves every request from results/e422/raw by its
receipted URL (attempt order preserved) and verifies each raw body against
its SHA-256 name.  Two passes per batch, each on a fresh manifest chained
batch to batch:

  SELF-CHECK  frozen labeller -> must reproduce the original checkpoint's
              rows, exclusions, counts, and manifest_checkpoint_sha256.
  CORRECTED   e427 labeller -> written to --out-dir as
              e422_batchKK_checkpoint.json (registered schema; receipts =
              original fetch receipts; e427_relabel provenance block;
              sha256 recomputed by the registered canonical rule).

Batch 1 (attested empty, e426) is copied byte-identical.  Any self-check
failure exits 2 (TYPED_FAILURE E427_REPLAY_SELF_CHECK_FAILED) before any
corrected checkpoint is written.

Usage: .venv/bin/python experiments/e427_replay_relabel.py \
           --out-dir results/e427/batches --report results/e427/replay_report.json
"""
from __future__ import annotations

import argparse
import collections
import datetime
import hashlib
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import e420_census as census  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e423_runner as runner  # noqa: E402
import e427_label_repair as repair  # noqa: E402

ORIG_DIR = os.path.join(ROOT, "results/e422/batches")
RAW_DIR = os.path.join(ROOT, "results/e422/raw")
T0 = runner.REGISTERED_T0
LABELLER_ID = "e427_label_repair (label_seq_id mapping + identity check)"


class ReplayError(Exception):
    pass


class ReplayTransport:
    """Serves receipted responses in attempt order, keyed by (method, url)."""

    def __init__(self, receipts, raw_dir):
        self.queue = collections.defaultdict(collections.deque)
        for r in receipts:
            self.queue[(r.get("method", "GET"), r["url"])].append(r)
        self.raw_dir = raw_dir
        self.served = 0

    def _serve(self, method, url):
        q = self.queue.get((method, url))
        if not q:
            raise ReplayError(f"no receipted response for {method} {url}")
        r = q.popleft()
        self.served += 1
        if r.get("status") is None:
            raise ConnectionError(r.get("error") or "receipted transport error")
        content = b""
        if r.get("raw_name"):
            with open(os.path.join(self.raw_dir, r["raw_name"]), "rb") as fh:
                content = fh.read()
            if hashlib.sha256(content).hexdigest() != r["raw_name"]:
                raise ReplayError(f"raw cache hash mismatch {r['raw_name']}")
        return census.HttpResponse(int(r["status"]), dict(r.get("headers") or {}), content)

    def get(self, url, timeout_s):
        return self._serve("GET", url)

    def post_json(self, url, payload, timeout_s):
        return self._serve("POST", url)

    def leftover(self):
        return sum(len(q) for q in self.queue.values())


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def replay(batch, orig, manifest_path, tmp, parse_fn=None, label_fn=None, labeller_id=None):
    out = tempfile.mkdtemp(dir=tmp)
    transport = ReplayTransport(orig["receipts"], RAW_DIR)
    ck = runner.run_batch(T0, batch, transport=transport, raw_dir=tempfile.mkdtemp(dir=tmp),
                          out_dir=out, manifest_path=manifest_path,
                          max_entries=orig.get("max_entries"), sleeper=lambda s: None,
                          now_iso=orig["completed_utc"], parse_mappings_fn=parse_fn,
                          label_triple_fn=label_fn, labeller_id=labeller_id)
    return ck, transport


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args(argv)
    if os.path.abspath(args.out_dir) == os.path.abspath(ORIG_DIR):
        print("TYPED_FAILURE E427_FROZEN_LEDGER_DIR: refusing to write into the "
              "original e422 ledger", file=sys.stderr)
        return 2
    tmp = tempfile.mkdtemp(prefix="e427_replay_")
    report = {"schema": "e427-replay-report-v1", "t0_utc": T0, "batches": {},
              "replay_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    corrected = {}
    man_old = os.path.join(tmp, "manifest_old.json")
    man_new = os.path.join(tmp, "manifest_new.json")
    for k in (2, 3, 4):
        orig = load(os.path.join(ORIG_DIR, f"e422_batch{k:02d}_checkpoint.json"))
        # --- self-check with the frozen labeller -------------------------
        ck_old, tr_old = replay(k, orig, man_old, tmp)
        checks = {
            "rows_equal": ck_old["rows"] == orig["rows"],
            "exclusions_equal": ck_old["exclusions"] == orig["exclusions"],
            "census_entries_equal": ck_old["census_entries"] == orig["census_entries"],
            "entries_with_mappings_equal": ck_old["entries_with_mappings"] == orig["entries_with_mappings"],
            "manifest_admitted_new_equal": ck_old["manifest_admitted_new"] == orig["manifest_admitted_new"],
            "verified_accessions_equal": ck_old["verified_accessions"] == orig["verified_accessions"],
            "manifest_checkpoint_sha256_equal":
                ck_old["manifest_checkpoint_sha256"] == orig["manifest_checkpoint_sha256"],
            "all_receipts_consumed": tr_old.leftover() == 0,
        }
        # --- corrected pass with the e427 labeller -----------------------
        ck_new, tr_new = replay(k, orig, man_new, tmp, repair.parse_mappings_labelspace,
                                repair.label_triple_labelspace, LABELLER_ID)
        n_fail = collections.Counter(f["code"] for f in ck_new["exclusions"]["triple_failures"])
        o_fail = collections.Counter(f["code"] for f in orig["exclusions"]["triple_failures"])
        report["batches"][str(k)] = {
            "original_sha256": orig["sha256"], "self_check": checks,
            "original_n_rows": orig["n_rows"], "corrected_n_rows": ck_new["n_rows"],
            "original_triple_failures": dict(o_fail), "corrected_triple_failures": dict(n_fail),
            "corrected_mapping_failures": len(ck_new["exclusions"]["mapping_failures"]),
            "original_mapping_failures": len(orig["exclusions"]["mapping_failures"]),
            "corrected_manifest_sha_equal": ck_new["manifest_checkpoint_sha256"] == orig["manifest_checkpoint_sha256"],
            "receipts_unconsumed_corrected": tr_new.leftover(),
        }
        if not all(checks.values()):
            report["verdict"] = f"SELF_CHECK_FAILED_BATCH_{k}"
            json.dump(report, open(args.report, "w"), indent=1, sort_keys=True)
            print(f"TYPED_FAILURE E427_REPLAY_SELF_CHECK_FAILED: batch {k}: "
                  f"{[c for c, v in checks.items() if not v]}", file=sys.stderr)
            return 2
        # registered schema + provenance; original fetch receipts
        ck_new["receipts"] = orig["receipts"]
        ck_new["completed_utc"] = orig["completed_utc"]
        ck_new["freeze_receipt_sha256"] = orig["freeze_receipt_sha256"]
        ck_new["e427_relabel"] = {
            "registration": "docs/e427_registration.md section 4.1",
            "original_checkpoint_sha256": orig["sha256"],
            "labeller": LABELLER_ID,
            "replay_utc": report["replay_utc"],
            "source": "results/e422/raw (hash-verified replay of the original receipts)",
        }
        ck_new.pop("sha256", None)
        ck_new["sha256"] = manifest_mod.canonical_sha256(ck_new)
        corrected[k] = ck_new
    os.makedirs(args.out_dir, exist_ok=True)
    shutil.copy2(os.path.join(ORIG_DIR, "e422_batch01_checkpoint.json"),
                 os.path.join(args.out_dir, "e422_batch01_checkpoint.json"))
    for k in (1, 2, 3, 4):
        src = os.path.join(ORIG_DIR, f"e422_batch{k:02d}_freeze.json")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(args.out_dir, f"e422_batch{k:02d}_freeze.json"))
    for k, ck in corrected.items():
        census._atomic_write_json(os.path.join(args.out_dir, f"e422_batch{k:02d}_checkpoint.json"), ck)
        report["batches"][str(k)]["corrected_sha256"] = ck["sha256"]
    report["verdict"] = "OK"
    json.dump(report, open(args.report, "w"), indent=1, sort_keys=True)
    shutil.rmtree(tmp, ignore_errors=True)
    print(json.dumps({k: {kk: v[kk] for kk in ("original_n_rows", "corrected_n_rows")}
                      for k, v in report["batches"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
