"""Re-fetch the raw third-party files behind the released derived data and
check them against the SHA-256 digests recorded when they were first fetched.

Raw files (PDB/PDBe mmCIF, PDBe SIFTS mapping JSON, AFDB model files and
metadata, RCSB search responses, ATLAS archives) are not redistributed in
this repository. Every original request was receipted with its URL, method,
status and the SHA-256 of the response body. This script replays the GET
receipts of one receipt set, stores each body under the raw directory the
registered code expects (named by its SHA-256, as the original pipeline did),
and reports whether the bytes are identical to the recorded ones.

  MATCH    body SHA-256 equals the recorded digest -> identical bytes; saved
           as <raw_dir>/<sha256> so the registered replay code can use it.
  CHANGED  the source now serves different bytes (e.g. a PDB entry was
           re-versioned or AFDB re-issued a file); saved under
           <raw_dir>/_changed/<recorded>.<new> and NOT under the recorded
           name, so any replay that needs it fails loudly instead of silently
           using different inputs.
  FAILED   network / HTTP error.

POST receipts (RCSB Search API queries) are reported but not replayed: their
request bodies are rebuilt by the registered code, and a later search returns
a different result set because the archive grows. Use the registered dispatch
code with a replay transport for those, or treat the receipt digests as the
record.

Receipt sets (--set):
  e422         results/e422/batches/e422_batch0{1..4}_checkpoint.json receipts
               -> <work>/results/e422/raw/   (preregistered census, original bytes)
  e421relabel  results/e427/e421_relabel/receipts.json
               -> <work>/data/e427/e421_raw/ (post-2022 relabel, October 2026 bytes)
  e420census   results/e420/source_availability_manifest.json requests
               -> <work>/data/e420/census_raw/
  e420labels   results/e420/label_run_receipt.json requests
               -> <work>/data/e420/labels_raw/
  9qj6         the three files behind Figure 1 (9qj6 mmCIF, 9qj6 SIFTS mapping,
               AF-K7PQ54-F1-model_v6) from the e422 receipts -> <work>/results/e422/raw/

Usage:
  python -I data/fetch/fetch_raw.py --set 9qj6 --work _work
  python -I data/fetch/fetch_raw.py --set e422 --work _work [--limit 50] [--dry-run]
Network: yes (public endpoints, no credentials). Be polite: default 0.5 s
between requests, one request at a time.
"""
import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

UA = {"User-Agent": "afdb-exposure-audit-refetch/1 (research reproduction)"}
NINE_QJ6 = ("https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/9qj6",
            "https://www.ebi.ac.uk/pdbe/entry-files/9qj6.cif")


def receipts(derived, which):
    d = Path(derived)
    if which in ("e422", "9qj6"):
        out = []
        for b in (1, 2, 3, 4):
            out += json.loads((d / f"results/e422/batches/e422_batch0{b}_checkpoint.json").read_text())["receipts"]
        if which == "9qj6":
            out = [r for r in out if r.get("url") in NINE_QJ6
                   or ("AF-K7PQ54-F1" in r.get("url", "") and r["url"].endswith(".cif"))]
        return out, "results/e422/raw"
    if which == "e421relabel":
        return json.loads((d / "results/e427/e421_relabel/receipts.json").read_text()), "data/e427/e421_raw"
    if which == "e420census":
        return json.loads((d / "results/e420/source_availability_manifest.json").read_text())["requests"], "data/e420/census_raw"
    if which == "e420labels":
        return json.loads((d / "results/e420/label_run_receipt.json").read_text())["requests"], "data/e420/labels_raw"
    raise SystemExit(f"unknown set {which}")


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.status, r.read()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["e422", "e421relabel", "e420census", "e420labels", "9qj6"])
    ap.add_argument("--derived", default="data/derived")
    ap.add_argument("--work", default="_work", help="work root (see code/tools/make_workroot.py)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--delay", type=float, default=0.5)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    recs, sub = receipts(a.derived, a.set)
    raw = Path(a.work) / sub
    seen, todo, skipped_post = set(), [], 0
    for r in recs:
        digest = r.get("raw_name") or r.get("sha256")
        if r.get("status") != 200 or not digest:
            continue
        if r.get("method", "GET") != "GET":
            skipped_post += 1
            continue
        if (r["url"], digest) in seen:
            continue
        seen.add((r["url"], digest))
        todo.append((r["url"], digest))
    if a.limit:
        todo = todo[:a.limit]
    print(f"{a.set}: {len(todo)} GET bodies to check; {skipped_post} POST receipts not replayed; raw dir {raw}")
    if a.dry_run:
        for u, dg in todo[:20]:
            print(f"  {dg[:16]}  {u}")
        return 0
    raw.mkdir(parents=True, exist_ok=True)
    tally = {"MATCH": 0, "CHANGED": 0, "FAILED": 0, "PRESENT": 0}
    for u, dg in todo:
        if (raw / dg).exists():
            tally["PRESENT"] += 1
            continue
        try:
            status, body = fetch(u)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            tally["FAILED"] += 1
            print(f"FAILED   {u}: {exc}")
            time.sleep(a.delay)
            continue
        h = hashlib.sha256(body).hexdigest()
        if h == dg:
            (raw / dg).write_bytes(body)
            tally["MATCH"] += 1
        else:
            (raw / "_changed").mkdir(exist_ok=True)
            (raw / "_changed" / f"{dg}.{h}").write_bytes(body)
            tally["CHANGED"] += 1
            print(f"CHANGED  {u}: recorded {dg[:16]}..., now {h[:16]}...")
        time.sleep(a.delay)
    print(json.dumps(tally))
    return 1 if tally["FAILED"] else 0


if __name__ == "__main__":
    sys.exit(main())
