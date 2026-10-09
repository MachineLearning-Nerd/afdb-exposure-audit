"""Download the PDBe *updated* mmCIF (residue-level SIFTS annotation) for every entry of the corrected
census ledger and compare each file with the SHA-256 recorded for the cross-check run.

PDBe regenerates updated mmCIF files weekly, so a later download can differ
from the bytes the recorded cross-check used (updated_mmcif_sha256.json); the
script reports MATCH / CHANGED per entry and keeps whatever it downloads.

Usage: python -I fetch_updated_mmcif.py --derived data/derived --out-dir _updated_mmcif
Network: yes (https://www.ebi.ac.uk/pdbe/entry-files/download/<entry>_updated.cif.gz), 0.3 s between requests.
"""
import argparse
import gzip
import hashlib
import json
import time
import urllib.request
from pathlib import Path

URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{e}_updated.cif.gz"
HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data/derived")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recorded = json.loads((HERE / "updated_mmcif_sha256.json").read_text())["files"]
    entries = set()
    for b in (2, 3, 4):
        ck = json.loads((Path(a.derived) / f"results/e427/batches/e422_batch{b:02d}_checkpoint.json").read_text())
        entries |= {r["entry"].lower() for r in ck["rows"]}
        entries |= {f["entry"].lower() for f in ck["exclusions"]["triple_failures"]}
    counts = {"MATCH": 0, "CHANGED": 0, "FAILED": 0, "NOT_RECORDED": 0}
    for e in sorted(entries):
        p = out / f"{e}_updated.cif"
        if not p.exists():
            try:
                req = urllib.request.Request(URL.format(e=e), headers={"User-Agent": "afdb-exposure-audit-crosscheck/1"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    p.write_bytes(gzip.decompress(r.read()))
            except Exception as exc:  # noqa: BLE001 - reported per entry
                counts["FAILED"] += 1
                print(f"FAILED  {e}: {exc}")
                continue
            time.sleep(0.3)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        want = recorded.get(f"{e}_updated.cif")
        state = "NOT_RECORDED" if want is None else "MATCH" if got == want else "CHANGED"
        counts[state] += 1
        if state != "MATCH":
            print(f"{state:8} {e}")
    print(counts)


if __name__ == "__main__":
    main()
