#!/usr/bin/env python3
"""Screen the frozen e422 checkpoint rows for high-confidence/low-lDDT cases.

Uses only Python's standard library and the batch checkpoint JSON files.
Run from outside the repository with ``python -I``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median


AUTHOR_FLAGS = ["Q9HWK6", "Q5JF22", "P07342", "A0A4Q0WMG2"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoints", type=Path, nargs="+", required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    by_accession = defaultdict(list)
    by_chain = defaultdict(list)
    input_meta = []
    n_rows = 0
    for path in args.checkpoints:
        raw = path.read_bytes()
        doc = json.loads(raw)
        input_meta.append({"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                           "n_rows": len(doc.get("rows", [])), "batch": doc.get("batch_index")})
        for row in doc.get("rows", []):
            if any(row.get(k) is None for k in ("accession", "entry", "chain", "plddt", "lddt")):
                continue
            accession = row["accession"]
            entry = row["entry"].lower()
            chain = row["chain"]
            item = {
                "entry": entry,
                "chain": chain,
                "plddt": float(row["plddt"]),
                "lddt": float(row["lddt"]),
                "d_kabsch": float(row["d_kabsch"]) if row.get("d_kabsch") is not None else None,
                "uniprot_pos": int(row["uniprot_pos"]) if row.get("uniprot_pos") is not None else None,
            }
            by_accession[accession].append(item)
            by_chain[(accession, entry, chain)].append(item)
            n_rows += 1

    accession_rows = []
    accession_suspects = []
    for accession in sorted(by_accession):
        rows = by_accession[accession]
        item = {
            "accession": accession,
            "n_rows": len(rows),
            "n_entries": len({r["entry"] for r in rows}),
            "n_entry_chains": len({(r["entry"], r["chain"]) for r in rows}),
            "median_plddt": round(median(r["plddt"] for r in rows), 4),
            "median_lddt": round(median(r["lddt"] for r in rows), 6),
        }
        accession_rows.append(item)
        if item["median_plddt"] >= 90.0 and item["median_lddt"] < 0.4:
            accession_suspects.append(item)

    chain_rows = []
    for key in sorted(by_chain):
        accession, entry, chain = key
        rows = by_chain[key]
        chain_rows.append({
            "accession": accession,
            "entry": entry,
            "chain": chain,
            "n_rows": len(rows),
            "median_plddt": round(median(r["plddt"] for r in rows), 4),
            "median_lddt": round(median(r["lddt"] for r in rows), 6),
            "median_d_kabsch_A": round(median(r["d_kabsch"] for r in rows if r["d_kabsch"] is not None), 4),
            "uniprot_range": [min(r["uniprot_pos"] for r in rows if r["uniprot_pos"] is not None),
                              max(r["uniprot_pos"] for r in rows if r["uniprot_pos"] is not None)],
        })

    chain_suspects = [r for r in chain_rows if r["n_rows"] >= 30
                      and r["median_plddt"] >= 90.0 and r["median_lddt"] < 0.4]
    by_entry = defaultdict(list)
    for row in chain_rows:
        by_entry[(row["accession"], row["entry"])].append(row)
    inconsistent_entries = []
    for (accession, entry), rows in sorted(by_entry.items()):
        if len(rows) < 2:
            continue
        eligible = [r for r in rows if r["n_rows"] >= 30 and r["median_plddt"] >= 90.0]
        if len(eligible) < 2:
            continue
        lo = min(eligible, key=lambda r: r["median_lddt"])
        hi = max(eligible, key=lambda r: r["median_lddt"])
        gap = hi["median_lddt"] - lo["median_lddt"]
        if gap >= 0.4:
            inconsistent_entries.append({
                "accession": accession,
                "entry": entry,
                "chain_median_lddt_gap": round(gap, 6),
                "low_chain": lo,
                "high_chain": hi,
            })

    author_flag_detail = {
        accession: {
            "accession_summary": next((r for r in accession_rows if r["accession"] == accession), None),
            "chains": [r for r in chain_rows if r["accession"] == accession],
            "inconsistent_entries": [r for r in inconsistent_entries if r["accession"] == accession],
        }
        for accession in AUTHOR_FLAGS
    }
    result = {
        "method": {
            "accession_signature": "median pLDDT >= 90 and median lDDT < 0.4 across all ledger rows for accession",
            "chain_signature": "at least 30 rows, median pLDDT >= 90, median lDDT < 0.4 for one entry/chain",
            "chain_inconsistency": "at least two >=30-row chains for the same accession/entry with median pLDDT >=90 and median lDDT gap >=0.4",
            "note": "The 0.4 gap is a transparent audit-screen threshold; flagged rows are candidates, not proof of mapping error.",
        },
        "inputs": input_meta,
        "total_valid_rows": n_rows,
        "n_accessions": len(accession_rows),
        "accession_suspects": accession_suspects,
        "chain_suspects": chain_suspects,
        "inconsistent_entry_suspects": inconsistent_entries,
        "author_flag_detail": author_flag_detail,
        "all_accession_summaries": accession_rows,
        "all_chain_summaries": chain_rows,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "n_rows": n_rows,
        "n_accessions": len(accession_rows),
        "accession_suspects": accession_suspects,
        "chain_suspects": chain_suspects,
        "inconsistent_entry_suspects": inconsistent_entries,
        "author_flag_detail": author_flag_detail,
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
