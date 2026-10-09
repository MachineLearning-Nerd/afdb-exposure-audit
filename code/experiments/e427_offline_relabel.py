"""e427 offline relabel harness (no network; reads only results/e422/raw).

For every (entry, chain, accession) triple that contributed rows to the e422
batch checkpoints, re-derive the labels twice from the cached raw bytes:
  OLD = frozen e423_runner.parse_mappings + label_triple (must reproduce the
        ledger rows exactly — a harness self-check), and
  NEW = e427_label_repair (label_seq_id mapping + identity check).
Writes a comparison JSON to --out (never under results/e422).

Usage: .venv/bin/python -I experiments/e427_offline_relabel.py --out <path>
       [--rows-out <path>]   # optional: full NEW rows (for a later,
                             # separately registered P5 rerun)
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import e423_runner as old  # noqa: E402
import e427_label_repair as new  # noqa: E402

BATCHES = (1, 2, 3, 4)


def load(raw_dir, name):
    with open(os.path.join(raw_dir, name), "rb") as fh:
        data = fh.read()
    if hashlib.sha256(data).hexdigest() != name:
        raise RuntimeError(f"raw cache hash mismatch: {name}")
    return data


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--rows-out", default=None)
    args = ap.parse_args(argv)
    raw_dir = os.path.join(ROOT, "results/e422/raw")

    url_raw, ledger = {}, collections.defaultdict(list)
    for b in BATCHES:
        ck = json.load(open(os.path.join(ROOT, f"results/e422/batches/e422_batch{b:02d}_checkpoint.json")))
        for r in ck.get("receipts", []):
            if r.get("status") == 200 and r.get("raw_name"):
                url_raw[r["url"]] = r["raw_name"]
        for row in ck.get("rows", []):
            ledger[(row["entry"], row["chain"], row["accession"])].append(row)

    af_url = {}
    for u in url_raw:
        if u.startswith("https://alphafold.ebi.ac.uk/files/AF-") and u.endswith(".cif"):
            acc = u.split("AF-")[1].split("-F1")[0]
            af_url[acc] = u

    entries = sorted({e for e, _c, _a in ledger})
    out = {"schema": "e427-offline-relabel-v1", "triples": [], "harness_self_check": {}}
    new_rows_all = []
    mism = 0
    for entry in entries:
        map_raw = url_raw[f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry}"]
        ref_raw = url_raw[f"https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"]
        payload = load(raw_dir, map_raw)
        ref_text = load(raw_dir, ref_raw).decode("utf-8", errors="replace")
        old_tr = old.parse_mappings(payload, entry)
        try:
            new_tr = new.parse_mappings_labelspace(payload, entry)
            new_err = None
        except Exception as exc:  # typed CensusError
            new_tr, new_err = [], f"{type(exc).__name__}: {exc}"
        for t in old_tr:
            key = (entry, t["chain_id"], t["accession"])
            if key not in ledger:
                continue
            acc = t["accession"]
            af_text = load(raw_dir, url_raw[af_url[acc]]).decode("utf-8", errors="replace")
            o = old.label_triple(t, entry, af_text, ref_text, "", "")
            # harness self-check: OLD must reproduce the ledger rows of this segment
            led = {(r["uniprot_pos"]): r for r in ledger[key]}
            ok_self = all(led.get(r["uniprot_pos"]) == r for r in o["rows"])
            mism += 0 if ok_self else 1
            # pair NEW segment by chain/accession and UniProt start
            nts = [x for x in new_tr if x["chain_id"] == t["chain_id"] and x["accession"] == acc
                   and x["segment"]["unp"][0] == min(t["res_map"].values())]
            n = new.label_triple_labelspace(nts[0], entry, af_text, ref_text, "", "") if nts else None
            if n and n["status"] == "ok":
                new_rows_all += n["rows"]

            def summ(res):
                if res is None:
                    return {"status": "NO_SEGMENT" if not new_err else new_err}
                rws = res["rows"]
                return {"status": res["status"], "n_rows": len(rws),
                        "n_identity_mismatch": res.get("n_identity_mismatch"),
                        "median_lddt": float(np.median([r["lddt"] for r in rws if r["lddt"] is not None])) if rws else None,
                        "median_d": float(np.median([r["d_kabsch"] for r in rws])) if rws else None,
                        "median_plddt": float(np.median([r["plddt"] for r in rws if r["plddt"] is not None])) if rws else None}
            out["triples"].append({"entry": entry, "chain": t["chain_id"], "accession": acc,
                                   "old_reproduces_ledger": ok_self, "old": summ(o), "new": summ(n)})
    out["harness_self_check"] = {"triples": len(out["triples"]), "old_mismatch_with_ledger": mism}
    json.dump(out, open(args.out, "w"), indent=1, sort_keys=True)
    if args.rows_out:
        json.dump(new_rows_all, open(args.rows_out, "w"))
    print(json.dumps(out["harness_self_check"]))


if __name__ == "__main__":
    main()
