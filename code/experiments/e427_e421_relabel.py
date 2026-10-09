"""e427 §4.5: corrected relabel of the exploratory e421 census.

Re-fetches the 168 qualifying entries (the e421 run kept no raw bytes) from
the same public lanes, persists every body under data/e427/e421_raw/<sha256>
with a receipt, and labels each entry's chain A / qualifying accession with
the e427 label-space mapping + identity check.  The e421 label rules are kept
(verbatim copies below): Kabsch, lDDT-Calpha (15 A, 0.5/1/2/4 A), robust z over
the ALIGNED residues' B-factors, E/F thresholds, four cells.  Deviation
(disclosed): altlocs other than '.', '?', '' and 'A' are dropped (the e421
parser kept them, which could duplicate a residue).

Outputs (results/e427/e421_relabel/): census_results.json (same per-entry
schema as results/e421/census_results.json), receipts.json, summary.json.
Usage: .venv/bin/python experiments/e427_e421_relabel.py
"""
from __future__ import annotations

import collections
import datetime
import hashlib
import json
import math
import os
import random
import sys
import time

import numpy as np
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import e421_pipeline as pipe  # noqa: E402
import e427_label_repair as repair  # noqa: E402

RAW = os.path.join(ROOT, "data/e427/e421_raw")
OUT = os.path.join(ROOT, "results/e427/e421_relabel")
UA = {"User-Agent": "moluq-e427-e421-relabel/1.0 (research; public endpoints)"}
CHAIN = "A"
MIN_ALIGNED = 30
receipts: list[dict] = []


def fetch(url: str) -> bytes | None:
    time.sleep(max(0.0, 0.25 + random.uniform(-0.1, 0.1)))
    for attempt in range(1, 4):
        try:
            r = requests.get(url, headers=UA, timeout=60)
            body = r.content
            sha = hashlib.sha256(body).hexdigest()
            receipts.append({"url": url, "status": r.status_code, "sha256": sha,
                             "bytes_len": len(body), "attempt": attempt,
                             "utc": datetime.datetime.now(datetime.timezone.utc).isoformat()})
            if r.status_code == 200:
                path = os.path.join(RAW, sha)
                if not os.path.exists(path):
                    with open(path, "wb") as fh:
                        fh.write(body)
                return body
            if r.status_code not in (429,) and r.status_code < 500:
                return None
        except requests.RequestException as exc:
            receipts.append({"url": url, "status": None, "error": type(exc).__name__,
                             "attempt": attempt,
                             "utc": datetime.datetime.now(datetime.timezone.utc).isoformat()})
        time.sleep(2 ** attempt)
    return None


# --- e421 label rules (verbatim from experiments/e421_census_fixed.py) ------
def kabsch(mobile, ref):
    M = np.array(mobile); R = np.array(ref)
    if M.shape != R.shape or M.shape[0] < 3: return None, None  # noqa: E701
    mc = np.array([math.fsum(M[:, j]) for j in range(3)]) / M.shape[0]
    rc = np.array([math.fsum(R[:, j]) for j in range(3)]) / R.shape[0]
    H = (M - mc).T @ (R - rc)
    U, S, Vt = np.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False)
    det = np.linalg.det(U @ Vt)
    sup = (M - mc) @ (U @ np.diag([1, 1, 1 if det > 0 else -1]) @ Vt) + rc
    return sup, np.linalg.norm(sup - R, axis=1)


def lddt(sup, ref, cutoff_sq=225.0, thresholds=(0.5, 1.0, 2.0, 4.0)):
    n = ref.shape[0]
    if n < 2: return np.full(n, np.nan)  # noqa: E701
    diff = ref[:, None, :] - ref[None, :, :]
    dsq = np.sum(diff ** 2, axis=-1); np.fill_diagonal(dsq, np.inf)  # noqa: E702
    nb = dsq < cutoff_sq
    sd = np.sqrt(np.sum((sup[:, None, :] - sup[None, :, :]) ** 2, axis=-1))
    rd = np.sqrt(dsq)
    pres = np.zeros((n, n, len(thresholds)), dtype=bool)
    for ti, t in enumerate(thresholds): pres[:, :, ti] = np.abs(sd - rd) < t  # noqa: E701
    pres &= nb[:, :, None]
    nc = nb.sum(axis=1).astype(float)
    return np.where(nc > 0, pres.sum(axis=2).sum(axis=1) / (nc * len(thresholds)), np.nan)


def bfactor_z(bfs):
    a = np.array(bfs); med = np.median(a)  # noqa: E702
    q25, q75 = np.percentile(a, [25, 75]); iqr = q75 - q25 or 1.0  # noqa: E702
    return (a - med) / iqr


def cell_of(l, d, zj):
    e = "correct" if (l >= 0.60 and d <= 4.0) else ("error" if (l < 0.60 and d > 4.0) else "U")
    f = "flexible" if zj >= 1.0 else ("ordered" if zj <= -1.0 else "U")
    cell = "U"
    if e != "U" and f != "U":
        cell = {("correct", "ordered"): "OC", ("correct", "flexible"): "A",
                ("error", "ordered"): "B", ("error", "flexible"): "O"}[(e, f)]
    return e, f, cell


def label_entry(pdb_id, acc, cif_url):
    row = {"entry": pdb_id, "accession": acc}
    mp = fetch(f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}")
    if mp is None:
        return {**row, "status": "mapping_fail"}
    try:
        triples = repair.parse_mappings_labelspace(mp, pdb_id)
    except Exception as exc:  # typed CensusError
        return {**row, "status": "mapping_parse_failed", "detail": str(exc)[:200]}
    segs = [t for t in triples if t["accession"] == acc and t["chain_id"] == CHAIN]
    if not segs:
        return {**row, "status": "mapping_acc_absent"}
    af = fetch(cif_url)
    if af is None:
        return {**row, "status": "af_fetch_fail"}
    ref = fetch(f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id}.cif")
    if ref is None:
        return {**row, "status": "ref_fetch_fail"}
    af_atoms = pipe.parse_cif_ca(af.decode("utf-8", errors="replace"))
    af_by_pos = {a["residue_number"]: a for a in af_atoms}
    ref_text = ref.decode("utf-8", errors="replace")
    pairs = []  # (ref_atom, af_atom, unp)
    n_pairs = n_mm = 0
    seen = set()
    for t in segs:
        res_map = t["res_map"]
        if res_map is None:
            (ls, le), (us, ue) = t["segment"]["label"], t["segment"]["unp"]
            ent = repair.parse_entity_poly_seq(ref_text, t["entity_id"])
            res_map = repair.align_map([(k, ent[k]) for k in range(ls, le + 1) if k in ent],
                                       [(p, af_by_pos[p]["aa"]) for p in range(us, ue + 1) if p in af_by_pos])
        for a in repair.parse_cif_ca_labelspace(ref_text, t["struct_asym_id"]):
            unp = res_map.get(a["label_seq_id"])
            if unp is None or unp in seen or unp not in af_by_pos:
                continue
            n_pairs += 1
            if af_by_pos[unp]["aa"] != a["aa"]:
                n_mm += 1
                continue
            seen.add(unp)
            pairs.append((a, af_by_pos[unp], unp))
    if n_pairs and n_mm / n_pairs > repair.MAX_MISMATCH_FRACTION:
        return {**row, "status": "sequence_identity_failed", "n_identity_mismatch": n_mm}
    if len(pairs) < MIN_ALIGNED:
        return {**row, "status": "insufficient_aligned", "n_aligned": len(pairs)}
    ref_ca = [[a["x"], a["y"], a["z"]] for a, _, _ in pairs]
    af_ca = [[b["x"], b["y"], b["z"]] for _, b, _ in pairs]
    sup, dists = kabsch(af_ca, ref_ca)
    if sup is None:
        return {**row, "status": "kabsch_fail"}
    lv = lddt(sup, np.array(ref_ca))
    z = bfactor_z([a["b_factor"] for a, _, _ in pairs])
    labels = []
    for j, (a, b, unp) in enumerate(pairs):
        l = lv[j] if not np.isnan(lv[j]) else -1
        e, f, cell = cell_of(l, dists[j], z[j])
        labels.append({"accession": acc, "entry": pdb_id, "uniprot_pos": unp, "e": e, "f": f,
                       "cell": cell, "lddt": round(float(lv[j]), 4), "dist": round(float(dists[j]), 4),
                       "b": a["b_factor"], "z": round(float(z[j]), 4), "plddt": b["b_factor"]})
    return {**row, "status": "ok", "labels": labels, "n_aligned": len(pairs),
            "n_identity_mismatch": n_mm,
            "n_cells": dict(collections.Counter(x["cell"] for x in labels))}


def summarise(results, labs=None):
    ok = [r for r in results.values() if r.get("status") == "ok"]
    if labs is None:
        labs = [x for r in ok for x in r["labels"]]
    cells = collections.Counter(x["cell"] for x in labs)
    uniq = {c: len({(x["accession"], x["uniprot_pos"]) for x in labs if x["cell"] == c}) for c in cells}
    bins = [(0, 50), (50, 70), (70, 80), (80, 90), (90, 100.01)]
    grad = []
    for lo, hi in bins:
        sel = [x for x in labs if x.get("lddt") is not None and lo <= x["plddt"] < hi
               and not math.isnan(x["lddt"])]
        grad.append({"bin": [lo, hi], "n": len(sel),
                     "err_rate_lddt_lt_0.60": (sum(x["lddt"] < 0.60 for x in sel) / len(sel)) if sel else None})
    from scipy.stats import spearmanr
    fin = [x for x in labs if x.get("lddt") is not None and not math.isnan(x["lddt"])]
    rho = spearmanr([x["plddt"] for x in fin], [x["lddt"] for x in fin]).statistic if fin else None
    return {"entries": len(results), "status_counts": dict(collections.Counter(r["status"] for r in results.values())),
            "ok_entries": len(ok), "ok_accessions": len({r["accession"] for r in ok}),
            "residue_rows": len(labs), "cells_rows": dict(cells), "cells_unique": uniq,
            "gradient": grad, "spearman_plddt_lddt": rho, "n_spearman": len(fin)}


def main():
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    qual = json.load(open(os.path.join(ROOT, "results/e421/qualifying_entries.json")))
    results = {}
    for i, q in enumerate(qual):
        try:
            results[q["pdb_id"]] = label_entry(q["pdb_id"], q["afdb_accession"], q["cif_url"])
        except Exception as exc:  # noqa: BLE001 - recorded, never silent
            results[q["pdb_id"]] = {"entry": q["pdb_id"], "accession": q["afdb_accession"],
                                    "status": "exception", "detail": f"{type(exc).__name__}: {exc}"[:200]}
        if (i + 1) % 25 == 0:
            print(f"{i + 1}/{len(qual)}", flush=True)
    old = json.load(open(os.path.join(ROOT, "results/e421/census_results.json")))
    transitions = collections.Counter((old.get(k, {}).get("status"), v["status"]) for k, v in results.items())
    summary = {"schema": "e427-e421-relabel-v1", "registration": "docs/e427_registration.md section 4.5",
               "corrected": summarise(results),
               "original": summarise(old, json.load(open(os.path.join(
                   ROOT, "results/e421/labeled_residues_with_plddt.json")))),
               "status_transitions": {f"{a} -> {b}": n for (a, b), n in sorted(transitions.items(), key=str)},
               "n_receipts": len(receipts)}
    json.dump(results, open(os.path.join(OUT, "census_results.json"), "w"), indent=1, default=str)
    json.dump(receipts, open(os.path.join(OUT, "receipts.json"), "w"), indent=1)
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=str)
    print(json.dumps({k: summary[k] for k in ("status_transitions",)}, indent=1))
    print(json.dumps({k: summary["corrected"][k] for k in ("ok_entries", "residue_rows", "cells_rows")}))


if __name__ == "__main__":
    main()
