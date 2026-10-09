"""Full corrected-ledger cross-check with a different parser, lDDT implementation and mapping source.

Every (entry, chain, accession) group of the corrected census ledger
(results/e427/batches, batches 2-4), plus the groups the study rejected for
too few reference CA atoms / aligned residues / sequence identity, is
recomputed with biotite's mmCIF parser and biotite.structure.lddt (CA atoms,
15 A inclusion radius, 0.5/1/2/4 A thresholds), using the residue-level
SIFTS UniProt mapping of the PDBe *updated* mmCIF (pdbx_sifts_xref_db_acc /
pdbx_sifts_xref_db_num) instead of the segment-level PDBe API mapping the
study used. No study or audit code is imported.

Reference coordinates and AFDB models are the census's own source bytes
(SHA-256 verified against the checkpoint receipts). Reference CA policy:
model 1, altloc '.', '?', '' or 'A', the 20 standard amino acids plus MSE.

Inputs:
  --derived   release data/derived (corrected checkpoints and receipts)
  --raw       directory with the census raw files named by SHA-256
              (python -I data/fetch/fetch_raw.py --set e422 --work _work
              -> _work/results/e422/raw)
  --updated   directory with <entry>_updated.cif (fetch_updated_mmcif.py)
  --out       output JSON (one record per group)

Output fields per group: bt_n_aligned / bt_n_mismatch (biotite pairs and
residue-identity mismatches), ledger_n, only_ledger / only_biotite (UniProt
positions present in one side only), max_abs_diff_segment and n_diff_gt_1e-4
(per-residue lDDT, ledger vs biotite, per SIFTS segment as in the study),
max_abs_diff_pooled_chain, and the chain medians of lDDT. No pLDDT value is
written; the AFDB B-factor column is not read.

Usage: python -I full_ledger_biotite.py --derived data/derived --raw _work/results/e422/raw \
           --updated _updated_mmcif --out full_ledger_out.json
"""
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as px

AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
FAIL_CODES = ("E422_REF_INSUFFICIENT_CA", "E422_INSUFFICIENT_ALIGNED", "E427_SEQUENCE_IDENTITY_FAILED")


def table(path):
    a = px.CIFFile.read(str(path)).block["atom_site"]
    return {k: a[k].as_array() for k in a.keys()}


def lddt_ca(R, M):
    arr = struc.AtomArray(len(R))
    arr.coord = R
    arr.res_id = np.arange(len(R)) + 1
    arr.chain_id[:] = "A"
    return struc.lddt(arr, M, aggregation="atom")


def receipts_and_ledger(derived):
    rc, ledger, fails = {}, defaultdict(list), []
    for b in range(2, 5):
        ck = json.loads((Path(derived) / f"results/e427/batches/e422_batch{b:02d}_checkpoint.json").read_text())
        for r in ck["receipts"]:
            if r.get("raw_name") and r.get("status") == 200:
                rc.setdefault(r["url"], r["raw_name"])
        for r in ck["rows"]:
            ledger[(r["entry"].lower(), r["chain"], r["accession"])].append(r)
        for f in ck["exclusions"]["triple_failures"]:
            if f["code"] in FAIL_CODES:
                fails.append(((f["entry"].lower(), f["chain"], f["accession"]), f["code"]))
    return rc, ledger, fails


def raw_file(raw, rc, url):
    name = rc[url]
    p = Path(raw) / name
    if hashlib.sha256(p.read_bytes()).hexdigest() != name:
        raise RuntimeError(f"SHA-256 mismatch for {url}")
    return p


def af_ca(raw, rc, acc):
    af = table(raw_file(raw, rc, f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif"))
    d = {}
    for i in np.nonzero((af["label_atom_id"] == "CA") & (af["pdbx_PDB_model_num"] == "1"))[0]:
        d[int(af["auth_seq_id"][i])] = (AA.get(af["label_comp_id"][i], "X"),
                                        np.array([float(af[c][i]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")]))
    return d


def sifts_residues(upd, chain, acc):
    """(label_asym_id, label_seq_id) -> UniProt position, model 1 CA atoms of the auth chain."""
    out = {}
    sel = np.nonzero((upd["label_atom_id"] == "CA") & (upd["pdbx_sifts_xref_db_acc"] == acc)
                     & (upd["auth_asym_id"] == chain) & (upd["pdbx_PDB_model_num"] == "1"))[0]
    for i in sel:
        if upd["pdbx_sifts_xref_db_num"][i] in (".", "?"):
            continue
        out[(upd["label_asym_id"][i], int(upd["label_seq_id"][i]))] = int(upd["pdbx_sifts_xref_db_num"][i])
    return out


def reference_ca(ref, chain):
    m = ((ref["label_atom_id"] == "CA") & (ref["pdbx_PDB_model_num"] == "1")
         & np.isin(ref["label_comp_id"], list(AA)) & np.isin(ref["label_alt_id"], [".", "?", "", "A"])
         & (ref["auth_asym_id"] == chain))
    ca = {}
    for i in np.nonzero(m)[0]:
        if ref["label_seq_id"][i] in (".", "?"):
            continue
        ca.setdefault((ref["label_asym_id"][i], int(ref["label_seq_id"][i])), i)
    return ca


def pair_residues(ref, ca, sifts, afd):
    pairs, mism = [], 0
    for key in sorted(ca):
        p = sifts.get(key)
        if p is None or p not in afd:
            continue
        i = ca[key]
        if AA[ref["label_comp_id"][i]] != afd[p][0]:
            mism += 1
            continue
        pairs.append((p, key, np.array([float(ref[c][i]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")])))
    return pairs, mism


def segment_lddt(pairs, segs, afd):
    """Per-residue lDDT computed per SIFTS segment (mirrors the study: one computation per segment)."""
    out = {}
    for s in segs:
        lo, hi = s["start"]["residue_number"], s["end"]["residue_number"]
        sp = [x for x in pairs if x[1][0] == s["struct_asym_id"] and lo <= x[1][1] <= hi]
        if len(sp) >= 2:
            R = np.array([x[2] for x in sp])
            M = np.array([afd[x[0]][1] for x in sp])
            for x, v in zip(sp, lddt_ca(R, M)):
                out[x[0]] = float(v)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data/derived")
    ap.add_argument("--raw", required=True)
    ap.add_argument("--updated", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rc, ledger, fails = receipts_and_ledger(a.derived)
    groups = [(k, "ok") for k in sorted(ledger)] + fails
    cache, afcache, results = {}, {}, []
    for (entry, chain, acc), status in groups:
        rec = {"entry": entry, "chain": chain, "acc": acc, "study_status": status}
        try:
            if entry not in cache:
                upd = Path(a.updated) / f"{entry}_updated.cif"
                mp = json.loads(raw_file(a.raw, rc, f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry}")
                                .read_text())
                cache = {entry: (table(raw_file(a.raw, rc, f"https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif")),
                                 table(upd) if upd.exists() else None, mp)}
            ref, upd, mp = cache[entry]
            if upd is None:
                rec["error"] = "no updated cif"
                results.append(rec)
                continue
            if acc not in afcache:
                afcache = {acc: af_ca(a.raw, rc, acc)}
            afd = afcache[acc]
            segs = [m for m in mp[entry]["UniProt"][acc]["mappings"] if str(m["chain_id"]) == chain]
            rec["n_segments"] = len(segs)
            pairs, mism = pair_residues(ref, reference_ca(ref, chain), sifts_residues(upd, chain, acc), afd)
            rec["bt_n_aligned"], rec["bt_n_mismatch"] = len(pairs), mism
            if len(pairs) < 2:
                results.append(rec)
                continue
            seg_l = segment_lddt(pairs, segs, afd)
            R = np.array([x[2] for x in pairs])
            M = np.array([afd[x[0]][1] for x in pairs])
            pooled = {x[0]: float(v) for x, v in zip(pairs, lddt_ca(R, M))}
            rows = ledger.get((entry, chain, acc), [])
            lp = {r["uniprot_pos"]: r["lddt"] for r in rows}
            rec["ledger_n"] = len(rows)
            rec["dup_positions_in_ledger"] = len(rows) - len(lp)
            common = set(lp) & set(seg_l)
            rec["only_ledger"] = len(set(lp) - set(seg_l))
            rec["only_biotite"] = len(set(seg_l) - set(lp))
            diffs = [abs(lp[p] - seg_l[p]) for p in common if lp[p] is not None]
            rec["max_abs_diff_segment"] = max(diffs) if diffs else None
            rec["n_diff_gt_1e-4"] = sum(d > 1.0001e-4 for d in diffs)
            pd = [abs(lp[p] - pooled[p]) for p in set(lp) & set(pooled) if lp[p] is not None]
            rec["max_abs_diff_pooled_chain"] = max(pd) if pd else None
            if rows:
                rec["ledger_median"] = float(np.median([r["lddt"] for r in rows if r["lddt"] is not None]))
                rec["bt_median_segment"] = float(np.median([seg_l[p] for p in common])) if common else None
                rec["bt_median_pooled"] = float(np.median(list(pooled.values())))
        except Exception as exc:  # noqa: BLE001 - recorded per group, counted in the summary
            rec["error"] = f"{type(exc).__name__}: {exc}"
        results.append(rec)
    Path(a.out).write_text(json.dumps(results, indent=0))
    ok = [r for r in results if r["study_status"] == "ok" and "error" not in r]
    exact = [r for r in ok if r.get("n_diff_gt_1e-4") == 0 and r.get("only_ledger") == 0 and r.get("only_biotite") == 0]
    print("groups", len(results), "ledger segments compared", len(ok), "errors", sum("error" in r for r in results))
    print("identical positions and |delta lDDT| <= 1e-4:", len(exact), "of", len(ok),
          "| rows", sum(r["ledger_n"] for r in exact), "of", sum(r["ledger_n"] for r in ok),
          "| max |delta| among them", max((r["max_abs_diff_segment"] for r in exact), default=None))
    for r in results:
        if "error" in r or r["study_status"] != "ok" or r.get("n_diff_gt_1e-4") or r.get("only_ledger") \
                or r.get("only_biotite") or r.get("dup_positions_in_ledger"):
            print(json.dumps(r))


if __name__ == "__main__":
    main()
