
"""e421 fixed pipeline: SIFTS-based residue mapping + E/F labels + four-cell census.

Uses PDBe mapping responses (already in census raw bytes) to align reference
CA atoms to AF model CA atoms via UniProt positions. Computes E (Kabsch +
lDDT) and F (B-factor z-score) per frozen label rules. Produces the
four-cell census and G1 verdict.
"""
import argparse, hashlib, json, math, os, sys, time, random
from collections import Counter, defaultdict
import numpy as np
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS_RAW = os.path.join(ROOT, "data/e420/census_raw")
LABELS_RAW = os.path.join(ROOT, "data/e420/labels_raw")
OUT_DIR = os.path.join(ROOT, "results/e421")
UA = {"User-Agent": "moluq-e421-pipeline/1.0"}
AA3_TO_1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
            "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
            "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
AA1 = set("ACDEFGHIKLMNPQRSTVWY")

def sha(b): return hashlib.sha256(b).hexdigest()

def tokenize_cif(text):
    tokens = []; i = 0; n = len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n": i += 1; continue
        if c == "#":
            while i < n and text[i] != "\n": i += 1
            continue
        if c == ";" and (i == 0 or text[i-1] == "\n"):
            end = text.find("\n;", i)
            if end == -1: break
            tokens.append(text[i+1:end]); i = end + 2; continue
        if c in "'\"":
            q = c; end = i + 1
            while True:
                end = text.find(q, end)
                if end == -1: break
                if text[end+1:end+2] in ("", " ", "\t", "\r", "\n"): break
                end += 1
            tokens.append(text[i+1:end]); i = end + 1; continue
        s = i
        while i < n and text[i] not in " \t\r\n": i += 1
        tokens.append(text[s:i])
    return tokens

def parse_cif_ca(text, chain_filter=None):
    """Extract CA atoms (residue_number, x, y, z, b_factor, comp_id) from mmCIF."""
    tokens = tokenize_cif(text)
    atoms = []; idx = 0
    while idx < len(tokens):
        if tokens[idx] == "loop_":
            idx += 1; tags = []
            while idx < len(tokens) and tokens[idx].startswith("_"):
                tags.append(tokens[idx]); idx += 1
            if any("_atom_site." in t for t in tags):
                col = {t.split(".",1)[1]: ci for ci,t in enumerate(tags)}
                for k in ("label_atom_id","label_comp_id","auth_seq_id","Cartn_x","Cartn_y","Cartn_z",
                          "B_iso_or_equiv","pdbx_PDB_model_num","label_alt_id","auth_asym_id"):
                    if k not in col: break
                else:
                    while idx < len(tokens):
                        t = tokens[idx]
                        if t in ("loop_","stop_") or t.startswith("_") or t.startswith("data_"): break
                        row = tokens[idx:idx+len(tags)]
                        if len(row) < len(tags): break
                        idx += len(tags)
                        if row[col["label_atom_id"]] != "CA": continue
                        if row[col["label_comp_id"]] not in AA3_TO_1: continue
                        if row[col["pdbx_PDB_model_num"]] != "1": continue
                        alt = row[col["label_alt_id"]]
                        if alt not in (".","?","","A"): continue
                        if chain_filter and row[col["auth_asym_id"]] != chain_filter: continue
                        try:
                            atoms.append({"residue_number": int(row[col["auth_seq_id"]]),
                                          "x": float(row[col["Cartn_x"]]), "y": float(row[col["Cartn_y"]]),
                                          "z": float(row[col["Cartn_z"]]),
                                          "b_factor": float(row[col["B_iso_or_equiv"]]),
                                          "comp_id": row[col["label_comp_id"]],
                                          "aa": AA3_TO_1[row[col["label_comp_id"]]]})
                        except (ValueError, KeyError): continue
                    continue
            continue
        if tokens[idx].startswith("_") and idx + 1 < len(tokens): idx += 2; continue
        idx += 1
    return atoms

def kabsch_superpose(mobile, ref):
    M = np.array(mobile, dtype=np.float64); R = np.array(ref, dtype=np.float64)
    if M.shape != R.shape or M.shape[0] < 3: return None, None
    mc = np.array([math.fsum(M[:,j]) for j in range(3)]) / M.shape[0]
    rc = np.array([math.fsum(R[:,j]) for j in range(3)]) / R.shape[0]
    H = (M - mc).T @ (R - rc)
    U, S, Vt = np.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False)
    det = np.linalg.det(U @ Vt)
    D = np.diag([1.0, 1.0, 1.0 if det > 0 else -1.0])
    sup = (M - mc) @ (U @ D @ Vt) + rc
    dists = np.linalg.norm(sup - R, axis=1)
    return sup, dists

def lddt_calc(sup, ref, cutoff_sq=225.0, thresholds=(0.5,1.0,2.0,4.0)):
    n = ref.shape[0]
    if n < 2: return np.full(n, np.nan), np.zeros(n)
    diff = ref[:,None,:] - ref[None,:,:]
    dist_sq = np.sum(diff**2, axis=-1)
    np.fill_diagonal(dist_sq, np.inf)
    nb_mask = dist_sq < cutoff_sq
    sup_diff = sup[:,None,:] - sup[None,:,:]
    sup_dist = np.sqrt(np.sum(sup_diff**2, axis=-1))
    ref_dist = np.sqrt(dist_sq)
    preserved = np.zeros((n, n, len(thresholds)), dtype=bool)
    for ti, t in enumerate(thresholds):
        preserved[:,:,ti] = np.abs(sup_dist - ref_dist) < t
    preserved &= nb_mask[:,:,None]
    nb_count = nb_mask.sum(axis=1).astype(np.float64)
    lddt = np.where(nb_count > 0, preserved.sum(axis=2).sum(axis=1) / (nb_count * len(thresholds)), np.nan)
    return lddt, nb_count

def bfactor_z(bfs):
    arr = np.array(bfs, dtype=np.float64)
    med = float(np.median(arr))
    q25, q75 = np.percentile(arr, [25,75], method="linear")
    iqr = q75 - q25 or 1.0
    return (arr - med) / iqr, med, iqr

def process_entry(pdb_id, af_acc, chain_id, census_raw_dir, raw_dir):
    """Process one entry. Returns {status, cell_counts, label_rows} or raises."""
    entry_lower = pdb_id.lower()
    # 1. Load PDBe mapping from census raw
    mapping = None
    for f in os.listdir(census_raw_dir):
        path = os.path.join(census_raw_dir, f)
        try:
            d = json.load(open(path))
            if isinstance(d, dict) and entry_lower in d:
                uniprot = d[entry_lower].get("UniProt", {})
                if af_acc in uniprot:
                    mapping = uniprot[af_acc]
                    break
        except (json.JSONDecodeError, KeyError): continue
    if mapping is None:
        return {"status": "mapping_not_found"}
    # 2. Extract residue mapping: author_residue_number -> unp position
    res_map = {}
    for m in mapping.get("mappings", []):
        unp_s, unp_e = m["unp_start"], m["unp_end"]
        s = m["start"]["author_residue_number"] or m["start"]["residue_number"]
        e = m["end"]["author_residue_number"] or m["end"]["residue_number"]
        for offset in range(e - s + 1):
            auth_num = s + offset
            unp_pos = unp_s + offset
            res_map[auth_num] = unp_pos
    # 3. Fetch AF model
    af_url = f"https://alphafold.ebi.ac.uk/files/AF-{af_acc}-F1-model_v6.cif"
    import requests as req_mod
    r = req_mod.get(af_url, headers=UA, timeout=60)
    if r.status_code != 200:
        return {"status": "af_fetch_fail", "code": r.status_code}
    af_atoms = parse_cif_ca(r.text)
    if len(af_atoms) < 3: return {"status": "af_insufficient_ca"}
    # AF model: residue_number = UniProt position (1-based)
    af_by_pos = {a["residue_number"]: a for a in af_atoms}
    # 4. Fetch reference structure
    ref_url = f"https://www.ebi.ac.uk/pdbe/entry-files/{entry_lower}.cif"
    r2 = req_mod.get(ref_url, headers=UA, timeout=60)
    if r2.status_code != 200:
        return {"status": "ref_fetch_fail", "code": r2.status_code}
    ref_atoms = parse_cif_ca(r2.text, chain_filter=chain_id)
    if len(ref_atoms) < 30: return {"status": "ref_insufficient_ca", "n": len(ref_atoms)}
    # 5. Align via mapping: reference auth_residue_number -> UniProt position -> AF model
    ref_ca, af_ca, unp_positions = [], [], []
    for a in ref_atoms:
        unp = res_map.get(a["residue_number"])
        if unp is None: continue
        af_atom = af_by_pos.get(unp)
        if af_atom is None: continue
        ref_ca.append([a["x"], a["y"], a["z"]])
        af_ca.append([af_atom["x"], af_atom["y"], af_atom["z"]])
        unp_positions.append(unp)
    if len(ref_ca) < 30:
        return {"status": "insufficient_aligned", "n": len(ref_ca)}
    # 6. Kabsch + lDDT
    sup, dists = kabsch_superpose(af_ca, ref_ca)
    if sup is None: return {"status": "kabsch_fail"}
    lddt, nbr = lddt_calc(sup, np.array(ref_ca))
    # 7. E labels
    e_labels = []
    for j in range(len(ref_ca)):
        d = dists[j]; l = lddt[j] if not np.isnan(lddt[j]) else -1
        if l >= 0.60 and d <= 4.0: e_labels.append("correct")
        elif l < 0.60 and d > 4.0: e_labels.append("error")
        else: e_labels.append("U")
    # 8. F labels from B-factor z-score
    ref_bf = [a["b_factor"] for a in ref_atoms]
    z, med, iqr = bfactor_z(ref_bf)
    f_labels = []
    for j, a in enumerate(ref_atoms):
        if z[j] >= 1.0: f_labels.append("flexible")
        elif z[j] <= -1.0: f_labels.append("ordered")
        else: f_labels.append("U")
    # 9. Combined labels
    labels = []
    for j in range(len(ref_ca)):
        e, f = e_labels[j], f_labels[j]
        if e == "U" or f == "U": cell = "U"
        elif e == "correct" and f == "ordered": cell = "OC"
        elif e == "correct" and f == "flexible": cell = "A"
        elif e == "error" and f == "ordered": cell = "B"
        elif e == "error" and f == "flexible": cell = "O"
        else: cell = "U"
        labels.append({"uniprot_pos": unp_positions[j], "e": e, "f": f, "cell": cell,
                       "lddt": round(float(lddt[j]), 4) if not np.isnan(lddt[j]) else None,
                       "dist": round(float(dists[j]), 4) if not np.isnan(dists[j]) else None,
                       "b_factor": ref_bf[j] if j < len(ref_bf) else None,
                       "z_score": round(float(z[j]), 4) if j < len(z) else None})
    return {"status": "ok", "labels": labels, "n_aligned": len(ref_ca),
            "entry": entry_lower, "accession": af_acc, "chain": chain_id}

def run(entries, raw_dir, out_dir, max_entries=None):
    os.makedirs(out_dir, exist_ok=True); os.makedirs(raw_dir, exist_ok=True)
    total = min(max_entries or len(entries), len(entries))
    all_labels = []; proteins = {}
    print(f"Processing {total}/{len(entries)} entries")
    census_raw = os.path.join(ROOT, "data/e420/census_raw")
    for i, entry in enumerate(entries[:total]):
        pdb_id = entry["pdb_id"]; acc = entry["afdb_accession"]
        chain = entry.get("chain_id") or entry.get("chain") or "A"
        try:
            result = process_entry(pdb_id, acc, chain, census_raw, raw_dir)
        except Exception as e:
            result = {"status": "exception", "detail": str(e)[:200]}
        status = result.get("status")
        if status == "ok":
            for lbl in result["labels"]:
                all_labels.append({"accession": acc, "entry": pdb_id, "chain": chain, **lbl})
            proteins[acc] = {"cell": None, "status": "ok", "n_aligned": result["n_aligned"]}
        else:
            proteins[acc] = {"cell": None, "status": status}
        if (i+1) % 10 == 0:
            print(f"  {i+1}/{total} | labels so far: {len(all_labels)}")
    # four-cell census
    cells = Counter(l["cell"] for l in all_labels if l["cell"] != "U")
    acc_cells = defaultdict(set)
    for l in all_labels:
        if l["cell"] in ("A","B","OC","O"):
            acc_cells[l["accession"]].add((l["accession"], l["uniprot_pos"], l["cell"]))
    protein_cells = {}
    for acc, res_set in acc_cells.items():
        for _, _, cell in res_set:
            protein_cells[acc] = cell; break
    cell_proteins = Counter(protein_cells.values())
    cell_groups = defaultdict(set)
    # fold assignment: simplified for census (without full group clustering)
    for l in all_labels:
        if l["cell"] in ("A","B","OC","O"):
            cell_groups[l["cell"]].add(l["accession"])
    floor_ok = (cell_proteins.get("A",0) >= 20 and cell_proteins.get("B",0) >= 20)
    verdict = "E421_H1_SUPPORT_PASS" if floor_ok and cells.get("A",0) >= 200 and cells.get("B",0) >= 200 else "E421_UNDERPOWERED_POOL"
    print(f"\n{'='*60}")
    print(f"FOUR-CELL CENSUS:")
    for cell in ("OC","A","B","O"):
        print(f"  {cell}: rows={cells.get(cell,0)} proteins={cell_proteins.get(cell,0)}")
    print(f"  U: {sum(1 for l in all_labels if l['cell']=='U')}")
    print(f"  Floors satisfied (A>=20, B>=20): {floor_ok}")
    print(f"  VERDICT: {verdict}")
    # save
    artifact = {"schema": "e421-four-cell-census-v1", "verdict": verdict,
                "cap_bound": True, "floors_satisfied": floor_ok,
                "cells": {c: cells.get(c,0) for c in ("OC","A","B","O","U")},
                "proteins_per_cell": dict(cell_proteins) if isinstance(cell_proteins, Counter) else {},
                "label_count": len(all_labels), "entry_count": total}
    with open(os.path.join(out_dir, "four_cell_census.json"), "w") as f:
        json.dump(artifact, f, indent=2, sort_keys=True); f.write("\n")
    with open(os.path.join(out_dir, "label_ledger.json"), "w") as f:
        json.dump({"schema": "e421-label-ledger-v1", "rows": all_labels,
                   "summary": {"total_labels": len(all_labels), "cells": dict(cells)}}, f, indent=2, sort_keys=True)
    return verdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-entries", type=int, default=None)
    args = parser.parse_args()
    qual = json.load(open(os.path.join(ROOT, "results/e421/qualifying_entries.json")))
    # add chain_id from census if available
    for q in qual:
        q.setdefault("chain_id", "A")
    verdict = run(qual, os.path.join(ROOT, "data/e420/labels_raw"), os.path.join(ROOT, "results/e421"),
                  max_entries=args.max_entries)
    print(f"\nFINAL VERDICT: {verdict}")
