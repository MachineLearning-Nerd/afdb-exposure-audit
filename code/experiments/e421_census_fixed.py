import json

import hashlib, json, math, os, sys, time, random
from collections import Counter, defaultdict
import numpy as np
import requests

ROOT = os.getcwd()
AA3_TO_1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
            "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
            "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
UA = {"User-Agent": "moluq-e421-pipeline/1.0"}

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

def parse_cif_ca(text, chain=None):
    tokens = tokenize_cif(text)
    atoms = []; idx = 0
    while idx < len(tokens):
        if tokens[idx] == "loop_":
            idx += 1; tags = []
            while idx < len(tokens) and tokens[idx].startswith("_"):
                tags.append(tokens[idx]); idx += 1
            if not any("_atom_site." in t for t in tags): continue
            col = {t.split(".",1)[1]: ci for ci,t in enumerate(tags)}
            for k in ("label_atom_id","label_comp_id","auth_seq_id","Cartn_x","Cartn_y","Cartn_z",
                      "B_iso_or_equiv","pdbx_PDB_model_num","auth_asym_id"):
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
                    if chain and row[col["auth_asym_id"]] != chain: continue
                    try:
                        atoms.append({"res": int(row[col["auth_seq_id"]]),
                                      "aa": AA3_TO_1[row[col["label_comp_id"]]],
                                      "b": float(row[col["B_iso_or_equiv"]]),
                                      "x": float(row[col["Cartn_x"]]),
                                      "y": float(row[col["Cartn_y"]]),
                                      "z": float(row[col["Cartn_z"]])})
                    except (ValueError, KeyError): continue
            continue
        if tokens[idx].startswith("_") and idx+1 < len(tokens): idx += 2; continue
        idx += 1
    return atoms

def kabsch(mobile, ref):
    M = np.array(mobile); R = np.array(ref)
    if M.shape != R.shape or M.shape[0] < 3: return None, None
    mc = np.array([math.fsum(M[:,j]) for j in range(3)]) / M.shape[0]
    rc = np.array([math.fsum(R[:,j]) for j in range(3)]) / R.shape[0]
    H = (M - mc).T @ (R - rc)
    U, S, Vt = np.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False)
    det = np.linalg.det(U @ Vt)
    sup = (M - mc) @ (U @ np.diag([1,1,1 if det>0 else -1]) @ Vt) + rc
    return sup, np.linalg.norm(sup - R, axis=1)

def lddt(sup, ref, cutoff_sq=225.0, thresholds=(0.5,1.0,2.0,4.0)):
    n = ref.shape[0]
    if n < 2: return np.full(n, np.nan)
    diff = ref[:,None,:] - ref[None,:,:]
    dsq = np.sum(diff**2, axis=-1); np.fill_diagonal(dsq, np.inf)
    nb = dsq < cutoff_sq
    sd = np.sqrt(np.sum((sup[:,None,:]-sup[None,:,:])**2, axis=-1))
    rd = np.sqrt(dsq)
    pres = np.zeros((n,n,len(thresholds)), dtype=bool)
    for ti, t in enumerate(thresholds): pres[:,:,ti] = np.abs(sd - rd) < t
    pres &= nb[:,:,None]
    nc = nb.sum(axis=1).astype(float)
    return np.where(nc > 0, pres.sum(axis=2).sum(axis=1)/(nc*len(thresholds)), np.nan)

def bfactor_z(bfs):
    a = np.array(bfs); med = np.median(a)
    q25, q75 = np.percentile(a, [25,75]); iqr = q75-q25 or 1.0
    return (a - med)/iqr

# Load qualifying entries
qual = json.load(open("results/e421/qualifying_entries.json"))
print(f"entries: {len(qual)}")

# Process each entry WITHOUT sequence-equality check
results = {}
req = 0
for i, entry in enumerate(qual):
    pdb_id = entry["pdb_id"]
    acc = entry["afdb_accession"]
    chain = "A"
    row = {"entry": pdb_id, "accession": acc}
    
    try:
        # 1. PDBe mapping
        time.sleep(max(0, 0.25 + random.uniform(-0.1, 0.1)))
        r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}", headers=UA, timeout=30)
        req += 1
        if r.status_code != 200:
            results[pdb_id] = {**row, "status": "mapping_fail"}; continue
        mdata = r.json()
        uniprot = mdata.get(pdb_id, {}).get("UniProt", {})
        if acc not in uniprot:
            results[pdb_id] = {**row, "status": "mapping_acc_absent"}; continue
        mobj = uniprot[acc]
        res_map = {}
        for m in mobj.get("mappings", []):
            s = m["start"]["author_residue_number"] or m["start"]["residue_number"]
            e = m["end"]["author_residue_number"] or m["end"]["residue_number"]
            for off in range(e - s + 1):
                res_map[s + off] = m["unp_start"] + off
        
        # 2. AF model (from the API-provided cifUrl)
        time.sleep(max(0, 0.25 + random.uniform(-0.1, 0.1)))
        r2 = requests.get(entry["cif_url"], headers=UA, timeout=60)
        req += 1
        if r2.status_code != 200:
            results[pdb_id] = {**row, "status": "af_fetch_fail", "code": r2.status_code}; continue
        af_atoms = parse_cif_ca(r2.text)
        af_by_pos = {a["res"]: a for a in af_atoms}
        
        # 3. Reference cif
        time.sleep(max(0, 0.25 + random.uniform(-0.1, 0.1)))
        r3 = requests.get(f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb_id}.cif", headers=UA, timeout=60)
        req += 1
        if r3.status_code != 200:
            results[pdb_id] = {**row, "status": "ref_fetch_fail", "code": r3.status_code}; continue
        ref_atoms = parse_cif_ca(r3.text, chain=chain)
        
        # 4. Align via SIFTS mapping (NO sequence check — mapping handles it)
        ref_ca = []; af_ca = []; unp_pos_l = []; ref_bf_l = []
        for a in ref_atoms:
            unp = res_map.get(a["res"])
            if unp is None: continue
            af_a = af_by_pos.get(unp)
            if af_a is None: continue
            # amino-acid check: reference AA must match AF model AA at this position
            if a["aa"] != af_a["aa"]: continue
            ref_ca.append([a["x"], a["y"], a["z"]])
            af_ca.append([af_a["x"], af_a["y"], af_a["z"]])
            unp_pos_l.append(unp)
            ref_bf_l.append(a["b"])
        
        if len(ref_ca) < 30:
            results[pdb_id] = {**row, "status": "insufficient_aligned", "n_aligned": len(ref_ca)}; continue
        
        # 5. Kabsch + lDDT
        sup, dists = kabsch(af_ca, ref_ca)
        if sup is None: results[pdb_id] = {**row, "status": "kabsch_fail"}; continue
        lddt_vals = lddt(sup, np.array(ref_ca))
        z_scores = bfactor_z(ref_bf_l)
        
        # 6. Labels
        labels = []
        for j in range(len(ref_ca)):
            l = lddt_vals[j] if not np.isnan(lddt_vals[j]) else -1
            d = dists[j]
            if l >= 0.60 and d <= 4.0: e = "correct"
            elif l < 0.60 and d > 4.0: e = "error"
            else: e = "U"
            zj = z_scores[j]
            if zj >= 1.0: f = "flexible"
            elif zj <= -1.0: f = "ordered"
            else: f = "U"
            cell = "U"
            if e != "U" and f != "U":
                if e == "correct" and f == "ordered": cell = "OC"
                elif e == "correct" and f == "flexible": cell = "A"
                elif e == "error" and f == "ordered": cell = "B"
                elif e == "error" and f == "flexible": cell = "O"
            labels.append({"accession": acc, "entry": pdb_id,
                           "uniprot_pos": unp_pos_l[j], "e": e, "f": f, "cell": cell,
                           "lddt": round(float(lddt_vals[j]),4), "dist": round(float(dists[j]),4),
                           "b": ref_bf_l[j], "z": round(float(zj),4)})
        
        results[pdb_id] = {**row, "status": "ok", "labels": labels, "n_aligned": len(ref_ca),
                           "n_cells": dict(Counter(l["cell"] for l in labels))}
    except Exception as e:
        results[pdb_id] = {**row, "status": "exception", "detail": str(e)[:200]}
    
    if (i+1) % 25 == 0:
        ok_count = sum(1 for r in results.values() if r.get("status") == "ok")
        total_labels = sum(len(r.get("labels",[])) for r in results.values() if r.get("status") == "ok")
        print(f"  {i+1}/{len(qual)} | ok: {ok_count} | labels: {total_labels}")

# Summary
ok_entries = [r for r in results.values() if r.get("status") == "ok"]
all_labels = [l for r in ok_entries for l in r.get("labels",[])]
cells = Counter(l["cell"] for l in all_labels if l["cell"] != "U")
acc_with_labels = set(l["accession"] for l in all_labels if l["cell"] != "U")
print(f"\n{'='*60}")
print(f"CENSUS COMPLETE: {len(results)} entries processed")
print(f"OK entries: {len(ok_entries)} | labels: {len(all_labels)}")
print(f"Cells: {dict(cells)}")
print(f"Accessions with labels: {len(acc_with_labels)}")

verdict = "E421_H1_SUPPORT_PASS" if all(cells.get(c,0) >= 200 for c in ("A","B")) else "E421_UNDERPOWERED_POOL"
print(f"VERDICT: {verdict}")

json.dump(results, open("results/e421/census_results.json","w"), indent=1, default=str)
print("saved census_results.json")
