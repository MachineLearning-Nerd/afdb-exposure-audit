"""Observed residues the corrected labeller excludes by policy: alternative locations other than A, and
modified residues other than selenomethionine.

For every (entry, auth chain, accession) present in the corrected census
ledger (results/e427/batches, batches 2-4), reads the PDBe updated mmCIF and
counts observed polymer residues (model 1, a CA atom, a label_seq_id and a
SIFTS UniProt mapping to that accession) that the labeller's policy drops:
  * modified: no CA record with a standard amino acid or MSE (counted by
    component, e.g. MLY, FME, TPO);
  * altloc-not-A: standard residue whose CA atoms all carry an alternative
    location other than '.', '?', '' or 'A'.
No coordinates, lDDT or pLDDT values are read.

Usage: python -I nonstd_residues.py --derived data/derived --updated _updated_mmcif --out nonstd_out.json
"""
import argparse
import collections
import json
from pathlib import Path

import numpy as np
import biotite.structure.io.pdbx as px

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split())
COLS = ("label_atom_id", "auth_asym_id", "pdbx_sifts_xref_db_acc", "pdbx_PDB_model_num", "label_comp_id",
        "label_alt_id", "label_asym_id", "label_seq_id")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data/derived")
    ap.add_argument("--updated", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    groups = collections.defaultdict(set)
    for b in (2, 3, 4):
        ck = json.loads((Path(a.derived) / f"results/e427/batches/e422_batch{b:02d}_checkpoint.json").read_text())
        for r in ck["rows"]:
            groups[r["entry"].lower()].add((r["chain"], r["accession"]))
    comp, comp_entries, alt_only, alt_entries, total = collections.Counter(), set(), 0, set(), 0
    per_entry = {}
    for e, keys in sorted(groups.items()):
        site = px.CIFFile.read(str(Path(a.updated) / f"{e}_updated.cif")).block["atom_site"]
        t = {k: site[k].as_array() for k in COLS}
        m = (t["label_atom_id"] == "CA") & (t["pdbx_PDB_model_num"] == "1")
        res = collections.defaultdict(list)
        for i in np.nonzero(m)[0]:
            if (t["auth_asym_id"][i], t["pdbx_sifts_xref_db_acc"][i]) in keys and t["label_seq_id"][i] not in ".?":
                res[(t["label_asym_id"][i], t["label_seq_id"][i])].append((t["label_comp_id"][i], t["label_alt_id"][i]))
        n_mod = n_alt = 0
        for v in res.values():
            total += 1
            if any(c in AA and al in (".", "?", "", "A") for c, al in v):
                continue
            if any(c not in AA for c, _ in v):
                comp[v[0][0]] += 1
                comp_entries.add(e)
                n_mod += 1
            else:
                alt_only += 1
                alt_entries.add(e)
                n_alt += 1
        if n_mod or n_alt:
            per_entry[e] = {"modified": n_mod, "altloc_not_A": n_alt}
    out = {"observed_sifts_mapped_ca_residues": total,
           "excluded_modified_residues": sum(comp.values()), "modified_by_component": dict(comp.most_common()),
           "entries_with_modified": len(comp_entries),
           "excluded_altloc_not_A_residues": alt_only, "entries_with_altloc_not_A": len(alt_entries),
           "per_entry": per_entry}
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "per_entry"}, indent=1))


if __name__ == "__main__":
    main()
