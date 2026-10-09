"""The one corrected chain segment the biotite cross-check does not reproduce: 21xg chain A (Q4KH59).

Recomputes chain A's per-residue lDDT with the biotite pipeline of
full_ledger_biotite.py (residue-level SIFTS from the PDBe updated mmCIF) and
compares it with the corrected ledger: positions present on one side only,
the number of lDDT values that differ by more than 1e-4 and the largest
difference, the chain medians, and the number of residues below lDDT 0.60
on each side. It also lists, for label_seq_id 1-40, the entity sequence,
the observed reference residue, the SIFTS UniProt position and the AFDB
residue at that position, which shows the alignment tie next to the two
crosslinked residues. No pLDDT value is read or written.

Usage: python -I analyse_21xg.py --derived data/derived --raw _work/results/e422/raw \
           --updated _updated_mmcif --out analyse_21xg_out.json
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import biotite.structure.io.pdbx as px

sys.path.insert(0, str(Path(__file__).resolve().parent))
import full_ledger_biotite as fl  # noqa: E402  (helpers of this cross-check, not study code)

ENTRY, CHAIN, ACC = "21xg", "A", "Q4KH59"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--derived", default="data/derived")
    ap.add_argument("--raw", required=True)
    ap.add_argument("--updated", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rc, ledger, _ = fl.receipts_and_ledger(a.derived)
    upd_path = Path(a.updated) / f"{ENTRY}_updated.cif"
    ref = fl.table(fl.raw_file(a.raw, rc, f"https://www.ebi.ac.uk/pdbe/entry-files/{ENTRY}.cif"))
    upd = fl.table(upd_path)
    afd = fl.af_ca(a.raw, rc, ACC)
    mp = json.loads(fl.raw_file(a.raw, rc, f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{ENTRY}").read_text())
    segs = [m for m in mp[ENTRY]["UniProt"][ACC]["mappings"] if str(m["chain_id"]) == CHAIN]
    sifts = fl.sifts_residues(upd, CHAIN, ACC)
    pairs, mism = fl.pair_residues(ref, fl.reference_ca(ref, CHAIN), sifts, afd)
    bt = fl.segment_lddt(pairs, segs, afd)
    rows = ledger[(ENTRY, CHAIN, ACC)]
    led = {r["uniprot_pos"]: r["lddt"] for r in rows}
    common = sorted(set(led) & set(bt))
    diff = [(p, led[p], round(bt[p], 4), round(bt[p] - led[p], 4)) for p in common if abs(bt[p] - led[p]) > 1.0001e-4]

    # alignment listing (entity sequence, observed residue, SIFTS UniProt number, AFDB residue there)
    blk = px.CIFFile.read(str(upd_path)).block
    ps = blk["entity_poly_seq"]
    ent = {int(n): c for n, c, e in zip(ps["num"].as_array(), ps["mon_id"].as_array(), ps["entity_id"].as_array())
           if e == "1"}
    m = (ref["label_atom_id"] == "CA") & (ref["label_asym_id"] == "A")
    obs = {int(l): c for l, c in zip(ref["label_seq_id"][m], ref["label_comp_id"][m])}
    mu = (upd["label_atom_id"] == "CA") & (upd["label_asym_id"] == "A")
    sif = {int(l): n for l, n in zip(upd["label_seq_id"][mu], upd["pdbx_sifts_xref_db_num"][mu])}
    listing = []
    for lab in range(1, 41):
        n = sif.get(lab)
        unp = int(n) if n not in (None, ".", "?") else None
        listing.append({"label_seq_id": lab, "entity_residue": ent.get(lab), "observed_reference_residue": obs.get(lab),
                        "sifts_uniprot_pos": unp, "afdb_residue_at_sifts_pos": afd[unp][0] if unp in afd else None,
                        "in_corrected_ledger_at_sifts_pos": unp in led if unp is not None else None})
    out = {
        "target": f"{ENTRY}:{CHAIN}:{ACC}",
        "biotite_pairs": len(pairs), "biotite_identity_mismatches": mism, "ledger_rows": len(rows),
        "positions_only_in_ledger": {p: led[p] for p in sorted(set(led) - set(bt))},
        "positions_only_in_biotite": {p: round(bt[p], 4) for p in sorted(set(bt) - set(led))},
        "n_common_positions": len(common),
        "n_lddt_values_differing_gt_1e-4": len(diff),
        "max_abs_lddt_difference": max(abs(x[3]) for x in diff) if diff else 0.0,
        "max_abs_lddt_difference_unrounded": max((abs(bt[p] - led[p]) for p in common), default=0.0),
        "median_lddt_ledger": float(np.median(list(led.values()))),
        "median_lddt_biotite": float(np.median(list(bt.values()))),
        "median_difference": float(np.median(list(led.values())) - np.median(list(bt.values()))),
        "n_below_0.60_ledger": sum(v < 0.60 for v in led.values()),
        "n_below_0.60_biotite": sum(v < 0.60 for v in bt.values()),
        "crossing_0.60": {"ledger_below_biotite_not": [p for p in common if led[p] < 0.60 <= bt[p]],
                          "biotite_below_ledger_not": [p for p in common if bt[p] < 0.60 <= led[p]]},
        "differing_values": [{"uniprot_pos": p, "ledger": l, "biotite": b, "delta": d} for p, l, b, d in diff],
        "alignment_label_seq_1_40": listing,
    }
    Path(a.out).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("differing_values", "alignment_label_seq_1_40")},
                     indent=1))


if __name__ == "__main__":
    main()
