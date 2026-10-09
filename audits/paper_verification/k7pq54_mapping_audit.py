#!/usr/bin/env python3
"""Independent SIFTS-to-coordinate audit for 9qj6 chains A/B.

This script uses only the cached PDBe mmCIF, PDBe UniProt/SIFTS mapping JSON,
the cached AFDB v6 mmCIF, the Python standard library, and NumPy. It does not
import project code. Run from outside the repository with ``python -I``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shlex
from pathlib import Path
from statistics import mean, median

import numpy as np

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def cif_ca(path: Path) -> dict[str, dict[int, dict]]:
    """Read CA records from atom_site loops, keyed by auth chain and label seq."""
    lines = path.read_text(errors="replace").splitlines()
    result: dict[str, dict[int, dict]] = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        tags = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            tags.append(lines[i].strip().split()[0])
            i += 1
        if not any(tag.startswith("_atom_site.") for tag in tags):
            continue
        col = {tag.split(".", 1)[1]: j for j, tag in enumerate(tags)}
        required = {
            "label_atom_id", "label_comp_id", "label_alt_id", "label_seq_id",
            "auth_seq_id", "auth_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
            "B_iso_or_equiv", "pdbx_PDB_model_num",
        }
        missing = sorted(required - set(col))
        if missing:
            raise ValueError(f"{path}: atom_site loop missing {missing}")
        while i < len(lines):
            line = lines[i].strip()
            if not line or line == "#":
                i += 1
                continue
            if line == "loop_" or line.startswith("_") or line.startswith("data_"):
                break
            fields = shlex.split(line, posix=True)
            i += 1
            if len(fields) != len(tags):
                raise ValueError(
                    f"{path}: atom_site row has {len(fields)} values, expected {len(tags)}"
                )
            if fields[col["label_atom_id"]] != "CA":
                continue
            comp = fields[col["label_comp_id"]]
            if comp not in AA3 or fields[col["pdbx_PDB_model_num"]] != "1":
                continue
            if fields[col["label_alt_id"]] not in (".", "?", "", "A"):
                continue
            try:
                label_seq = int(fields[col["label_seq_id"]])
                auth_seq = int(fields[col["auth_seq_id"]])
                atom = {
                    "label_seq": label_seq,
                    "auth_seq": auth_seq,
                    "aa": AA3[comp],
                    "xyz": [float(fields[col[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                    "bfactor": float(fields[col["B_iso_or_equiv"]]),
                    "alt": fields[col["label_alt_id"]],
                }
            except (ValueError, KeyError):
                continue
            chain = fields[col["auth_asym_id"]]
            existing = result.setdefault(chain, {}).get(label_seq)
            # Prefer a nonalternate record, then alt A. In this structure the
            # CA atom is unique at each chain/label position.
            rank = {".": 0, "?": 0, "": 0, "A": 1}
            if existing is None or rank[atom["alt"]] < rank[existing["alt"]]:
                result[chain][label_seq] = atom
    return result


def af_ca(path: Path) -> dict[int, dict]:
    """Read AFDB C-alpha coordinates indexed by the embedded UniProt index."""
    lines = path.read_text(errors="replace").splitlines()
    result = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        tags = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            tags.append(lines[i].strip().split()[0])
            i += 1
        if not any(tag.startswith("_atom_site.") for tag in tags):
            continue
        col = {tag.split(".", 1)[1]: j for j, tag in enumerate(tags)}
        required = {
            "label_atom_id", "label_comp_id", "label_alt_id", "label_seq_id",
            "auth_seq_id", "auth_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
            "B_iso_or_equiv", "pdbx_PDB_model_num", "pdbx_sifts_xref_db_num",
            "pdbx_sifts_xref_db_acc", "pdbx_sifts_xref_db_name",
        }
        missing = sorted(required - set(col))
        if missing:
            raise ValueError(f"{path}: atom_site loop missing {missing}")
        while i < len(lines):
            line = lines[i].strip()
            if not line or line == "#":
                i += 1
                continue
            if line == "loop_" or line.startswith("_") or line.startswith("data_"):
                break
            fields = shlex.split(line, posix=True)
            i += 1
            if len(fields) != len(tags):
                raise ValueError(
                    f"{path}: atom_site row has {len(fields)} values, expected {len(tags)}"
                )
            if fields[col["label_atom_id"]] != "CA":
                continue
            comp = fields[col["label_comp_id"]]
            if comp not in AA3 or fields[col["pdbx_PDB_model_num"]] != "1":
                continue
            if fields[col["label_alt_id"]] not in (".", "?", "", "A"):
                continue
            if (fields[col["pdbx_sifts_xref_db_acc"]] != "K7PQ54"
                    or fields[col["pdbx_sifts_xref_db_name"]] != "UNP"):
                continue
            try:
                unp = int(fields[col["pdbx_sifts_xref_db_num"]])
                atom = {
                    "uniprot_pos": unp,
                    "auth_seq": int(fields[col["auth_seq_id"]]),
                    "aa": AA3[comp],
                    "xyz": [float(fields[col[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                    "bfactor": float(fields[col["B_iso_or_equiv"]]),
                }
            except (ValueError, KeyError):
                continue
            result[unp] = atom
    return result


def kabsch_distances(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Least-squares AF-to-reference fit, with a proper (nonreflecting) rotation."""
    x0 = mobile.mean(axis=0)
    y0 = reference.mean(axis=0)
    x = mobile - x0
    y = reference - y0
    u, _s, vt = np.linalg.svd(x.T @ y, full_matrices=False)
    fix = np.eye(3)
    fix[2, 2] = 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0
    rotation = u @ fix @ vt
    fitted = x @ rotation + y0
    return np.linalg.norm(fitted - reference, axis=1)


def local_lddt(af: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-position lDDT-Ca; neighbors selected by reference distance < 15 A."""
    d_ref = np.linalg.norm(ref[:, None, :] - ref[None, :, :], axis=-1)
    d_af = np.linalg.norm(af[:, None, :] - af[None, :, :], axis=-1)
    np.fill_diagonal(d_ref, np.inf)
    neighbor = d_ref < 15.0
    delta = np.abs(d_af - d_ref)
    preserved = np.stack([delta < t for t in (0.5, 1.0, 2.0, 4.0)], axis=2)
    counts = neighbor.sum(axis=1)
    values = np.full(len(ref), np.nan)
    valid = counts > 0
    values[valid] = (preserved & neighbor[:, :, None]).sum(axis=(1, 2))[valid] / (4.0 * counts[valid])
    return values, counts


def aligned_metrics(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    rows.sort(key=lambda row: row["uniprot_pos"])
    ref = np.asarray([r["pdb_xyz"] for r in rows], dtype=float)
    af = np.asarray([r["af_xyz"] for r in rows], dtype=float)
    return kabsch_distances(af, ref), local_lddt(af, ref)[0]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--sifts", type=Path, required=True)
    parser.add_argument("--af", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, nargs="*", default=[])
    args = parser.parse_args()

    pdb_atoms = cif_ca(args.pdb)
    af_atoms = af_ca(args.af)
    sifts_doc = json.loads(args.sifts.read_text())
    maps = sifts_doc["9qj6"]["UniProt"]["K7PQ54"]["mappings"]
    map_by_chain = {m["chain_id"]: m for m in maps}
    all_chains = {}
    checkpoint_matches = []
    mapping_reproduction = {}
    legacy_maps = {}
    for chain in ("A", "B"):
        mapping = map_by_chain[chain]
        start = mapping["start"]
        use_author_number = start.get("author_residue_number") is not None
        cif_space = "auth_seq" if use_author_number else "label_seq"
        pdb_chain = pdb_atoms.get(chain, {})
        matched = []
        seen_unp = set()
        for atom in pdb_chain.values():
            pdb_number = atom[cif_space]
            map_start = (start["author_residue_number"] if use_author_number
                         else start["residue_number"])
            unp_pos = int(mapping["unp_start"]) + (pdb_number - int(map_start))
            if not (int(mapping["unp_start"]) <= unp_pos <= int(mapping["unp_end"])):
                continue
            if unp_pos in seen_unp or unp_pos not in af_atoms:
                continue
            seen_unp.add(unp_pos)
            ref_atom = atom
            af_atom = af_atoms[unp_pos]
            matched.append({
                "uniprot_pos": unp_pos,
                "pdb_label_seq": ref_atom["label_seq"],
                "pdb_auth_seq": ref_atom["auth_seq"],
                "aa_pdb": ref_atom["aa"],
                "aa_af": af_atom["aa"],
                "pdb_xyz": ref_atom["xyz"],
                "af_xyz": af_atom["xyz"],
                "plddt": af_atom["bfactor"],
            })
        matched.sort(key=lambda x: x["uniprot_pos"])
        d_kabsch, lddt = aligned_metrics(matched)
        _lddt_again, neighbor_count = local_lddt(
            np.asarray([r["af_xyz"] for r in matched], dtype=float),
            np.asarray([r["pdb_xyz"] for r in matched], dtype=float),
        )
        per_residue = []
        for row, dist, score, nbr in zip(matched, d_kabsch, lddt, neighbor_count):
            per_residue.append({
                "uniprot_pos": row["uniprot_pos"],
                "pdb_label_seq": row["pdb_label_seq"],
                "pdb_auth_seq": row["pdb_auth_seq"],
                "aa": row["aa_pdb"],
                "aa_matches_af": row["aa_pdb"] == row["aa_af"],
                "plddt": round(float(row["plddt"]), 4),
                "d_kabsch_A": round(float(dist), 4),
                "lddt_ca": None if math.isnan(float(score)) else round(float(score), 6),
                "reference_neighbors_lt15A": int(nbr),
            })
        low_group = [r for r in per_residue if r["uniprot_pos"] < 70]
        summary = {
            "chain": chain,
            "mapping_source": "PDBe SIFTS-derived UniProt mapping API response",
            "mapping_coordinate_space": cif_space,
            "sifts": mapping,
            "n_ca_reference": len(pdb_chain),
            "n_aligned": len(per_residue),
            "uniprot_range_aligned": [per_residue[0]["uniprot_pos"], per_residue[-1]["uniprot_pos"]],
            "n_sequence_mismatches": sum(not r["aa_matches_af"] for r in per_residue),
            "median_plddt": round(median(r["plddt"] for r in per_residue), 4),
            "median_d_kabsch_A": round(median(r["d_kabsch_A"] for r in per_residue), 4),
            "mean_d_kabsch_A": round(mean(r["d_kabsch_A"] for r in per_residue), 4),
            "median_lddt_ca": round(median(r["lddt_ca"] for r in per_residue if r["lddt_ca"] is not None), 6),
            "mean_lddt_ca": round(mean(r["lddt_ca"] for r in per_residue if r["lddt_ca"] is not None), 6),
            "first_10": per_residue[:10],
            "last_10": per_residue[-10:],
            "residues_before_70": {
                "n": len(low_group),
                "median_plddt": round(median(r["plddt"] for r in low_group), 4) if low_group else None,
                "median_d_kabsch_A": round(median(r["d_kabsch_A"] for r in low_group), 4) if low_group else None,
                "median_lddt_ca": round(median(r["lddt_ca"] for r in low_group), 6) if low_group else None,
            },
            "per_residue": per_residue,
        }
        all_chains[chain] = summary

        # Recreate the runner's number-space error independently: it builds
        # a mapping from start.residue_number (label numbering here) but then
        # applies that map to auth_seq_id. This intentionally wrong alignment
        # should reproduce the stored checkpoint rows if it is causal.
        legacy = []
        map_start = int(mapping["start"]["residue_number"])
        for atom in pdb_chain.values():
            legacy_unp = int(mapping["unp_start"]) + (atom["auth_seq"] - map_start)
            if not (int(mapping["unp_start"]) <= legacy_unp <= int(mapping["unp_end"])):
                continue
            if legacy_unp not in af_atoms:
                continue
            legacy.append({
                "uniprot_pos": legacy_unp,
                "pdb_xyz": atom["xyz"],
                "af_xyz": af_atoms[legacy_unp]["xyz"],
                "plddt": af_atoms[legacy_unp]["bfactor"],
            })
        legacy_d, legacy_l = aligned_metrics(legacy)
        legacy_by_pos = {
            row["uniprot_pos"]: {
                "plddt": row["plddt"],
                "d_kabsch": float(dist),
                "lddt": None if math.isnan(float(score)) else float(score),
            }
            for row, dist, score in zip(legacy, legacy_d, legacy_l)
        }
        legacy_maps[chain] = legacy_by_pos
        mapping_reproduction[chain] = {
            "legacy_alignment_n": len(legacy),
            "legacy_uniprot_range": [min(legacy_by_pos), max(legacy_by_pos)],
            "legacy_median_plddt": round(median(row["plddt"] for row in legacy_by_pos.values()), 4),
            "legacy_median_d_kabsch_A": round(median(row["d_kabsch"] for row in legacy_by_pos.values()), 4),
            "legacy_median_lddt_ca": round(median(row["lddt"] for row in legacy_by_pos.values()), 6),
        }

    for path in args.checkpoints:
        checkpoint = json.loads(path.read_text())
        for row in checkpoint.get("rows", []):
            if row.get("accession") == "K7PQ54" and row.get("entry", "").lower() == "9qj6":
                checkpoint_matches.append({
                    "checkpoint": str(path), "chain": row.get("chain"),
                    "uniprot_pos": row.get("uniprot_pos"), "plddt": row.get("plddt"),
                    "d_kabsch": row.get("d_kabsch"), "lddt": row.get("lddt"),
                })

    e422_summary = {}
    for chain in ("A", "B"):
        entries = [r for r in checkpoint_matches if r["chain"] == chain]
        e422_summary[chain] = {
            "n_rows": len(entries),
            "uniprot_range": [min(r["uniprot_pos"] for r in entries), max(r["uniprot_pos"] for r in entries)] if entries else None,
            "median_plddt": round(median(r["plddt"] for r in entries), 4) if entries else None,
            "median_d_kabsch_A": round(median(r["d_kabsch"] for r in entries), 4) if entries else None,
            "median_lddt_ca": round(median(r["lddt"] for r in entries), 6) if entries else None,
            "first_rows": entries[:3],
        }
        stored_by_pos = {r["uniprot_pos"]: r for r in entries}
        differences = []
        for pos, calculated in legacy_maps[chain].items():
            observed = stored_by_pos.get(pos)
            if observed is None:
                continue
            differences.extend([
                abs(calculated["plddt"] - observed["plddt"]),
                abs(calculated["d_kabsch"] - observed["d_kabsch"]),
                abs(calculated["lddt"] - observed["lddt"]),
            ])
        mapping_reproduction[chain]["checkpoint_common_positions"] = len(set(stored_by_pos) & set(legacy_maps[chain]))
        mapping_reproduction[chain]["checkpoint_max_abs_raw_difference"] = round(max(differences), 6) if differences else None
    output = {
        "method": {
            "mapping": "PDBe UniProt mapping start/end with label_seq_id because SIFTS API has null author_residue_number",
            "af_index": "AFDB atom_site.pdbx_sifts_xref_db_num (UniProt residue number)",
            "distance": "per-residue C-alpha Euclidean residual after least-squares Kabsch AF-to-PDB fit per chain",
            "lddt": "pairwise C-alpha lDDT; native reference distance <15 A; strict delta thresholds <0.5, <1, <2, <4 A",
            "numpy_version": np.__version__,
        },
        "inputs": {
            "pdbe_cif": {"path": str(args.pdb), "sha256": digest(args.pdb)},
            "pdbe_sifts_json": {"path": str(args.sifts), "sha256": digest(args.sifts)},
            "alphafold_cif": {"path": str(args.af), "sha256": digest(args.af)},
        },
        "chains": all_chains,
        "runner_mapping_error_reproduction": mapping_reproduction,
        "e422_checkpoint_matches_summary": e422_summary,
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        chain: {k: v for k, v in summary.items() if k != "per_residue" and k not in ("first_10", "last_10")}
        for chain, summary in all_chains.items()
    }, indent=2))
    print("e422 9qj6 checkpoint comparison:", json.dumps(e422_summary))
    print("output:", args.output)


if __name__ == "__main__":
    main()
