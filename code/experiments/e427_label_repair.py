"""e427 label-mapping repair (PROVISIONAL E-number, pending main's allocation).

Defect (confirmed 2026-10-07, paper/arxiv/verification/audit/AUDIT_REVISION2.md
F1): the frozen e422 runner and its e423 dispatch copy build the SIFTS
residue map from ``author_residue_number or residue_number``.  PDBe gives
``residue_number`` in mmCIF ``label_seq_id`` space, while the reference CA
atoms are keyed by ``auth_seq_id`` (e421_pipeline.parse_cif_ca).  When a
mapped segment starts on an unobserved residue, ``author_residue_number`` is
null, the map silently switches to label numbering, and every reference
residue is paired with the wrong UniProt position (9qj6 / K7PQ54: a 34-residue
register shift; true lDDT 0.99, ledger lDDT 0.20).  Nothing caught it because
the label path has no per-residue amino-acid identity check.

Repair (this module; the frozen e422_runner.py / e423_runner.py bytes are
NOT modified):
  1. ``parse_mappings_labelspace`` keys every residue map on ``residue_number``
     (label_seq_id) and binds the segment to its ``struct_asym_id``
     (label_asym_id).  A segment without a usable residue_number or
     struct_asym_id is the typed MAPPING_PARSE_FAILED — never a fallback to
     author numbering.  A segment whose label and UniProt spans differ
     (internal insertion/deletion) is not offset-mapped: its residue map is
     built at label time by global alignment of the reference entity
     sequence (_entity_poly_seq, label numbering) to the UniProt sequence
     carried by the AF model, keeping identical aligned pairs only.
  2. ``parse_cif_ca_labelspace`` reads reference CA atoms keyed by
     (label_asym_id, label_seq_id), with the same model / altloc / residue
     filters as e421_pipeline.parse_cif_ca.  The AF model keeps the e421 rule
     (auth_seq_id = UniProt position for AFDB F1 models).
  3. ``label_triple_labelspace`` adds a per-residue identity check: a mapped
     pair whose reference and AF residues differ is excluded and counted; if
     more than MAX_MISMATCH_FRACTION of mapped pairs differ, the whole triple
     is the typed E427_SEQUENCE_IDENTITY_FAILED (a register error, not a few
     engineered point mutations).

Output rows keep the e423 schema; triples additionally carry
``mapping_space`` and ``n_identity_mismatch``.
"""
from __future__ import annotations

import json
import os
import sys
from enum import Enum

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census  # noqa: E402
import e421_pipeline as pipe  # noqa: E402
from e423_runner import MIN_ALIGNED, RunnerFailureCode  # noqa: E402

MAPPING_SPACE = "label_seq_id"
# A register shift mismatches most residues; engineered point mutations
# mismatch a handful.  10% separates the two with a wide margin.
MAX_MISMATCH_FRACTION = 0.10


class E427FailureCode(str, Enum):
    SEQUENCE_IDENTITY_FAILED = "E427_SEQUENCE_IDENTITY_FAILED"
    REF_NO_LABEL_ASYM = "E427_REF_NO_LABEL_ASYM"


def _int_or_none(v):
    return v if isinstance(v, int) and not isinstance(v, bool) else None


def parse_mappings_labelspace(payload: bytes, entry: str) -> list[dict]:
    """PDBe whole-entry UniProt mapping -> triples keyed in label_seq_id space.

    Returns [{chain_id, struct_asym_id, entity_id, accession, res_map}] where
    res_map maps label_seq_id -> UniProt position.  Every structural surprise
    is the typed MAPPING_PARSE_FAILED (same envelope as e423).
    """
    try:
        doc = json.loads(payload.decode("utf-8"))
        uniprot = doc[entry.lower()]["UniProt"]
    except Exception as exc:  # noqa: BLE001 - typed
        raise census.CensusError(
            RunnerFailureCode.MAPPING_PARSE_FAILED.value,
            f"{entry}: {type(exc).__name__}") from exc
    if not isinstance(uniprot, dict):
        raise census.CensusError(
            RunnerFailureCode.MAPPING_PARSE_FAILED.value,
            f"{entry}: UniProt block is {type(uniprot).__name__}, not an object")
    triples: list[dict] = []
    for accession in sorted(uniprot):
        try:
            acc_block = uniprot[accession]
            mappings = (acc_block.get("mappings", [])
                        if isinstance(acc_block, dict) else None)
            if not isinstance(mappings, list):
                raise TypeError("accession value is not an object with a mappings list")
            for mapping in mappings:
                chain_id = mapping["chain_id"]
                asym = mapping.get("struct_asym_id")
                if not isinstance(asym, str) or not asym:
                    raise TypeError("missing struct_asym_id")
                unp_s = _int_or_none(mapping["unp_start"])
                unp_e = _int_or_none(mapping["unp_end"])
                s = _int_or_none(mapping["start"]["residue_number"])
                e = _int_or_none(mapping["end"]["residue_number"])
                if None in (unp_s, unp_e, s, e):
                    raise TypeError("non-integer label-space residue bounds")
                if e < s or unp_e < unp_s:
                    raise TypeError("inverted residue bounds")
                seg = {"label": (s, e), "unp": (unp_s, unp_e)}
                if (e - s) == (unp_e - unp_s):
                    # contiguous in both spaces: SIFTS offset map is exact
                    res_map = {s + off: unp_s + off for off in range(e - s + 1)}
                    method = "sifts_offset"
                else:
                    # internal insertion/deletion: an offset would guess;
                    # resolved per residue by sequence alignment at label time
                    res_map = None
                    method = "needs_alignment"
                triples.append({"chain_id": chain_id, "struct_asym_id": asym,
                                "entity_id": mapping["entity_id"],
                                "accession": accession, "res_map": res_map,
                                "segment": seg, "map_method": method})
        except (KeyError, TypeError) as exc:
            raise census.CensusError(
                RunnerFailureCode.MAPPING_PARSE_FAILED.value,
                f"{entry}/{accession}: malformed mapping item "
                f"({type(exc).__name__}: {exc})") from exc
    return triples


def parse_entity_poly_seq(text: str, entity_id) -> dict[int, str]:
    """label_seq_id -> one-letter residue for one entity (_entity_poly_seq)."""
    tokens = pipe.tokenize_cif(text)
    seq: dict[int, str] = {}
    idx = 0
    while idx < len(tokens):
        if tokens[idx] == "loop_":
            idx += 1
            tags = []
            while idx < len(tokens) and tokens[idx].startswith("_"):
                tags.append(tokens[idx])
                idx += 1
            if not any(t.startswith("_entity_poly_seq.") for t in tags):
                continue
            col = {t.split(".", 1)[1]: ci for ci, t in enumerate(tags)}
            while idx < len(tokens):
                t = tokens[idx]
                if t in ("loop_", "stop_") or t.startswith("_") or t.startswith("data_"):
                    break
                row = tokens[idx:idx + len(tags)]
                if len(row) < len(tags):
                    break
                idx += len(tags)
                if row[col["entity_id"]] != str(entity_id):
                    continue
                try:
                    num = int(row[col["num"]])
                except ValueError:
                    continue
                if num not in seq:  # first of any microheterogeneity
                    seq[num] = pipe.AA3_TO_1.get(row[col["mon_id"]], "X")
            continue
        idx += 1
    return seq


def align_map(ref_seq: list[tuple[int, str]], unp_seq: list[tuple[int, str]]) -> dict[int, int]:
    """Global Needleman-Wunsch (match +2, mismatch -1, gap -2) of two
    (number, residue) lists; returns ref number -> unp number for aligned
    IDENTICAL pairs only."""
    n, m = len(ref_seq), len(unp_seq)
    if n == 0 or m == 0:
        return {}
    gap = -2
    sc = np.zeros((n + 1, m + 1))
    sc[:, 0] = gap * np.arange(n + 1)
    sc[0, :] = gap * np.arange(m + 1)
    tb = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left
    tb[1:, 0] = 1
    tb[0, 1:] = 2
    ua = np.array([c for _, c in unp_seq])
    for i in range(1, n + 1):
        match = np.where(ua == ref_seq[i - 1][1], 2.0, -1.0)
        diag = sc[i - 1, :-1] + match
        up = sc[i - 1, 1:] + gap
        best = np.maximum(diag, up)
        choice = np.where(diag >= up, 0, 1)
        row = sc[i]
        for j in range(1, m + 1):  # left moves are sequential
            left = row[j - 1] + gap
            if best[j - 1] >= left:
                row[j] = best[j - 1]
                tb[i, j] = choice[j - 1]
            else:
                row[j] = left
                tb[i, j] = 2
    out: dict[int, int] = {}
    i, j = n, m
    while i > 0 and j > 0:
        t = tb[i, j]
        if t == 0:
            if ref_seq[i - 1][1] == unp_seq[j - 1][1]:
                out[ref_seq[i - 1][0]] = unp_seq[j - 1][0]
            i, j = i - 1, j - 1
        elif t == 1:
            i -= 1
        else:
            j -= 1
    return out


def parse_cif_ca_labelspace(text: str, label_asym: str) -> list[dict]:
    """Reference CA atoms of one label_asym_id, keyed by label_seq_id."""
    tokens = pipe.tokenize_cif(text)
    need = ("label_atom_id", "label_comp_id", "label_seq_id", "label_asym_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv",
            "pdbx_PDB_model_num", "label_alt_id")
    atoms: list[dict] = []
    idx = 0
    while idx < len(tokens):
        if tokens[idx] == "loop_":
            idx += 1
            tags = []
            while idx < len(tokens) and tokens[idx].startswith("_"):
                tags.append(tokens[idx])
                idx += 1
            if not any(t.startswith("_atom_site.") for t in tags):
                continue
            col = {t.split(".", 1)[1]: ci for ci, t in enumerate(tags)}
            if any(k not in col for k in need):
                continue
            while idx < len(tokens):
                t = tokens[idx]
                if t in ("loop_", "stop_") or t.startswith("_") or t.startswith("data_"):
                    break
                row = tokens[idx:idx + len(tags)]
                if len(row) < len(tags):
                    break
                idx += len(tags)
                if row[col["label_atom_id"]] != "CA":
                    continue
                if row[col["label_comp_id"]] not in pipe.AA3_TO_1:
                    continue
                if row[col["pdbx_PDB_model_num"]] != "1":
                    continue
                if row[col["label_alt_id"]] not in (".", "?", "", "A"):
                    continue
                if row[col["label_asym_id"]] != label_asym:
                    continue
                try:
                    atoms.append({
                        "label_seq_id": int(row[col["label_seq_id"]]),
                        "x": float(row[col["Cartn_x"]]),
                        "y": float(row[col["Cartn_y"]]),
                        "z": float(row[col["Cartn_z"]]),
                        "b_factor": float(row[col["B_iso_or_equiv"]]),
                        "aa": pipe.AA3_TO_1[row[col["label_comp_id"]]],
                    })
                except (ValueError, KeyError):
                    continue
            continue
        idx += 1
    return atoms


def label_triple_labelspace(triple: dict, entry: str, af_cif_text: str,
                            ref_cif_text: str, af_sha: str, ref_sha: str) -> dict:
    """e423.label_triple with label-space mapping and an identity check."""
    out = {"entry": entry, "chain": triple["chain_id"],
           "accession": triple["accession"], "rows": [], "status": None,
           "af_sha256": af_sha, "ref_sha256": ref_sha,
           "mapping_space": MAPPING_SPACE, "n_identity_mismatch": 0}
    af_atoms = pipe.parse_cif_ca(af_cif_text)
    if len(af_atoms) < 3:
        out["status"] = RunnerFailureCode.AF_INSUFFICIENT_CA.value
        return out
    af_by_pos = {a["residue_number"]: a for a in af_atoms}
    asym = triple.get("struct_asym_id")
    if not asym:
        out["status"] = E427FailureCode.REF_NO_LABEL_ASYM.value
        return out
    ref_atoms = parse_cif_ca_labelspace(ref_cif_text, asym)
    if len(ref_atoms) < MIN_ALIGNED:
        out["status"] = RunnerFailureCode.REF_INSUFFICIENT_CA.value
        return out
    res_map = triple["res_map"]
    out["map_method"] = triple.get("map_method", "sifts_offset")
    if res_map is None:
        (ls, le), (us, ue) = triple["segment"]["label"], triple["segment"]["unp"]
        ent = parse_entity_poly_seq(ref_cif_text, triple["entity_id"])
        ref_seq = [(k, ent[k]) for k in range(ls, le + 1) if k in ent]
        unp_seq = [(p, af_by_pos[p]["aa"]) for p in range(us, ue + 1) if p in af_by_pos]
        res_map = align_map(ref_seq, unp_seq)
    ref_ca, af_ca, unp_positions, ref_idx = [], [], [], []
    seen_pos: set[int] = set()
    n_pairs = n_mismatch = 0
    for i, atom in enumerate(ref_atoms):
        unp = res_map.get(atom["label_seq_id"])
        if unp is None or unp in seen_pos:
            continue
        af_atom = af_by_pos.get(unp)
        if af_atom is None:
            continue
        n_pairs += 1
        if af_atom["aa"] != atom["aa"]:
            n_mismatch += 1
            continue
        seen_pos.add(unp)
        ref_ca.append([atom["x"], atom["y"], atom["z"]])
        af_ca.append([af_atom["x"], af_atom["y"], af_atom["z"]])
        unp_positions.append(unp)
        ref_idx.append(i)
    out["n_identity_mismatch"] = n_mismatch
    if n_pairs and n_mismatch / n_pairs > MAX_MISMATCH_FRACTION:
        out["status"] = E427FailureCode.SEQUENCE_IDENTITY_FAILED.value
        return out
    if len(ref_ca) < MIN_ALIGNED:
        out["status"] = RunnerFailureCode.INSUFFICIENT_ALIGNED.value
        return out
    sup, dists = pipe.kabsch_superpose(af_ca, ref_ca)
    if sup is None:
        out["status"] = RunnerFailureCode.KABSCH_FAILED.value
        return out
    lddt, _nbr = pipe.lddt_calc(sup, np.array(ref_ca))
    # F channel: robust z over ALL CA b-factors of the mapped chain (e421
    # rule), indexed back to the aligned subset via ref_idx.
    z, _med, _iqr = pipe.bfactor_z([a["b_factor"] for a in ref_atoms])
    af_bf = {a["residue_number"]: a["b_factor"] for a in af_atoms}
    rows = []
    for j, unp in enumerate(unp_positions):
        plddt = af_bf.get(unp)
        lv = None if np.isnan(lddt[j]) else float(lddt[j])
        rows.append({
            "entry": entry, "chain": triple["chain_id"],
            "accession": triple["accession"], "uniprot_pos": int(unp),
            "plddt": round(float(plddt), 4) if plddt is not None else None,
            "lddt": round(lv, 4) if lv is not None else None,
            "d_kabsch": round(float(dists[j]), 4),
            "b_factor_z": round(float(z[ref_idx[j]]), 4),
        })
    out["rows"] = rows
    out["status"] = "ok"
    return out
