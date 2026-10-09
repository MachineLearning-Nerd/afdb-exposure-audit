#!/usr/bin/env python3
"""Independent, no-project-import audit of corrected e427 checkpoints.

Run from outside the repository with the repository Python in isolated mode.
The only output is a JSON artifact in this audit directory.
"""
from __future__ import annotations

import copy
import argparse
import hashlib
import json
import random
import shlex
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RAW = ROOT / "results/e422/raw"
ORIG = ROOT / "results/e422/batches"
CORR = ROOT / "results/e427/batches"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canon_digest(doc: dict) -> str:
    return sha(json.dumps({k: v for k, v in doc.items() if k != "sha256"},
                          sort_keys=True, separators=(",", ":")).encode())


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def cif_tokens(text: str):
    """Small independent mmCIF tokenizer for quoted and semicolon values."""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c == "#":
            end = text.find("\n", i)
            i = n if end < 0 else end + 1
        elif c == ";" and (i == 0 or text[i - 1] == "\n"):
            end = text.find("\n;", i)
            if end < 0:
                raise ValueError("unterminated mmCIF semicolon text")
            yield text[i + 1:end]
            i = end + 2
        elif c in "'\"":
            q, end = c, i + 1
            while True:
                end = text.find(q, end)
                if end < 0:
                    raise ValueError("unterminated mmCIF quoted value")
                if end + 1 == n or text[end + 1].isspace():
                    yield text[i + 1:end]
                    i = end + 1
                    break
                end += 1
        else:
            end = i
            while end < n and not text[end].isspace():
                end += 1
            yield text[i:end]
            i = end


def cif_loops(text: str):
    toks = list(cif_tokens(text))
    i = 0
    while i < len(toks):
        if toks[i] != "loop_":
            i += 1
            continue
        i += 1
        tags = []
        while i < len(toks) and toks[i].startswith("_"):
            tags.append(toks[i])
            i += 1
        if not tags:
            continue
        width = len(tags)
        while i < len(toks):
            t = toks[i]
            if t in ("loop_", "stop_") or t.startswith("_") or t.startswith("data_"):
                break
            row = toks[i:i + width]
            if len(row) != width:
                raise ValueError("incomplete mmCIF loop row")
            yield tags, row
            i += width


def parse_ref_full(text: str, asym_id: str, entity_id):
    atoms = []
    sequence = {}
    variants = defaultdict(set)
    for tags, row in cif_loops(text):
        names = [t.split(".", 1)[1] for t in tags]
        if any(t.startswith("_entity_poly_seq.") for t in tags):
            c = {name: i for i, name in enumerate(names)}
            if not {"entity_id", "num", "mon_id"} <= set(c):
                continue
            try:
                if row[c["entity_id"]] == str(entity_id):
                    pos = int(row[c["num"]])
                    aa = AA3.get(row[c["mon_id"]], "X")
                    variants[pos].add((row[c["mon_id"]], aa))
                    sequence.setdefault(pos, aa)
            except (ValueError, KeyError):
                pass
            continue
        if not any(t.startswith("_atom_site.") for t in tags):
            continue
        c = {name: i for i, name in enumerate(names)}
        needed = {"label_atom_id", "label_comp_id", "label_alt_id", "label_seq_id",
                  "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
                  "B_iso_or_equiv", "pdbx_PDB_model_num"}
        if not needed <= set(c):
            continue
        if row[c["label_atom_id"]] != "CA" or row[c["label_asym_id"]] != asym_id:
            continue
        comp, alt = row[c["label_comp_id"]], row[c["label_alt_id"]]
        if comp not in AA3 or row[c["pdbx_PDB_model_num"]] != "1" or alt not in (".", "?", "", "A"):
            continue
        try:
            atoms.append({
                "label_seq": int(row[c["label_seq_id"]]), "aa": AA3[comp],
                "xyz": [float(row[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                "bfactor": float(row[c["B_iso_or_equiv"]]), "alt": alt,
                "ins": row[c["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in c else "?",
                "auth_asym": row[c["auth_asym_id"]] if "auth_asym_id" in c else "?",
                "comp": comp,
            })
        except (ValueError, KeyError):
            continue
    return atoms, sequence, variants


def parse_af(text: str):
    out = {}
    for tags, row in cif_loops(text):
        if not any(t.startswith("_atom_site.") for t in tags):
            continue
        names = [t.split(".", 1)[1] for t in tags]
        c = {name: i for i, name in enumerate(names)}
        needed = {"label_atom_id", "label_comp_id", "label_alt_id", "auth_seq_id",
                  "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv", "pdbx_PDB_model_num"}
        if not needed <= set(c):
            continue
        if row[c["label_atom_id"]] != "CA" or row[c["pdbx_PDB_model_num"]] != "1":
            continue
        comp, alt = row[c["label_comp_id"]], row[c["label_alt_id"]]
        if comp not in AA3 or alt not in (".", "?", "", "A"):
            continue
        try:
            pos = int(row[c["auth_seq_id"]])
            out[pos] = {"aa": AA3[comp], "xyz": [float(row[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                        "bfactor": float(row[c["B_iso_or_equiv"]])}
        except (ValueError, KeyError):
            continue
    return out


def align_identical(ref_seq, af_seq):
    n, m = len(ref_seq), len(af_seq)
    if n == 0 or m == 0:
        return {}
    score = np.zeros((n + 1, m + 1), dtype=np.float64)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    score[:, 0] = -2 * np.arange(n + 1)
    score[0, :] = -2 * np.arange(m + 1)
    trace[1:, 0], trace[0, 1:] = 1, 2
    aa = np.asarray([c for _, c in af_seq])
    for i in range(1, n + 1):
        diag = score[i - 1, :-1] + np.where(aa == ref_seq[i - 1][1], 2, -1)
        up = score[i - 1, 1:] - 2
        best = np.maximum(diag, up)
        choice = np.where(diag >= up, 0, 1)
        for j in range(1, m + 1):
            left = score[i, j - 1] - 2
            if best[j - 1] >= left:
                score[i, j], trace[i, j] = best[j - 1], choice[j - 1]
            else:
                score[i, j], trace[i, j] = left, 2
    pairs, i, j = {}, n, m
    while i and j:
        step = trace[i, j]
        if step == 0:
            if ref_seq[i - 1][1] == af_seq[j - 1][1]:
                pairs[ref_seq[i - 1][0]] = af_seq[j - 1][0]
            i, j = i - 1, j - 1
        elif step == 1:
            i -= 1
        else:
            j -= 1
    return pairs


def recompute(mapping, ref_bytes, af_bytes):
    ref_text = ref_bytes.decode("utf-8", errors="replace")
    af_text = af_bytes.decode("utf-8", errors="replace")
    ref_atoms, ref_seq, variants = parse_ref_full(ref_text, mapping["struct_asym_id"], mapping["entity_id"])
    af_atoms = parse_af(af_text)
    ls, le = int(mapping["start"]["residue_number"]), int(mapping["end"]["residue_number"])
    us, ue = int(mapping["unp_start"]), int(mapping["unp_end"])
    if le - ls == ue - us:
        resmap = {ls + k: us + k for k in range(le - ls + 1)}
        map_method = "sifts_offset"
    else:
        ref_sub = [(p, ref_seq[p]) for p in range(ls, le + 1) if p in ref_seq]
        af_sub = [(p, af_atoms[p]["aa"]) for p in range(us, ue + 1) if p in af_atoms]
        resmap = align_identical(ref_sub, af_sub)
        map_method = "sequence_alignment"
    pairs, seen, n_pairs, n_mismatch = [], set(), 0, 0
    ref_index = []
    for i, atom in enumerate(ref_atoms):
        pos = resmap.get(atom["label_seq"])
        if pos is None or pos in seen or pos not in af_atoms:
            continue
        n_pairs += 1
        if af_atoms[pos]["aa"] != atom["aa"]:
            n_mismatch += 1
            continue
        seen.add(pos)
        pairs.append((pos, atom, af_atoms[pos]))
        ref_index.append(i)
    if n_pairs and n_mismatch / n_pairs > 0.10:
        return {"status": "E427_SEQUENCE_IDENTITY_FAILED", "rows": [], "n_pairs": n_pairs,
                "n_identity_mismatch": n_mismatch, "map_method": map_method,
                "ref_atoms": ref_atoms, "resmap": resmap, "variants": variants}
    if len(pairs) < 30:
        return {"status": "E422_INSUFFICIENT_ALIGNED", "rows": [], "n_pairs": n_pairs,
                "n_identity_mismatch": n_mismatch, "map_method": map_method,
                "ref_atoms": ref_atoms, "resmap": resmap, "variants": variants}
    mobile = np.asarray([af["xyz"] for _p, _r, af in pairs], dtype=np.float64)
    reference = np.asarray([r["xyz"] for _p, r, _af in pairs], dtype=np.float64)
    mc = np.asarray([np.sum(mobile[:, j], dtype=np.float64) / len(mobile) for j in range(3)])
    rc = np.asarray([np.sum(reference[:, j], dtype=np.float64) / len(reference) for j in range(3)])
    h = (mobile - mc).T @ (reference - rc)
    u, _s, vt = np.linalg.svd(h, full_matrices=False)
    d = np.diag([1.0, 1.0, 1.0 if np.linalg.det(u @ vt) > 0 else -1.0])
    sup = (mobile - mc) @ (u @ d @ vt) + rc
    dist = np.linalg.norm(sup - reference, axis=1)
    diff_ref = reference[:, None, :] - reference[None, :, :]
    dist_ref_sq = np.sum(diff_ref * diff_ref, axis=-1)
    np.fill_diagonal(dist_ref_sq, np.inf)
    neighbor = dist_ref_sq < 225.0
    diff_sup = sup[:, None, :] - sup[None, :, :]
    delta = np.abs(np.sqrt(np.sum(diff_sup * diff_sup, axis=-1)) - np.sqrt(dist_ref_sq))
    preserved = np.zeros((*delta.shape, 4), dtype=bool)
    for j, threshold in enumerate((0.5, 1.0, 2.0, 4.0)):
        preserved[:, :, j] = delta < threshold
    counts = neighbor.sum(axis=1)
    lddt = np.full(len(pairs), np.nan)
    good = counts > 0
    lddt[good] = (preserved & neighbor[:, :, None]).sum(axis=(1, 2))[good] / (4.0 * counts[good])
    bf = np.asarray([a["bfactor"] for a in ref_atoms], dtype=np.float64)
    med = np.median(bf)
    q25, q75 = np.percentile(bf, [25, 75], method="linear")
    iqr = q75 - q25 or 1.0
    bz = (bf - med) / iqr
    rows = []
    for j, (pos, atom, af) in enumerate(pairs):
        rows.append({
            "accession": "", "entry": "", "chain": "", "uniprot_pos": int(pos),
            "plddt": round(float(af["bfactor"]), 4),
            "lddt": None if np.isnan(lddt[j]) else round(float(lddt[j]), 4),
            "d_kabsch": round(float(dist[j]), 4),
            "b_factor_z": round(float(bz[ref_index[j]]), 4),
        })
    return {"status": "ok", "rows": rows, "n_pairs": n_pairs,
            "n_identity_mismatch": n_mismatch, "map_method": map_method,
            "ref_atoms": ref_atoms, "resmap": resmap, "variants": variants}


def load_receipt_bytes(checkpoint, url):
    matches = [r for r in checkpoint.get("receipts", []) if r.get("url") == url and r.get("raw_name")]
    for r in matches:
        path = RAW / r["raw_name"]
        if not path.is_file():
            continue
        data = path.read_bytes()
        if sha(data) != r["raw_name"] or r.get("sha256") not in (None, r["raw_name"]):
            raise AssertionError(f"receipt/raw hash mismatch for {url}")
        return data, r
    return None, None


def group_rows(ck):
    out = defaultdict(list)
    for row in ck.get("rows", []):
        out[(row["entry"].lower(), row["chain"], row["accession"])].append(row)
    for rows in out.values():
        rows.sort(key=lambda r: (r["uniprot_pos"], r["lddt"] if r["lddt"] is not None else -1))
    return dict(out)


def checkpoint_integrity(rerun_dir: Path | None):
    details, all_good = {}, True
    for batch in range(1, 5):
        op = ORIG / f"e422_batch{batch:02d}_checkpoint.json"
        cp = CORR / f"e422_batch{batch:02d}_checkpoint.json"
        ob, cb = op.read_bytes(), cp.read_bytes()
        old, new = json.loads(ob), json.loads(cb)
        digest_ok = canon_digest(new) == new.get("sha256")
        receipt_equal = old.get("receipts") == new.get("receipts")
        row = {
            "batch": batch, "corrected_sha256": new.get("sha256"),
            "canonical_digest_valid": digest_ok,
            "receipts_equal_original": receipt_equal,
            "manifest_sha_equal_original": new.get("manifest_checkpoint_sha256") == old.get("manifest_checkpoint_sha256"),
            "corrected_rows": len(new.get("rows", [])),
            "original_rows": len(old.get("rows", [])),
            "batch1_byte_identical": (ob == cb) if batch == 1 else None,
            "provenance_present": batch == 1 or isinstance(new.get("e427_relabel"), dict),
            "original_sha_record_matches": batch == 1 or new.get("e427_relabel", {}).get("original_checkpoint_sha256") == old.get("sha256"),
        }
        if batch == 1:
            row["canonical_digest_valid"] = digest_ok
        if rerun_dir is not None:
            rp = rerun_dir / cp.name
            if rp.exists():
                rerun = json.loads(rp.read_text())
                normalized_current, normalized_rerun = copy.deepcopy(new), copy.deepcopy(rerun)
                for doc in (normalized_current, normalized_rerun):
                    doc.pop("sha256", None)
                    prov = doc.get("e427_relabel")
                    if isinstance(prov, dict):
                        prov.pop("replay_utc", None)
                row["rerun_present"] = True
                row["rerun_byte_identical"] = cb == rp.read_bytes()
                row["rerun_equal_after_normalizing_replay_time_and_digest"] = normalized_current == normalized_rerun
                row["rerun_sha256"] = rerun.get("sha256")
            else:
                row["rerun_present"] = False
                row["rerun_byte_identical"] = False
                row["rerun_equal_after_normalizing_replay_time_and_digest"] = False
        all_good &= digest_ok and receipt_equal and row["manifest_sha_equal_original"] and row["provenance_present"] and row["original_sha_record_matches"]
        if batch == 1:
            all_good &= row["batch1_byte_identical"]
        details[str(batch)] = row
    report = load_json(ROOT / "results/e427/replay_report.json")
    report_checks = all(all(v for v in b["self_check"].values()) and
                        b["corrected_manifest_sha_equal"] and
                        b["receipts_unconsumed_corrected"] == 0
                        for b in report["batches"].values())
    details["replay_report"] = {"verdict": report.get("verdict"), "all_self_checks_true": report_checks,
                                "batches": report["batches"]}
    all_good &= report.get("verdict") == "OK" and report_checks
    return details, all_good


def segment_list(checkpoint, accession, entry, chain):
    url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry.lower()}"
    raw, receipt = load_receipt_bytes(checkpoint, url)
    if raw is None:
        return [], None
    doc = json.loads(raw)
    mappings = doc[entry.lower()]["UniProt"].get(accession, {}).get("mappings", [])
    return [m for m in mappings if str(m.get("chain_id")) == str(chain)], receipt


def sample_segments():
    original, corrected, rows_orig, rows_corr = {}, {}, {}, {}
    for batch in range(1, 5):
        original[batch] = load_json(ORIG / f"e422_batch{batch:02d}_checkpoint.json")
        corrected[batch] = load_json(CORR / f"e422_batch{batch:02d}_checkpoint.json")
        rows_orig[batch] = group_rows(original[batch])
        rows_corr[batch] = group_rows(corrected[batch])
    targets = [
        ("9qj6", "A", "K7PQ54"), ("9qj6", "B", "K7PQ54"),
        ("9sq1", "A", "Q5JF22"), ("9sq1", "B", "Q5JF22"),
        ("9sq1", "C", "Q5JF22"), ("9sq1", "D", "Q5JF22"),
        ("9y8y", "A", "P07342"), ("9y8y", "B", "P07342"),
        ("9wvj", "A", "A0A4Q0WMG2"), ("9wvj", "B", "A0A4Q0WMG2"),
        ("31rm", "A", "Q9HWK6"), ("31rm", "B", "Q9HWK6"),
        ("31rm", "C", "Q9HWK6"), ("31rm", "D", "Q9HWK6"),
        ("31al", "A", "D1MPT3"), ("31al", "B", "D1MPT3"),
        ("31jt", "A", "D1MPT3"), ("31jt", "B", "D1MPT3"),
        ("9qj7", "A", "K7PQ54"), ("9qj7", "B", "K7PQ54"),
        ("9qj8", "A", "K7PQ54"), ("9qj8", "B", "K7PQ54"),
        ("9qjd", "A", "K7PQ54"), ("9qjd", "B", "K7PQ54"),
        ("9qjf", "A", "K7PQ54"), ("9qjf", "B", "K7PQ54"),
    ]
    groups_by_batch = {}
    for b in range(1, 5):
        groups_by_batch[b] = rows_corr[b]
    chosen = []
    used = set()

    def add(entry, chain, accession, reason):
        key = (entry, chain, accession)
        batch_matches = [b for b in range(1, 5) if key in groups_by_batch[b]]
        if len(batch_matches) != 1:
            raise AssertionError(f"expected one corrected checkpoint group for {key}, got {batch_matches}")
        batch = batch_matches[0]
        mappings, map_receipt = segment_list(corrected[batch], accession, entry, chain)
        if len(mappings) != 1:
            raise AssertionError(f"expected one SIFTS segment for {key}, got {len(mappings)}")
        ref_url = f"https://www.ebi.ac.uk/pdbe/entry-files/{entry.lower()}.cif"
        af_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif"
        ref_bytes, ref_receipt = load_receipt_bytes(corrected[batch], ref_url)
        af_bytes, af_receipt = load_receipt_bytes(corrected[batch], af_url)
        if ref_bytes is None or af_bytes is None:
            raise AssertionError(f"missing receipted reference/AF bytes for {key}")
        if key in used:
            raise AssertionError(f"duplicate sample group {key}")
        chosen.append({"batch": batch, "entry": entry, "chain": chain, "accession": accession,
                       "mapping": mappings[0], "segment_index": 0, "map_receipt": map_receipt,
                       "ref_bytes": ref_bytes, "ref_receipt": ref_receipt,
                       "af_bytes": af_bytes, "af_receipt": af_receipt,
                       "reason": reason, "expected_rows": groups_by_batch[batch][key]})
        used.add(key)

    for entry, chain, accession in targets:
        add(entry, chain, accession, "required accession / diagnostic target")

    candidates = []
    for batch in range(1, 5):
        for key, new_rows in rows_corr[batch].items():
            old_rows = rows_orig[batch].get(key, [])
            if old_rows != new_rows or key in used:
                continue
            entry, chain, accession = key
            mappings, map_receipt = segment_list(corrected[batch], accession, entry, chain)
            if len(mappings) != 1:
                continue
            ref_url = f"https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"
            af_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif"
            ref_bytes, ref_receipt = load_receipt_bytes(corrected[batch], ref_url)
            af_bytes, af_receipt = load_receipt_bytes(corrected[batch], af_url)
            if ref_bytes is None or af_bytes is None:
                continue
            candidates.append((batch, key, mappings[0], map_receipt,
                               ref_bytes, ref_receipt, af_bytes, af_receipt, new_rows))
    candidates.sort(key=lambda x: (x[1][0], x[1][1], x[1][2]))
    if len(candidates) < 10:
        raise AssertionError(f"only {len(candidates)} unchanged groups available for random sample")
    rng = random.Random(427)
    for item in rng.sample(candidates, 10):
        batch, (entry, chain, accession), mapping, map_receipt, ref_bytes, ref_receipt, af_bytes, af_receipt, expected = item
        chosen.append({"batch": batch, "entry": entry, "chain": chain, "accession": accession,
                       "mapping": mapping, "segment_index": 0, "map_receipt": map_receipt,
                       "ref_bytes": ref_bytes, "ref_receipt": ref_receipt,
                       "af_bytes": af_bytes, "af_receipt": af_receipt,
                       "reason": "seed-427 sample from unchanged groups", "expected_rows": expected})

    results, all_match = [], True
    for item in chosen:
        mapping = item["mapping"]
        calc = recompute(mapping, item["ref_bytes"], item["af_bytes"])
        got = copy.deepcopy(calc["rows"])
        for row in got:
            row.update({"entry": item["entry"], "chain": item["chain"], "accession": item["accession"]})
        got.sort(key=lambda r: r["uniprot_pos"])
        expected = item["expected_rows"]
        exp_by_pos = {r["uniprot_pos"]: r for r in expected}
        got_by_pos = {r["uniprot_pos"]: r for r in got}
        shared = sorted(set(exp_by_pos) & set(got_by_pos))
        mismatches = []
        max_abs = 0.0
        for pos in shared:
            for field in ("plddt", "lddt", "d_kabsch", "b_factor_z"):
                a, b = exp_by_pos[pos].get(field), got_by_pos[pos].get(field)
                if a is None or b is None:
                    if a != b:
                        mismatches.append({"pos": pos, "field": field, "expected": a, "independent": b})
                else:
                    delta = abs(float(a) - float(b))
                    max_abs = max(max_abs, delta)
                    if delta > 0.0001:
                        mismatches.append({"pos": pos, "field": field, "expected": a, "independent": b, "abs_diff": delta})
        row_match = (calc["status"] == "ok" and len(expected) == len(got) and
                     set(exp_by_pos) == set(got_by_pos) and not mismatches)
        all_match &= row_match
        ref_atoms = calc.get("ref_atoms", [])
        dup_seq = sum(1 for n in Counter(a["label_seq"] for a in ref_atoms).values() if n > 1)
        seg = mapping
        results.append({
            "batch": item["batch"], "entry": item["entry"], "chain": item["chain"],
            "accession": item["accession"], "segment_index": item["segment_index"], "reason": item["reason"],
            "struct_asym_id": mapping.get("struct_asym_id"), "entity_id": mapping.get("entity_id"),
            "label_span": [seg["start"].get("residue_number"), seg["end"].get("residue_number")],
            "uniprot_span": [seg.get("unp_start"), seg.get("unp_end")],
            "mapping_method": calc.get("map_method"), "status": calc["status"],
            "n_identity_pairs": calc.get("n_pairs"), "n_identity_mismatch": calc.get("n_identity_mismatch"),
            "corrected_checkpoint_rows": len(expected), "independent_rows": len(got),
            "row_values_match": row_match, "max_abs_numeric_difference": max_abs,
            "mismatches_first_10": mismatches[:10],
            "ref_atom_sites": len(ref_atoms), "duplicate_label_seq_positions": dup_seq,
            "accepted_altloc_counts": dict(Counter(a["alt"] for a in ref_atoms)),
            "mse_ca_count": sum(a["comp"] == "MSE" for a in ref_atoms),
            "nontrivial_insertion_code_ca_count": sum(a["ins"] not in (".", "?", "", "null") for a in ref_atoms),
            "entity_poly_seq_duplicate_positions": sum(len(v) > 1 for v in calc.get("variants", {}).values()),
            "entity_poly_seq_variant_positions": sum(len({aa for _mon, aa in v}) > 1 for v in calc.get("variants", {}).values()),
            "auth_asym_ids": sorted({a["auth_asym"] for a in ref_atoms}),
            "raw_hashes": {"mapping": item["map_receipt"]["raw_name"],
                           "reference": item["ref_receipt"]["raw_name"],
                           "af": item["af_receipt"]["raw_name"]},
        })
    return {"n_segments": len(results), "targeted_segments": sum(x["reason"] != "seed-427 sample from unchanged groups" for x in results),
            "random_unchanged_segments": sum(x["reason"] == "seed-427 sample from unchanged groups" for x in results),
            "all_rows_match": all_match, "segments": results}


def nine_zxa_check():
    for batch in range(1, 5):
        ck = load_json(CORR / f"e422_batch{batch:02d}_checkpoint.json")
        entry = "9zxa"
        mapping_url = f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry}"
        map_raw, map_receipt = load_receipt_bytes(ck, mapping_url)
        if map_raw is None:
            continue
        mapping_doc = json.loads(map_raw)
        ref_url = f"https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"
        successful_ref = [r for r in ck.get("receipts", []) if r.get("url") == ref_url and r.get("raw_name")]
        failures = [f for f in ck.get("exclusions", {}).get("triple_failures", []) if f.get("entry", "").lower() == entry]
        q92802 = [f for f in failures if f.get("accession") == "Q92802"]
        return {
            "batch": batch, "mapping_receipt_status": map_receipt.get("status"),
            "mapping_response_has_entry": entry in mapping_doc,
            "q92802_mapping_count": len(mapping_doc.get(entry, {}).get("UniProt", {}).get("Q92802", {}).get("mappings", [])),
            "reference_cif_successful_receipt_count": len(successful_ref),
            "q92802_failure_codes": sorted(Counter(f.get("code") for f in q92802).items()),
            "q92802_failure_records": q92802,
        }
    return {"found": False}


def mapping_structure_summary():
    groups = defaultdict(list)
    receipts_checked = set()
    for batch in range(1, 5):
        ck = load_json(ORIG / f"e422_batch{batch:02d}_checkpoint.json")
        for receipt in ck.get("receipts", []):
            url = receipt.get("url", "")
            if "/mappings/uniprot/" not in url or not receipt.get("raw_name"):
                continue
            raw, _used = load_receipt_bytes(ck, url)
            if raw is None:
                continue
            receipts_checked.add(receipt["raw_name"])
            doc = json.loads(raw)
            for entry, block in doc.items():
                if not isinstance(block, dict) or not isinstance(block.get("UniProt"), dict):
                    continue
                for accession, acc_block in block["UniProt"].items():
                    if not isinstance(acc_block, dict):
                        continue
                    for mapping in acc_block.get("mappings", []):
                        try:
                            key = (entry.lower(), accession, str(mapping["chain_id"]))
                            start = int(mapping["start"]["residue_number"])
                            end = int(mapping["end"]["residue_number"])
                            us, ue = int(mapping["unp_start"]), int(mapping["unp_end"])
                        except (KeyError, TypeError, ValueError):
                            continue
                        groups[key].append({"struct_asym_id": mapping.get("struct_asym_id"),
                                            "label": [start, end], "uniprot": [us, ue]})
    multi = []
    overlaps = []
    multi_asym = []
    for key, segments in sorted(groups.items()):
        if len(segments) > 1:
            item = {"entry": key[0], "accession": key[1], "chain": key[2],
                    "n_segments": len(segments), "struct_asym_ids": sorted({str(s["struct_asym_id"]) for s in segments}),
                    "segments": segments}
            multi.append(item)
            if len({str(s["struct_asym_id"]) for s in segments}) > 1:
                multi_asym.append(item)
            for i, a in enumerate(segments):
                for b in segments[i + 1:]:
                    label_overlap = max(a["label"][0], b["label"][0]) <= min(a["label"][1], b["label"][1])
                    uniprot_overlap = max(a["uniprot"][0], b["uniprot"][0]) <= min(a["uniprot"][1], b["uniprot"][1])
                    if label_overlap or uniprot_overlap:
                        overlaps.append({"entry": key[0], "accession": key[1], "chain": key[2],
                                         "label_overlap": label_overlap, "uniprot_overlap": uniprot_overlap,
                                         "segments": [a, b]})
    return {"verified_mapping_response_hashes": len(receipts_checked),
            "unique_accession_entry_chain_groups": len(groups),
            "groups_with_multiple_segments": len(multi),
            "groups_with_multiple_struct_asym_ids": len(multi_asym),
            "overlapping_segment_pairs": len(overlaps),
            "multi_segment_examples_first_20": multi[:20],
            "multi_asym_examples_first_20": multi_asym[:20],
            "overlap_examples_first_20": overlaps[:20]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun-dir", type=Path, default=None)
    args = ap.parse_args()
    details, checkpoints_ok = checkpoint_integrity(args.rerun_dir)
    sample = sample_segments()
    audit = {
        "method": "Independent standard-library mmCIF/SIFTS parser plus NumPy Kabsch/lDDT; no moluq/e427 module imports",
        "checkpoint_integrity": details,
        "checkpoint_integrity_pass": checkpoints_ok,
        "independent_segment_sample": sample,
        "mapping_structure_summary": mapping_structure_summary(),
        "nine_zxa": nine_zxa_check(),
        "verdict": "PASS" if checkpoints_ok and sample["all_rows_match"] else "FAIL",
    }
    out = HERE / "independent_label_sample.json"
    out.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verdict": audit["verdict"], "checkpoint_integrity_pass": checkpoints_ok,
                      "sample_segments": sample["n_segments"], "targeted": sample["targeted_segments"],
                      "random_unchanged": sample["random_unchanged_segments"],
                      "all_rows_match": sample["all_rows_match"], "nine_zxa": audit["nine_zxa"],
                      "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
