#!/usr/bin/env python3
"""Independent checks of e427's e421 relabel outputs; imports no repo code."""
from __future__ import annotations

import collections
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RAW = ROOT / "data/e427/e421_raw"
NEW_PATH = ROOT / "results/e427/e421_relabel/census_results.json"
OLD_PATH = ROOT / "results/e421/census_results.json"
RECEIPTS_PATH = ROOT / "results/e427/e421_relabel/receipts.json"
SUMMARY_PATH = ROOT / "results/e427/e421_relabel/summary.json"
QUAL_PATH = ROOT / "results/e421/qualifying_entries.json"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def load(path):
    return json.loads(path.read_text())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def tokenizer(text):
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
                raise ValueError("unterminated mmCIF semicolon block")
            yield text[i + 1:end]
            i = end + 2
        elif c in "'\"":
            quote, end = c, i + 1
            while True:
                end = text.find(quote, end)
                if end < 0:
                    raise ValueError("unterminated mmCIF quote")
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


def loops(text):
    tok = list(tokenizer(text))
    i = 0
    while i < len(tok):
        if tok[i] != "loop_":
            i += 1
            continue
        i += 1
        tags = []
        while i < len(tok) and tok[i].startswith("_"):
            tags.append(tok[i])
            i += 1
        if not tags:
            continue
        width = len(tags)
        while i < len(tok):
            if tok[i] in ("loop_", "stop_") or tok[i].startswith(("_", "data_")):
                break
            row = tok[i:i + width]
            if len(row) != width:
                raise ValueError("truncated mmCIF loop")
            yield tags, row
            i += width


def parse_ref(text, asym, entity_id):
    atoms, sequence = [], {}
    duplicate_keys = collections.Counter()
    for tags, row in loops(text):
        names = [t.split(".", 1)[1] for t in tags]
        col = {name: i for i, name in enumerate(names)}
        if any(t.startswith("_entity_poly_seq.") for t in tags):
            if not {"entity_id", "num", "mon_id"} <= col.keys():
                continue
            try:
                if row[col["entity_id"]] == str(entity_id):
                    sequence[int(row[col["num"]])] = AA3.get(row[col["mon_id"]], "X")
            except (ValueError, KeyError):
                pass
            continue
        if not any(t.startswith("_atom_site.") for t in tags):
            continue
        need = {"label_atom_id", "label_comp_id", "label_alt_id", "label_seq_id",
                "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z",
                "B_iso_or_equiv", "pdbx_PDB_model_num"}
        if not need <= col.keys():
            continue
        if (row[col["label_atom_id"]] != "CA" or row[col["label_asym_id"]] != asym
                or row[col["pdbx_PDB_model_num"]] != "1"
                or row[col["label_comp_id"]] not in AA3
                or row[col["label_alt_id"]] not in (".", "?", "", "A")):
            continue
        try:
            atom = {
                "label_seq": int(row[col["label_seq_id"]]),
                "aa": AA3[row[col["label_comp_id"]]],
                "xyz": [float(row[col[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                "b": float(row[col["B_iso_or_equiv"]]),
                "alt": row[col["label_alt_id"]],
            }
        except (ValueError, KeyError):
            continue
        atoms.append(atom)
        duplicate_keys[atom["label_seq"]] += 1
    return atoms, sequence, {k: v for k, v in duplicate_keys.items() if v > 1}


def parse_af(text):
    af = {}
    for tags, row in loops(text):
        if not any(t.startswith("_atom_site.") for t in tags):
            continue
        names = [t.split(".", 1)[1] for t in tags]
        col = {name: i for i, name in enumerate(names)}
        need = {"label_atom_id", "label_comp_id", "auth_seq_id", "Cartn_x", "Cartn_y",
                "Cartn_z", "B_iso_or_equiv", "pdbx_PDB_model_num", "label_alt_id"}
        if not need <= col.keys():
            continue
        if (row[col["label_atom_id"]] != "CA" or row[col["pdbx_PDB_model_num"]] != "1"
                or row[col["label_comp_id"]] not in AA3
                or row[col["label_alt_id"]] not in (".", "?", "", "A")):
            continue
        try:
            af[int(row[col["auth_seq_id"]])] = {
                "aa": AA3[row[col["label_comp_id"]]],
                "xyz": [float(row[col[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")],
                "b": float(row[col["B_iso_or_equiv"]]),
            }
        except (ValueError, KeyError):
            continue
    return af


def align_map(ref_seq, af_seq):
    n, m = len(ref_seq), len(af_seq)
    if not n or not m:
        return {}
    score = np.zeros((n + 1, m + 1), dtype=np.int64)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    score[:, 0] = -2 * np.arange(n + 1)
    score[0, :] = -2 * np.arange(m + 1)
    trace[1:, 0], trace[0, 1:] = 1, 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + (2 if ref_seq[i - 1][1] == af_seq[j - 1][1] else -1)
            up, left = score[i - 1, j] - 2, score[i, j - 1] - 2
            best = max(diag, up, left)
            score[i, j] = best
            trace[i, j] = 0 if diag == best else (1 if up == best else 2)
    out, i, j = {}, n, m
    while i and j:
        step = trace[i, j]
        if step == 0:
            if ref_seq[i - 1][1] == af_seq[j - 1][1]:
                out[ref_seq[i - 1][0]] = af_seq[j - 1][0]
            i, j = i - 1, j - 1
        elif step == 1:
            i -= 1
        else:
            j -= 1
    return out


def entry_receipt_index(receipts):
    by_url = collections.defaultdict(list)
    for receipt in receipts:
        by_url[receipt["url"]].append(receipt)
    return by_url


def raw_for_url(url, by_url):
    found = [r for r in by_url.get(url, []) if r.get("status") == 200]
    if not found:
        return None
    # Use the first successful response, which is the request stream's first success.
    r = found[0]
    body = (RAW / r["sha256"]).read_bytes()
    assert sha(body) == r["sha256"] and len(body) == r["bytes_len"]
    return body


def independent_entry(pdb, accession, cif_url, by_url, compare_labels=False):
    base = f"https://www.ebi.ac.uk/pdbe/"
    mb = raw_for_url(f"{base}api/v2/mappings/uniprot/{pdb}", by_url)
    if mb is None:
        return {"status": "mapping_fail"}
    obj = json.loads(mb)
    mappings = obj[pdb]["UniProt"].get(accession, {}).get("mappings", [])
    segments = [m for m in mappings if m.get("chain_id") == "A"]
    if not segments:
        return {"status": "mapping_acc_absent", "n_segments": 0}
    af_bytes = raw_for_url(cif_url, by_url)
    if af_bytes is None:
        return {"status": "af_fetch_fail", "n_segments": len(segments)}
    ref_bytes = raw_for_url(f"{base}entry-files/{pdb}.cif", by_url)
    if ref_bytes is None:
        return {"status": "ref_fetch_fail", "n_segments": len(segments)}
    af = parse_af(af_bytes.decode("utf-8", "replace"))
    ref_text = ref_bytes.decode("utf-8", "replace")
    pairs, seen, n_pair, n_mismatch = [], set(), 0, 0
    segment_info = []
    dup_groups = 0
    for seg in segments:
        ls, le = int(seg["start"]["residue_number"]), int(seg["end"]["residue_number"])
        us, ue = int(seg["unp_start"]), int(seg["unp_end"])
        asym, entity = seg["struct_asym_id"], str(seg["entity_id"])
        ref_atoms, ref_seq, dups = parse_ref(ref_text, asym, entity)
        dup_groups += len(dups)
        if le - ls == ue - us:
            resmap = {ls + k: us + k for k in range(le - ls + 1)}
            method = "offset"
        else:
            rs = [(p, ref_seq[p]) for p in range(ls, le + 1) if p in ref_seq]
            usq = [(p, af[p]["aa"]) for p in range(us, ue + 1) if p in af]
            resmap = align_map(rs, usq)
            method = "alignment"
        before_pairs, before_mismatch = n_pair, n_mismatch
        for atom in ref_atoms:
            unp = resmap.get(atom["label_seq"])
            if unp is None or unp in seen or unp not in af:
                continue
            n_pair += 1
            if af[unp]["aa"] != atom["aa"]:
                n_mismatch += 1
                continue
            seen.add(unp)
            pairs.append((atom, af[unp], unp))
        segment_info.append({"asym": asym, "entity": entity, "map_method": method,
                             "n_pairs": n_pair - before_pairs,
                             "n_mismatch": n_mismatch - before_mismatch,
                             "duplicate_label_seq_groups": len(dups)})
    if n_pair and n_mismatch / n_pair > 0.10:
        return {"status": "sequence_identity_failed", "n_pairs": n_pair,
                "n_identity_mismatch": n_mismatch, "segments": segment_info,
                "duplicate_label_seq_groups": dup_groups}
    if len(pairs) < 30:
        return {"status": "insufficient_aligned", "n_aligned": len(pairs),
                "n_pairs": n_pair, "n_identity_mismatch": n_mismatch,
                "segments": segment_info, "duplicate_label_seq_groups": dup_groups}
    answer = {"status": "ok", "n_aligned": len(pairs), "n_pairs": n_pair,
              "n_identity_mismatch": n_mismatch, "segments": segment_info,
              "duplicate_label_seq_groups": dup_groups}
    if not compare_labels:
        return answer
    mobile = np.asarray([af_atom["xyz"] for _ref, af_atom, _pos in pairs], dtype=np.float64)
    reference = np.asarray([atom["xyz"] for atom, _af, _pos in pairs], dtype=np.float64)
    mc = np.asarray([math.fsum(mobile[:, j]) for j in range(3)]) / len(pairs)
    rc = np.asarray([math.fsum(reference[:, j]) for j in range(3)]) / len(pairs)
    u, _s, vt = np.linalg.svd((mobile - mc).T @ (reference - rc), full_matrices=False)
    rot = np.diag([1.0, 1.0, 1.0 if np.linalg.det(u @ vt) > 0 else -1.0])
    fitted = (mobile - mc) @ (u @ rot @ vt) + rc
    dists = np.linalg.norm(fitted - reference, axis=1)
    ref_diff = reference[:, None, :] - reference[None, :, :]
    ref_sq = np.sum(ref_diff * ref_diff, axis=-1)
    np.fill_diagonal(ref_sq, np.inf)
    neighbors = ref_sq < 225.0
    fitted_diff = fitted[:, None, :] - fitted[None, :, :]
    delta = np.abs(np.sqrt(np.sum(fitted_diff * fitted_diff, axis=-1)) - np.sqrt(ref_sq))
    preserved = np.zeros((*delta.shape, 4), dtype=bool)
    for k, threshold in enumerate((0.5, 1.0, 2.0, 4.0)):
        preserved[:, :, k] = delta < threshold
    n_neighbors = neighbors.sum(axis=1)
    lddt = np.full(len(pairs), np.nan)
    valid = n_neighbors > 0
    lddt[valid] = (preserved & neighbors[:, :, None]).sum(axis=(1, 2))[valid] / (4 * n_neighbors[valid])
    bfs = np.asarray([atom["b"] for atom, _af, _pos in pairs], dtype=np.float64)
    med = np.median(bfs)
    q25, q75 = np.percentile(bfs, [25, 75], method="linear")
    iqr = q75 - q25 or 1.0
    z = (bfs - med) / iqr
    rows = []
    for i, (atom, af_atom, pos) in enumerate(pairs):
        lv = float(lddt[i]) if not np.isnan(lddt[i]) else -1.0
        d = float(dists[i])
        e = "correct" if lv >= 0.60 and d <= 4.0 else ("error" if lv < 0.60 and d > 4.0 else "U")
        f = "flexible" if z[i] >= 1.0 else ("ordered" if z[i] <= -1.0 else "U")
        cell = "U"
        if e != "U" and f != "U":
            cell = {("correct", "ordered"): "OC", ("correct", "flexible"): "A",
                    ("error", "ordered"): "B", ("error", "flexible"): "O"}[(e, f)]
        rows.append({"uniprot_pos": int(pos), "e": e, "f": f, "cell": cell,
                     "lddt": round(float(lddt[i]), 4), "dist": round(d, 4),
                     "b": atom["b"], "z": round(float(z[i]), 4),
                     "plddt": af_atom["b"]})
    answer["labels"] = rows
    return answer


def average_ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = ((i + 1) + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = rank
        i = j
    return ranks


def spearman(xs, ys):
    if len(xs) < 2:
        return None
    x, y = average_ranks(xs), average_ranks(ys)
    mx, my = math.fsum(x) / len(x), math.fsum(y) / len(y)
    xx = math.fsum((v - mx) ** 2 for v in x)
    yy = math.fsum((v - my) ** 2 for v in y)
    xy = math.fsum((x[i] - mx) * (y[i] - my) for i in range(len(x)))
    return xy / math.sqrt(xx * yy)


def summarize(results, labs=None):
    ok = [r for r in results.values() if r.get("status") == "ok"]
    if labs is None:
        labs = [label for r in ok for label in r["labels"]]
    cells = collections.Counter(row["cell"] for row in labs)
    cells_unique = {cell: len({(r["accession"], r["uniprot_pos"])
                               for r in labs if r["cell"] == cell}) for cell in cells}
    bins, gradient = [(0, 50), (50, 70), (70, 80), (80, 90), (90, 100.01)], []
    for lo, hi in bins:
        selected = [r for r in labs if r.get("lddt") is not None
                    and not (isinstance(r["lddt"], str) and r["lddt"] == "nan")
                    and not math.isnan(float(r["lddt"])) and lo <= r["plddt"] < hi]
        gradient.append({"bin": [lo, hi], "n": len(selected),
                         "err_rate_lddt_lt_0.60": (sum(float(r["lddt"]) < 0.60 for r in selected) / len(selected))
                         if selected else None})
    finite = [r for r in labs if r.get("lddt") is not None
              and not (isinstance(r["lddt"], str) and r["lddt"] == "nan")
              and not math.isnan(float(r["lddt"]))]
    return {"entries": len(results),
            "status_counts": dict(collections.Counter(r["status"] for r in results.values())),
            "ok_entries": len(ok), "ok_accessions": len({r["accession"] for r in ok}),
            "residue_rows": len(labs), "cells_rows": dict(cells), "cells_unique": cells_unique,
            "gradient": gradient,
            "spearman_plddt_lddt": spearman([r["plddt"] for r in finite], [float(r["lddt"]) for r in finite]),
            "n_spearman": len(finite)}


def main():
    new, old = load(NEW_PATH), load(OLD_PATH)
    receipts, declared = load(RECEIPTS_PATH), load(SUMMARY_PATH)
    raw_paths = list(RAW.iterdir())
    bad_files, raw_shas = [], set()
    for path in raw_paths:
        if not path.is_file():
            continue
        data = path.read_bytes()
        if path.name != sha(data):
            bad_files.append(path.name)
        raw_shas.add(path.name)
    receipt_checks = {"n_receipts": len(receipts), "status_counts": {}, "http_200": 0,
                      "verified_success_bodies": 0, "missing_success_bodies": [],
                      "bad_hash_or_length": [], "raw_file_count": len(raw_paths),
                      "unreceipted_raw_files": sorted(raw_shas - {r.get("sha256") for r in receipts}),
                      "bad_raw_filename_hashes": bad_files, "url_with_multiple_success_hashes": {}}
    by_url = entry_receipt_index(receipts)
    for rec in receipts:
        status = str(rec.get("status"))
        receipt_checks["status_counts"][status] = receipt_checks["status_counts"].get(status, 0) + 1
        if rec.get("status") == 200:
            receipt_checks["http_200"] += 1
            p = RAW / str(rec.get("sha256"))
            if not p.is_file():
                receipt_checks["missing_success_bodies"].append(rec.get("url"))
                continue
            body = p.read_bytes()
            if sha(body) != rec.get("sha256") or len(body) != rec.get("bytes_len"):
                receipt_checks["bad_hash_or_length"].append(rec.get("url"))
            else:
                receipt_checks["verified_success_bodies"] += 1
    for url, items in by_url.items():
        hashes = sorted({r.get("sha256") for r in items if r.get("status") == 200})
        if len(hashes) > 1:
            receipt_checks["url_with_multiple_success_hashes"][url] = hashes

    new_labs = [row for entry in new.values() if entry.get("status") == "ok" for row in entry["labels"]]
    old_labs_doc = load(ROOT / "results/e421/labeled_residues_with_plddt.json")
    old_labs = old_labs_doc if isinstance(old_labs_doc, list) else old_labs_doc.get("rows", [])
    calc_new, calc_old = summarize(new), summarize(old, old_labs)
    transitions = collections.Counter((old.get(k, {}).get("status"), v.get("status")) for k, v in new.items())
    transition_doc = {f"{a} -> {b}": n for (a, b), n in sorted(transitions.items(), key=str)}

    # Exercise status changes and calculate exact residue rows for representative successful transitions.
    qualifying = {row["pdb_id"]: row for row in load(QUAL_PATH)}
    sample_ids = ["10af", "10ah", "10qf", "10cy", "10dv", "11ns", "11ny", "11pu", "10pa", "10ps"]
    spot_checks = {}
    for pdb in sample_ids:
        acc = new[pdb]["accession"]
        check = independent_entry(pdb, acc, qualifying[pdb]["cif_url"], by_url,
                                  compare_labels=new[pdb]["status"] == "ok")
        target = new[pdb]
        check_result = {k: v for k, v in check.items() if k != "labels"}
        check_result["matches_status"] = check.get("status") == target.get("status")
        if check.get("status") == "ok" and target.get("status") == "ok":
            expected = {row["uniprot_pos"]: row for row in target["labels"]}
            got = {row["uniprot_pos"]: row for row in check["labels"]}
            fields = ("e", "f", "cell", "lddt", "dist", "b", "z", "plddt")
            mismatches = []
            if expected.keys() != got.keys():
                mismatches.append({"row_keys": [len(expected), len(got)]})
            for pos in expected.keys() & got.keys():
                for field in fields:
                    a, b = expected[pos].get(field), got[pos].get(field)
                    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                        if not (math.isnan(float(a)) and math.isnan(float(b))) and abs(float(a) - float(b)) > 1e-9:
                            mismatches.append({"pos": pos, "field": field, "expected": a, "got": b})
                    elif a != b:
                        mismatches.append({"pos": pos, "field": field, "expected": a, "got": b})
            check_result["row_count"] = len(got)
            check_result["label_mismatch_count"] = len(mismatches)
            check_result["first_label_mismatches"] = mismatches[:5]
            check_result["matches_labels"] = not mismatches
        spot_checks[pdb] = check_result

    summary_fields = ("entries", "status_counts", "ok_entries", "ok_accessions", "residue_rows",
                      "cells_rows", "cells_unique", "gradient", "spearman_plddt_lddt", "n_spearman")
    def same_metric(a, b):
        if isinstance(a, float) or isinstance(b, float):
            return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
        return a == b
    new_summary_match = {field: same_metric(calc_new[field], declared["corrected"][field])
                         for field in summary_fields}
    old_summary_match = {field: same_metric(calc_old[field], declared["original"][field])
                         for field in summary_fields}

    # Count accepted alternate-location collisions across all selected chain-A segments.
    duplicate_report = {"entries_scanned": 0, "selected_asym_segment_count": 0,
                        "entries_with_duplicate_label_seq_ca": [], "duplicate_group_count": 0}
    for pdb, result in new.items():
        mb = raw_for_url(f"https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb}", by_url)
        if mb is None:
            continue
        obj = json.loads(mb)
        access = result["accession"]
        segs = [m for m in obj[pdb]["UniProt"].get(access, {}).get("mappings", [])
                if m.get("chain_id") == "A"]
        if not segs:
            continue
        rb = raw_for_url(f"https://www.ebi.ac.uk/pdbe/entry-files/{pdb}.cif", by_url)
        if rb is None:
            continue
        duplicate_report["entries_scanned"] += 1
        entry_dups = 0
        text = rb.decode("utf-8", "replace")
        for seg in segs:
            duplicate_report["selected_asym_segment_count"] += 1
            _atoms, _seq, dups = parse_ref(text, seg["struct_asym_id"], seg["entity_id"])
            entry_dups += len(dups)
        if entry_dups:
            duplicate_report["entries_with_duplicate_label_seq_ca"].append({"entry": pdb, "groups": entry_dups})
            duplicate_report["duplicate_group_count"] += entry_dups

    out = {
        "method": "independent parser and recomputation; no e427 or moluq imports",
        "receipt_checks": receipt_checks,
        "recomputed_corrected_summary": calc_new,
        "declared_corrected_summary_match": new_summary_match,
        "recomputed_original_summary": calc_old,
        "declared_original_summary_match": old_summary_match,
        "recomputed_status_transitions": transition_doc,
        "declared_status_transitions_match": transition_doc == declared["status_transitions"],
        "spot_checks": spot_checks,
        "accepted_altloc_duplicate_scan": duplicate_report,
    }
    (HERE / "independent_e421_audit.json").write_text(json.dumps(out, indent=2, allow_nan=True) + "\n")
    print(json.dumps({"receipt_checks": receipt_checks,
                      "new_summary_matches": new_summary_match,
                      "old_summary_matches": old_summary_match,
                      "transitions_match": transition_doc == declared["status_transitions"],
                      "spot_checks": spot_checks,
                      "altloc_scan": duplicate_report}, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
