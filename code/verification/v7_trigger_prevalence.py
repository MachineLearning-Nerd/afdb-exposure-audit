"""V7: trigger prevalence of the e422 numbering-fallback defect (label-free core).

See REGISTRATION.md, Addendum A4 (committed de0e7a32 before this script ran).

Core (registered, label-free): over the PDBe SIFTS UniProt mapping responses
that the original e422 census fetched (receipts in the original batch 1-4
checkpoints, bytes read from results/e422/raw and SHA-256 verified), count the
mapped segments for the submitted roster (accessions with >= 1 row in the
original cumulative ledger) and the segments whose start has a null
author_residue_number. For those segments, report the label->author offset
(auth_seq_id - label_seq_id) at the first observed residue of the segment in
the reference mmCIF (_atom_site, label_asym_id == struct_asym_id, first model).

Link to F-scope (requested by the caller, NOT part of the A4 V7 text; see the
`deviations` list): joins triggered segments to per-chain-segment median lDDT
changes between the original and corrected batch-4 cumulative ledgers. This
section reads only entry/chain/accession/lddt; pLDDT is dropped on load and is
never used.

Usage: python -I v7_trigger_prevalence.py <repo_root> <out_json>
"""
import hashlib
import json
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

MAPPING_ROLE = "E422_PDBE_MAPPINGS_V2"
CIF_ROLE = "E422_PDBE_ENTRY_FILES_V1"
DELTA = 0.2
TOKEN = re.compile(r"""'(?:[^']|'(?!\s|$))*'|"(?:[^"]|"(?!\s|$))*"|\S+""")


def load_rows(path):
    """Ledger rows with the pLDDT field dropped on load (registered blinding)."""
    d = json.loads(path.read_text())
    rows = []
    for r in d["rows"]:
        r = dict(r)
        r.pop("plddt", None)
        rows.append(r)
    return d, rows


def read_raw(raw_dir, rec):
    b = (raw_dir / rec["raw_name"]).read_bytes()
    h = hashlib.sha256(b).hexdigest()
    if h != rec["sha256"] or h != rec["raw_name"]:
        raise RuntimeError(f"SHA-256 mismatch for {rec['url']}: {h}")
    return b


def atom_site_first_model(text):
    """Return list of dicts (label_asym_id, label_seq_id, auth_seq_id, ins) for model 1."""
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "loop_" and i + 1 < n and lines[i + 1].startswith("_atom_site."):
            break
        i += 1
    if i >= n:
        raise RuntimeError("no _atom_site loop")
    i += 1
    fields = []
    while i < n and lines[i].startswith("_atom_site."):
        fields.append(lines[i].split()[0][len("_atom_site."):])
        i += 1
    idx = {f: k for k, f in enumerate(fields)}
    toks = []
    while i < n:
        ln = lines[i]
        if ln.startswith("_") or ln.startswith("loop_") or ln.startswith("#") or ln.startswith("data_"):
            break
        toks.extend(TOKEN.findall(ln))
        i += 1
    nf = len(fields)
    if len(toks) % nf:
        raise RuntimeError(f"_atom_site token count {len(toks)} not a multiple of {nf}")
    out, first_model = [], None
    for k in range(0, len(toks), nf):
        rec = toks[k:k + nf]
        model = rec[idx["pdbx_PDB_model_num"]] if "pdbx_PDB_model_num" in idx else "1"
        if first_model is None:
            first_model = model
        if model != first_model:
            continue
        ins = rec[idx["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in idx else "?"
        out.append((rec[idx["label_asym_id"]], rec[idx["label_seq_id"]], rec[idx["auth_seq_id"]],
                    ins if ins not in ("?", ".") else ""))
    return out


def to_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def first_observed_offset(atoms, asym, lab_s, lab_e):
    """Offset auth - label at the first observed residue of [lab_s, lab_e] in asym."""
    res = {}
    for a_asym, lab, auth, ins in atoms:
        if a_asym != asym:
            continue
        li = to_int(lab)
        if li is None or not (lab_s <= li <= lab_e):
            continue
        res.setdefault(li, (auth, ins))
    if not res:
        return {"status": "NO_OBSERVED_RESIDUE_IN_SEGMENT"}
    first = min(res)
    auth, ins = res[first]
    ai = to_int(auth)
    if ai is None:
        return {"status": "NON_INTEGER_AUTH_SEQ_ID", "first_label_seq_id": first, "auth_seq_id": auth}
    offs = {to_int(a) - li for li, (a, _) in res.items() if to_int(a) is not None}
    return {"status": "DETERMINED", "first_label_seq_id": first, "first_auth_seq_id": ai,
            "first_ins_code": ins, "offset": ai - first, "n_observed_residues": len(res),
            "offset_constant_across_segment": len(offs) == 1}


def median_by_key(rows):
    acc = defaultdict(list)
    for r in rows:
        y = r.get("lddt")
        if y is None or y != y:
            continue
        acc[(r["entry"], r["chain"], r["accession"])].append(float(y))
    return {k: statistics.median(v) for k, v in acc.items()}


def main(root, out_path):
    root = Path(root)
    raw_dir = root / "results/e422/raw"
    deviations = []
    orig_rows, receipts = [], []
    for b in (1, 2, 3, 4):
        d, rows = load_rows(root / f"results/e422/batches/e422_batch0{b}_checkpoint.json")
        orig_rows += rows
        for rec in d["receipts"]:
            receipts.append((b, rec))
    roster = sorted({r["accession"] for r in orig_rows})

    # --- receipts: mapping responses and reference CIFs (dedupe by URL) ---
    def collect(role):
        by_url, dup, conflict = {}, 0, []
        for b, rec in receipts:
            if rec["role"] != role or rec.get("status") != 200 or not rec.get("raw_name"):
                continue
            if rec["url"] in by_url:
                dup += 1
                if by_url[rec["url"]][1]["sha256"] != rec["sha256"]:
                    conflict.append(rec["url"])
                continue
            by_url[rec["url"]] = (b, rec)
        return by_url, dup, conflict

    maps, map_dup, map_conf = collect(MAPPING_ROLE)
    cifs, cif_dup, cif_conf = collect(CIF_ROLE)
    cif_by_entry = {u.rsplit("/", 1)[1].removesuffix(".cif").lower(): v for u, v in cifs.items()}
    n_non200 = sum(1 for _, r in receipts if r["role"] == MAPPING_ROLE and r.get("status") != 200)

    segments = []
    for url, (b, rec) in sorted(maps.items()):
        entry = url.rsplit("/", 1)[1].lower()
        doc = json.loads(read_raw(raw_dir, rec).decode("utf-8"))
        uni = (doc.get(entry) or {}).get("UniProt") or {}
        for acc in sorted(uni):
            if acc not in roster:
                continue
            for k, m in enumerate(uni[acc].get("mappings", [])):
                st, en = m["start"], m["end"]
                segments.append({
                    "entry": entry, "batch": b, "accession": acc, "segment_index": k,
                    "chain_id": m.get("chain_id"), "struct_asym_id": m.get("struct_asym_id"),
                    "entity_id": m.get("entity_id"), "unp_start": m["unp_start"], "unp_end": m["unp_end"],
                    "start_label": st["residue_number"], "start_author": st["author_residue_number"],
                    "end_label": en["residue_number"], "end_author": en["author_residue_number"],
                    "start_author_null": st["author_residue_number"] is None,
                    "start_author_zero": st["author_residue_number"] == 0,
                    "end_author_null": en["author_residue_number"] is None,
                })

    # --- offsets at the first observed residue for triggered segments ---
    atoms_cache = {}
    for s in segments:
        if not s["start_author_null"]:
            continue
        e = s["entry"]
        if e not in cif_by_entry:
            s["offset_check"] = {"status": "NO_REFERENCE_CIF_RECEIPT"}
            continue
        if e not in atoms_cache:
            atoms_cache.clear()  # segments are grouped by entry; keep one parsed CIF in memory
            _, rec = cif_by_entry[e]
            atoms_cache[e] = atom_site_first_model(read_raw(raw_dir, rec).decode("utf-8", "replace"))
        s["offset_check"] = first_observed_offset(atoms_cache[e], s["struct_asym_id"],
                                                  s["start_label"], s["end_label"])

    trig = [s for s in segments if s["start_author_null"]]
    det = [s for s in trig if s["offset_check"]["status"] == "DETERMINED"]
    offs = Counter(s["offset_check"]["offset"] for s in det)
    nonzero = [s for s in det if s["offset_check"]["offset"] != 0]

    per_acc = {}
    for a in roster:
        ss = [s for s in segments if s["accession"] == a]
        tt = [s for s in ss if s["start_author_null"]]
        dd = [s for s in tt if s["offset_check"]["status"] == "DETERMINED"]
        per_acc[a] = {"segments": len(ss), "entries": len({s["entry"] for s in ss}),
                      "triggered": len(tt),
                      "triggered_offset_zero": sum(s["offset_check"]["offset"] == 0 for s in dd),
                      "triggered_offset_nonzero": sum(s["offset_check"]["offset"] != 0 for s in dd),
                      "triggered_offset_undetermined": len(tt) - len(dd)}
    acc_with_seg = [a for a in roster if per_acc[a]["segments"]]

    core = {
        "roster_accessions": len(roster),
        "roster_accessions_with_mapped_segments": len(acc_with_seg),
        "mapping_responses_used": len(maps),
        "mapping_receipts_duplicate_urls": map_dup,
        "mapping_receipts_sha_conflicts": map_conf,
        "mapping_receipts_non200": n_non200,
        "reference_cif_receipts": len(cifs),
        "reference_cif_duplicate_urls": cif_dup,
        "reference_cif_sha_conflicts": cif_conf,
        "entries_with_roster_segments": len({s["entry"] for s in segments}),
        "segments": {
            "mapped": len(segments),
            "start_author_null_TRIGGER": len(trig),
            "start_author_null_share": len(trig) / len(segments) if segments else None,
            "offset_status_counts": dict(Counter(s["offset_check"]["status"] for s in trig)),
            "offset_distribution": {str(k): v for k, v in sorted(offs.items())},
            "offset_zero": sum(v for k, v in offs.items() if k == 0),
            "offset_nonzero": sum(v for k, v in offs.items() if k != 0),
            "offset_constant_across_segment_false": sum(
                not s["offset_check"]["offset_constant_across_segment"] for s in det),
            "supplementary_start_author_zero_also_falsy": sum(s["start_author_zero"] for s in segments),
            "supplementary_end_author_null": sum(s["end_author_null"] for s in segments),
        },
        "accessions": {
            "with_any_trigger": sum(per_acc[a]["triggered"] > 0 for a in roster),
            "with_nonzero_offset_trigger": sum(per_acc[a]["triggered_offset_nonzero"] > 0 for a in roster),
            "with_only_zero_offset_triggers": sum(per_acc[a]["triggered"] > 0
                                                  and per_acc[a]["triggered_offset_nonzero"] == 0
                                                  and per_acc[a]["triggered_offset_undetermined"] == 0
                                                  for a in roster),
        },
    }

    # --- link to F-scope (reads lddt only; plddt dropped on load) ---
    corr_rows = []
    for b in (1, 2, 3, 4):
        corr_rows += load_rows(root / f"results/e427/batches/e422_batch0{b}_checkpoint.json")[1]
    mo, mc = median_by_key(orig_rows), median_by_key(corr_rows)
    both = sorted(set(mo) & set(mc))
    damaged = {k for k in both if abs(mc[k] - mo[k]) > DELTA}
    key = lambda s: (s["entry"], s["chain_id"], s["accession"])
    trig_nz_keys = {key(s) for s in nonzero}
    trig_zero_keys = {key(s) for s in det if s["offset_check"]["offset"] == 0}
    trig_keys = {key(s) for s in trig}
    link = {
        "ledger_segment_key": "(entry, chain [= SIFTS chain_id / auth chain], accession)",
        "segments_in_both_ledgers": len(both),
        "segments_original_only": len(set(mo) - set(mc)),
        "segments_corrected_only": len(set(mc) - set(mo)),
        "segments_abs_delta_gt_0.2": len(damaged),
        "accessions_with_abs_delta_gt_0.2": len({k[2] for k in damaged}),
        "triggered_nonzero_offset_segments": len(nonzero),
        "triggered_nonzero_offset_segments_with_ledger_key_in_both": sum(key(s) in set(both) for s in nonzero),
        "triggered_nonzero_offset_segments_in_damaged_set": sum(key(s) in damaged for s in nonzero),
        "triggered_zero_offset_segments_in_damaged_set": sum(
            key(s) in damaged for s in det if s["offset_check"]["offset"] == 0),
        "triggered_zero_offset_in_damaged_set_with_offset_changing_within_segment": sum(
            key(s) in damaged and not s["offset_check"]["offset_constant_across_segment"]
            for s in det if s["offset_check"]["offset"] == 0),
        "damaged_keys_with_a_triggered_nonzero_segment": len(damaged & trig_nz_keys),
        "damaged_keys_with_only_zero_offset_triggers": len((damaged & trig_zero_keys) - trig_nz_keys),
        "damaged_keys_without_any_trigger": len(damaged - trig_keys),
        "note": "Mapping segments and ledger keys are not 1:1: one (entry, chain, accession) key can hold "
                "several SIFTS segments; an untriggered key can be absent from the ledger (e.g. failed triple).",
    }

    out = {
        "schema": "v7-trigger-prevalence-v1",
        "registration": "paper/arxiv/verification/REGISTRATION.md Addendum A4 (commit de0e7a32)",
        "trigger_definition": "SIFTS mapping start.author_residue_number is null",
        "inputs": "results/e422/batches/e422_batch0{1..4}_checkpoint.json receipts; bytes from results/e422/raw "
                  "(SHA-256 verified against receipt and raw_name); link section also reads "
                  "results/e427/batches rows (lddt only)",
        "core_label_free": core,
        "per_accession": per_acc,
        "triggered_segments": [
            {k: s[k] for k in ("entry", "accession", "chain_id", "struct_asym_id", "unp_start", "unp_end",
                               "start_label", "end_label", "end_author")} | {"offset_check": s["offset_check"]}
            for s in trig],
        "link_to_fscope_NOT_IN_A4_V7_TEXT": link,
        "deviations": [
            "Trigger->damage link (link_to_fscope_NOT_IN_A4_V7_TEXT) reads lDDT, which the A4 V7 text says V7 "
            "does not read; it was requested by the caller after registration and is reported as a separate, "
            "non-pre-specified section. The registered core (core_label_free/per_accession) reads no labels.",
            "Supplementary counts not named in A4 are added: start author_residue_number == 0 (also falsy for "
            "the old `or` fallback) and null end author_residue_number; the offset constancy flag.",
        ],
    }
    Path(out_path).write_text(json.dumps(out, indent=1))
    print(json.dumps({"core": core, "link": link}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
