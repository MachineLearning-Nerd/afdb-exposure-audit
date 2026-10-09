"""V5: template channel, measured (REGISTRATION.md Addendum A3).

Reads the AFDB model files the e422 census fetched (results/e422/raw, via the
batch-checkpoint receipts; no re-download), records software/version metadata and
every `_ma_template_ref_db_details` template with db_name = PDB, then fetches
PDBe SIFTS UniProt mappings and RCSB initial release dates for each template
(V2's cached, hash-receipted Fetcher; cache under cache/v5/).
Only the `accession` field of batch rows is read (labels and pLDDT are not read).
Usage: python -I v5_templates.py <repo_root> <out_json>
"""
import collections
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from v2_prior_exposure import Fetcher  # noqa: E402  (module has no import-time side effects)

SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{id}"
RCSB = "https://data.rcsb.org/rest/v1/core/entry/{ID}"
TEMPLATE_LIMIT = "2021-02-15"
BATCHES = (1, 2, 3, 4)
AFDB_RE = re.compile(r"^https://alphafold\.ebi\.ac\.uk/files/AF-(.+?)-F(\d+)-model_v(\d+)\.cif$")
TOK = re.compile(r"""'(?:[^']|'(?=\S))*'|"(?:[^"]|"(?=\S))*"|\S+""")


def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parse_cif(text):
    """Minimal mmCIF reader -> {category: [ {item: value} ]} (single data block)."""
    toks = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(";"):
            buf = [ln[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            toks.append(("V", "\n".join(buf).strip()))
            i += 1
            continue
        s = ln
        if s.lstrip().startswith("#"):
            i += 1
            continue
        for t in TOK.findall(s):
            if t[0] in "'\"" and len(t) >= 2 and t[-1] == t[0]:
                toks.append(("V", t[1:-1]))
            else:
                toks.append(("T", t))
        i += 1
    out = collections.defaultdict(list)
    j = 0
    while j < len(toks):
        kind, t = toks[j]
        if kind == "T" and t == "loop_":
            j += 1
            names = []
            while j < len(toks) and toks[j][0] == "T" and toks[j][1].startswith("_"):
                names.append(toks[j][1])
                j += 1
            vals = []
            while j < len(toks) and not (toks[j][0] == "T" and (toks[j][1].startswith("_") or toks[j][1] == "loop_"
                                                                   or toks[j][1].startswith("data_"))):
                vals.append(toks[j][1])
                j += 1
            cat = names[0].split(".")[0]
            for k in range(0, len(vals) - len(vals) % len(names), len(names)):
                out[cat].append({n.split(".", 1)[1]: vals[k + m] for m, n in enumerate(names)})
            continue
        if kind == "T" and t.startswith("_") and "." in t and j + 1 < len(toks):
            cat, item = t.split(".", 1)
            if not out[cat]:
                out[cat].append({})
            out[cat][0][item] = toks[j + 1][1]
            j += 2
            continue
        j += 1
    return out


def roster(root):
    accs, cifs, inputs = set(), collections.defaultdict(list), {}
    for b in BATCHES:
        p = root / f"results/e422/batches/e422_batch0{b}_checkpoint.json"
        inputs[str(p.relative_to(root))] = sha256_file(p)
        d = json.loads(p.read_text())
        accs |= {r["accession"] for r in d["rows"]}  # accession field only
        for r in d["receipts"]:
            m = AFDB_RE.match(r["url"])
            if m:
                cifs[m.group(1)].append({"batch": b, "url": r["url"], "raw_name": r["raw_name"],
                                         "receipt_sha256": r["sha256"], "status": r["status"],
                                         "utc": r.get("utc")})
    return sorted(accs), cifs, inputs


def sifts(f, pid):
    st, txt = f.get("GET", SIFTS.format(id=pid.lower()))
    if st != 200:
        return {"ok": False, "status": st}
    try:
        d = json.loads(txt)
        u = d.get(pid.lower(), {}).get("UniProt", {})
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "status": st, "error": repr(e)}
    chains = {acc: sorted({m.get("chain_id") for m in v.get("mappings", [])}) for acc, v in u.items()}
    return {"ok": True, "status": st, "uniprot": sorted(u), "chains": chains}


def rcsb(f, pid):
    st, txt = f.get("GET", RCSB.format(ID=pid.upper()))
    if st != 200:
        return {"ok": False, "status": st}
    try:
        info = json.loads(txt).get("rcsb_accession_info", {})
        rel = info.get("initial_release_date")
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "status": st, "error": repr(e)}
    if not rel:
        return {"ok": False, "status": st, "error": "no initial_release_date"}
    return {"ok": True, "status": st, "initial_release_date": rel[:10],
            "deposit_date": (info.get("deposit_date") or "")[:10] or None}


def same(acc, key):
    return key == acc or key.split("-")[0] == acc


def main(root, out):
    root = Path(root)
    accs, cifs, inputs = roster(root)
    deviations, notes = [], []
    assert len(accs) == 138, len(accs)
    p = root / "paper/arxiv/verification/out/v2_prior_exposure.json"
    inputs[str(p.relative_to(root))] = sha256_file(p)
    v2 = json.loads(p.read_text())["accessions"]

    per_acc, file_checks, sw_counter, mg_counter = {}, [], collections.Counter(), collections.Counter()
    for acc in accs:
        recs = cifs.get(acc, [])
        names = sorted({r["raw_name"] for r in recs})
        rec = {"cif_receipts": recs, "raw_names": names}
        if len(names) != 1:
            rec["error"] = f"expected exactly one AFDB model raw file, found {len(names)}"
            per_acc[acc] = rec
            continue
        rp = root / "results/e422/raw" / names[0]
        got = sha256_file(rp) if rp.exists() else None
        ok = got is not None and all(got == r["receipt_sha256"] for r in recs)
        file_checks.append({"accession": acc, "raw_name": names[0], "sha256": got, "verified": ok})
        rec["sha256_verified"] = ok
        if not ok:
            rec["error"] = "raw file missing or sha256 mismatch"
            per_acc[acc] = rec
            continue
        cif = parse_cif(rp.read_text())
        tgt = (cif.get("_ma_target_ref_db_details") or [{}])[0].get("db_accession")
        rec["target_db_accession"] = tgt
        rec["software"] = [{k: s.get(k) for k in ("name", "version", "type", "description")}
                           for s in cif.get("_software", [])]
        rec["model_group_name"] = [m.get("model_group_name") for m in cif.get("_ma_model_list", [])]
        rec["protocol_steps"] = [s.get("method_type") for s in cif.get("_ma_protocol_step", [])]
        rec["ma_software_group"] = cif.get("_ma_software_group", [])
        revs = cif.get("_pdbx_audit_revision_history", [])
        rec["revision_dates"] = [r.get("revision_date") for r in revs]
        tdet = {t.get("template_id"): t for t in cif.get("_ma_template_details", [])}
        tmpl, other_db = [], []
        for t in cif.get("_ma_template_ref_db_details", []):
            if t.get("db_name") == "PDB":
                tmpl.append({"template_id": t.get("template_id"), "pdb_id": t.get("db_accession_code").upper(),
                             "template_auth_asym_id": tdet.get(t.get("template_id"), {}).get("template_auth_asym_id")})
            else:
                other_db.append(t)
        rec["templates"] = tmpl
        rec["non_pdb_template_refs"] = other_db
        for s in rec["software"]:
            sw_counter[f"{s['name']} {s['version']}"] += 1
        for m in rec["model_group_name"]:
            mg_counter[m] += 1
        per_acc[acc] = rec

    # template lookups (V2 Fetcher, sequential, DELAY 0.3 s => < 4 req/s)
    f = Fetcher(HERE / "cache/v5")
    pids = sorted({t["pdb_id"] for r in per_acc.values() for t in r.get("templates", [])})
    look = {}
    for i, pid in enumerate(pids):
        look[pid] = {"sifts": sifts(f, pid), "rcsb": rcsb(f, pid)}
        print(i + 1, len(pids), pid, look[pid]["sifts"]["status"], look[pid]["rcsb"].get("initial_release_date"),
              flush=True)

    classes, classes_strict, chain_cls = {}, {}, {}
    for acc, r in per_acc.items():
        if "error" in r:
            classes[acc] = classes_strict[acc] = chain_cls[acc] = "UNRESOLVED"
            continue
        ts = r["templates"]
        hit = hit_chain = failed = failed_any = False
        for t in ts:
            L = look[t["pdb_id"]]
            s = L["sifts"]
            t["sifts_ok"], t["rcsb_ok"] = s["ok"], L["rcsb"]["ok"]
            t["initial_release_date"] = L["rcsb"].get("initial_release_date")
            if not s["ok"]:
                failed = failed_any = True
                t["maps_to_accession"] = None
                continue
            if not L["rcsb"]["ok"]:
                failed_any = True
            keys = [k for k in s["uniprot"] if same(acc, k)]
            t["maps_to_accession"] = bool(keys)
            t["isoform_only_match"] = bool(keys) and acc not in keys
            ch = set()
            for k in keys:
                ch |= set(s["chains"].get(k, []))
            t["template_chain_maps_to_accession"] = t["template_auth_asym_id"] in ch if keys else False
            hit |= bool(keys)
            hit_chain |= t["template_chain_maps_to_accession"]
        r["latest_template_release"] = max((t["initial_release_date"] for t in ts if t.get("initial_release_date")),
                                           default=None)

        def cls(h, fl):
            if not ts:
                return "NONE"
            if h:
                return "SAME_ACC"
            return "UNRESOLVED" if fl else "OTHER"
        classes[acc] = cls(hit, failed)
        classes_strict[acc] = cls(hit, failed_any)
        chain_cls[acc] = cls(hit_chain, failed)
        r["template_class"] = classes[acc]

    v2cls = {a: v2[a]["training"]["class"] if a in v2 else "NOT_IN_V2" for a in accs}
    xt = collections.defaultdict(collections.Counter)
    for a in accs:
        xt[classes[a]][v2cls[a]] += 1
    rel = {pid: L["rcsb"].get("initial_release_date") for pid, L in look.items()}
    dated = {k: v for k, v in rel.items() if v}
    latest = max(dated.values(), default=None)
    after = sorted([{"pdb_id": k, "initial_release_date": v,
                     "accessions": sorted(a for a, r in per_acc.items()
                                          if any(t["pdb_id"] == k for t in r.get("templates", [])))}
                    for k, v in dated.items() if v > TEMPLATE_LIMIT], key=lambda x: x["initial_release_date"])
    fail_sifts = sorted(k for k, L in look.items() if not L["sifts"]["ok"])
    fail_rcsb = sorted(k for k, L in look.items() if not L["rcsb"]["ok"])
    ntemp = [len(r.get("templates", [])) for r in per_acc.values()]
    tgt_mismatch = sorted(a for a, r in per_acc.items() if r.get("target_db_accession") not in (None, a))
    if tgt_mismatch:
        notes.append(f"_ma_target_ref_db_details.db_accession differs from roster accession for {tgt_mismatch}")
    notes.append("Template class uses entry-level SIFTS mapping as registered; a template counts as SAME_ACC if any "
                 "SIFTS UniProt key equals the accession or is an isoform of it (key.split('-')[0]). "
                 "Isoform-only matches are flagged per template.")
    notes.append("UNRESOLVED (primary) = a SIFTS lookup failed and no SAME_ACC; the RCSB release-date lookup does not "
                 "determine the class. `template_class_counts_strict` also treats an RCSB failure as a failed lookup.")
    notes.append("`template_class_counts_chain_level_NOT_PRESPECIFIED`: SAME_ACC only if the template's "
                 "_ma_template_details.template_auth_asym_id is a SIFTS chain_id mapped to the accession. Not registered.")
    notes.append("Batch rows: only the `accession` field is read; CIF pLDDT (_ma_qa_metric_local) is not parsed.")
    notes.append("Any non-200 SIFTS response (including HTTP 404, which PDBe also returns for an entry with no "
                 "UniProt mapping) is typed as a failed lookup, never as 'maps to nothing'; the RCSB status of each "
                 "such entry is in template_lookups. V2's Fetcher caches only 200/204 bodies, so non-200 responses "
                 "are re-requested on rerun (their receipts are kept in `receipts`).")

    res = {
        "schema": "paper-verification-v5-templates",
        "registration": "REGISTRATION.md Addendum A3",
        "utc": datetime.now(timezone.utc).isoformat(),
        "inputs_sha256": inputs,
        "raw_file_checks": {"n": len(file_checks), "n_verified": sum(c["verified"] for c in file_checks),
                            "files": file_checks},
        "n_accessions": len(accs),
        "results": {
            "template_class_counts": dict(collections.Counter(classes.values())),
            "template_class_counts_strict": dict(collections.Counter(classes_strict.values())),
            "template_class_counts_chain_level_NOT_PRESPECIFIED": dict(collections.Counter(chain_cls.values())),
            "crosstab_template_class_x_v2_training_class": {k: dict(v) for k, v in sorted(xt.items())},
            "n_unique_template_pdb_ids": len(pids),
            "templates_per_accession": {"min": min(ntemp), "max": max(ntemp),
                                        "mean": sum(ntemp) / len(ntemp), "n_with_zero": ntemp.count(0)},
            "latest_template_release_overall": latest,
            "earliest_template_release_overall": min(dated.values(), default=None),
            "documented_template_limit": TEMPLATE_LIMIT,
            "n_templates_released_after_limit": len(after),
            "templates_released_after_limit": after,
            "software_versions": dict(sw_counter),
            "model_group_names": dict(mg_counter),
            "sifts_failed_pdb_ids": fail_sifts,
            "rcsb_failed_pdb_ids": fail_rcsb,
            "unresolved_accessions": sorted(a for a, c in classes.items() if c == "UNRESOLVED"),
            "per_accession": {a: per_acc[a] | {"v2_training_class": v2cls[a]} for a in accs},
            "template_lookups": look,
        },
        "receipts_summary": {"n": len(f.receipts), "n_cached": sum(r["cached"] for r in f.receipts),
                             "status_counts": dict(collections.Counter(str(r["status"]) for r in f.receipts))},
        "receipts": f.receipts,
        "deviations": deviations,
        "notes": notes,
    }
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: res["results"][k] for k in ("template_class_counts", "template_class_counts_strict",
                                                      "crosstab_template_class_x_v2_training_class",
                                                      "latest_template_release_overall",
                                                      "n_templates_released_after_limit", "software_versions",
                                                      "sifts_failed_pdb_ids", "rcsb_failed_pdb_ids")}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
