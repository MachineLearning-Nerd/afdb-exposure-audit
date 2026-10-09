"""e428: exposure audit of two published post-cutoff evaluations.

See docs/e428_registration.md (registered and timestamped before this script
was run). T = Terwilliger et al. 2024 (102 structures, published Ca r.m.s.d.);
A = AlphaFlow ATLAS test split (82 chains). Every HTTP response is cached under
data/e428/http/<sha256(request)>.json; receipts are written with the results.

Usage: python -I experiments/e428_exposure_audit.py <repo_root>
"""
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

TRAIN_CUTOFF = "2018-04-30"
ESMFOLD_CUTOFF = "2020-05-01"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
RCSB_ENTITY = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}"
SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
UA = {"User-Agent": "moluq-e428/1 (research; public endpoints)"}
DELAY = 0.3
XLSX_SHA = "7b3e1473b99e1575019d8c095d292d6a96de702c36580241cf85932e7d69d143"


class Fetcher:
    def __init__(self, cache):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.receipts = []

    def get(self, method, url, body=None):
        k = hashlib.sha256(f"{method}\n{url}\n{body or ''}".encode()).hexdigest()
        p = self.cache / f"{k}.json"
        if p.exists():
            rec = json.loads(p.read_text())
            self.receipts.append({x: rec[x] for x in ("url", "status", "sha256", "utc")} | {"cached": True})
            return rec["status"], rec["body"]
        status, text = None, None
        for attempt in range(3):
            try:
                hdr = dict(UA)
                if body:
                    hdr["Content-Type"] = "application/json"
                req = urllib.request.Request(url, data=body.encode() if body else None, headers=hdr, method=method)
                with urllib.request.urlopen(req, timeout=90) as r:
                    status, text = r.status, r.read().decode()
                break
            except urllib.error.HTTPError as e:
                status, text = e.code, e.read().decode(errors="replace")
                if e.code < 500 and e.code != 429:
                    break
            except Exception as e:  # network error -> retry
                status, text = -1, repr(e)
            time.sleep(DELAY * 2 ** (attempt + 1))
        time.sleep(DELAY)
        rec = {"method": method, "url": url, "request_body": body, "status": status, "body": text,
               "sha256": hashlib.sha256((text or "").encode()).hexdigest(),
               "utc": datetime.now(timezone.utc).isoformat()}
        if status in (200, 204, 404):
            p.write_text(json.dumps(rec))
        self.receipts.append({x: rec[x] for x in ("url", "status", "sha256", "utc")} | {"cached": False})
        return status, text


# ----------------------------------------------------------------------------- RCSB search
def search_ids(f, nodes):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes} if len(nodes) > 1 else nodes[0],
         "return_type": "entry", "request_options": {"return_all_hits": True, "results_verbosity": "compact"}}
    st, txt = f.get("POST", SEARCH, json.dumps(q, sort_keys=True))
    if st == 204:
        return []
    if st != 200:
        return None
    return sorted(json.loads(txt).get("result_set", []))


def date_node(d):
    return {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_accession_info.initial_release_date", "operator": "less_or_equal",
        "value": f"{d}T23:59:59Z"}}


def acc_node(acc):
    return {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}}


def seq_node(seq, ident):
    return {"type": "terminal", "service": "sequence", "parameters": {
        "evalue_cutoff": 0.1, "identity_cutoff": ident, "sequence_type": "protein", "value": seq}}


def exposure(f, accs, seq):
    """Registered classes: SAME_ACCESSION > SEQ95 > SEQ30 > NOVEL, entries released <= 2018-04-30."""
    out = {"accessions": accs, "no_accession": not accs}
    same = []
    failed = False
    for a in accs:
        ids = search_ids(f, [acc_node(a), date_node(TRAIN_CUTOFF)])
        if ids is None:
            failed = True
        else:
            same += ids
    out["same_accession"] = len(set(same))
    for lab, ident in (("seq95", 0.95), ("seq30", 0.30)):
        ids = search_ids(f, [seq_node(seq, ident), date_node(TRAIN_CUTOFF)])
        if ids is None:
            failed = True
        out[lab] = None if ids is None else len(ids)
        out[lab + "_ids"] = (ids or [])[:20]
    if failed:
        out["class"] = "QUERY_FAILED"
    elif out["same_accession"] > 0:
        out["class"] = "SAME_ACCESSION"
    elif out["seq95"] > 0:
        out["class"] = "SEQ95"
    elif out["seq30"] > 0:
        out["class"] = "SEQ30"
    else:
        out["class"] = "NOVEL"
    return out


def first_protein_entity(f, pdb):
    st, txt = f.get("GET", RCSB_ENTRY.format(pdb=pdb))
    if st != 200:
        return None
    entry = json.loads(txt)
    eids = sorted(int(e) for e in entry["rcsb_entry_container_identifiers"].get("polymer_entity_ids", []))
    for eid in eids:
        st, txt = f.get("GET", RCSB_ENTITY.format(pdb=pdb, eid=eid))
        if st != 200:
            return None
        ent = json.loads(txt)
        if ent.get("entity_poly", {}).get("type") != "polypeptide(L)":
            continue
        refs = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers") or []
        accs = sorted({r["database_accession"] for r in refs if r.get("database_name") == "UniProt"})
        return {"entity_id": eid, "seq": ent["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", ""),
                "accessions": accs, "release": entry["rcsb_accession_info"]["initial_release_date"][:10]}
    return None


def sifts_accessions(f, pdb, chain):
    st, txt = f.get("GET", SIFTS.format(pdb=pdb.lower()))
    if st != 200:
        return None
    d = json.loads(txt).get(pdb.lower(), {}).get("UniProt", {})
    return sorted(a for a, v in d.items() if any(m.get("chain_id") == chain for m in v.get("mappings", [])))


# ----------------------------------------------------------------------------- workbook
def read_xlsx(path):
    z = zipfile.ZipFile(path)
    wb = z.read("xl/workbook.xml").decode()
    rels = z.read("xl/_rels/workbook.xml.rels").decode()
    rid2file = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))
    rid2file.update({a: b for b, a in re.findall(r'Target="([^"]+)"[^>]*Id="(rId\d+)"', rels)})
    sheets = {n: "xl/" + rid2file[r].lstrip("/").replace("xl/", "")
              for n, r in re.findall(r'<sheet name="([^"]+)"[^>]*r:id="(rId\d+)"', wb)}
    ss = z.read("xl/sharedStrings.xml").decode()
    strs = [re.sub("<[^>]+>", "", m) for m in re.findall(r"<si>(.*?)</si>", ss, re.S)]

    def cells(name):
        t = z.read(sheets[name]).decode()
        out = {}
        for ref, attrs, v in re.findall(r'<c r="([A-Z]+\d+)"([^>]*)>(?:<f[^>]*/>|<f[^>]*>.*?</f>)?(?:<v>(.*?)</v>)?', t, re.S):
            if v == "":
                continue
            out[ref] = strs[int(v)] if 't="s"' in attrs else v
        return out
    return cells


def load_terwilliger(root):
    p = root / "data/e428/terwilliger/AlphaFoldCrystal_M.xlsx"
    assert hashlib.sha256(p.read_bytes()).hexdigest() == XLSX_SHA, "workbook hash changed"
    cells = read_xlsx(p)
    pr = cells("PRED-RMSDLIST")
    rows = []
    for r in range(8, 110):
        rows.append({"pdb": pr[f"AP{r}"].strip().upper(), "rmsd": float(pr[f"AR{r}"])})
    sm = cells("Summary")
    plddt = {}
    for r in range(21, 236):
        if f"B{r}" in sm and f"H{r}" in sm:
            try:
                plddt[sm[f"B{r}"].strip().upper()] = float(sm[f"H{r}"])
            except ValueError:
                pass
    return rows, plddt


# ----------------------------------------------------------------------------- statistics
def diff_median(a, b):
    return float(np.median(a) - np.median(b))


def boot_ci(a, b, n, seed):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a), np.asarray(b)
    d = np.median(a[rng.integers(0, a.size, (n, a.size))], axis=1) - \
        np.median(b[rng.integers(0, b.size, (n, b.size))], axis=1)
    return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]


def perm_p(a, b, n, seed):
    rng = np.random.default_rng(seed)
    x = np.concatenate([a, b])
    obs = abs(diff_median(a, b))
    hits = 0
    for _ in range(n):
        rng.shuffle(x)
        hits += abs(np.median(x[:len(a)]) - np.median(x[len(a):])) >= obs - 1e-12
    return float((hits + 1) / (n + 1))


def compare(vals, close, label, seed_b, seed_p):
    a = [v for v, c in zip(vals, close) if c]
    b = [v for v, c in zip(vals, close) if not c]
    out = {"n_close": len(a), "n_other": len(b),
           "median_close": float(np.median(a)) if a else None, "median_other": float(np.median(b)) if b else None}
    if len(a) < 10 or len(b) < 10:
        out["tested"] = False
        return out
    out.update({"tested": True, "diff_median_close_minus_other": diff_median(a, b),
                "ci95": boot_ci(a, b, 10_000, seed_b), "perm_p_two_sided": perm_p(np.array(a), np.array(b), 10_000, seed_p)})
    return out


# ----------------------------------------------------------------------------- main
def main(root):
    root = Path(root)
    f = Fetcher(root / "data/e428/http")
    out_dir = root / "results/e428"
    out_dir.mkdir(parents=True, exist_ok=True)

    # T
    rows, plddt = load_terwilliger(root)
    t = []
    for r in rows:
        ent = first_protein_entity(f, r["pdb"])
        rec = dict(r, plddt=plddt.get(r["pdb"]))
        if ent is None:
            rec["class"] = "QUERY_FAILED"
        else:
            rec.update(entity_id=ent["entity_id"], release=ent["release"], seqlen=len(ent["seq"]))
            rec.update(exposure(f, ent["accessions"], ent["seq"]))
        t.append(rec)
    ok = [r for r in t if r["class"] != "QUERY_FAILED"]
    close = [r["class"] in ("SAME_ACCESSION", "SEQ95") for r in ok]
    stats_t = {
        "n": len(t), "n_query_failed": len(t) - len(ok),
        "class_counts": {c: sum(r["class"] == c for r in t)
                         for c in ("SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL", "QUERY_FAILED")},
        "n_no_accession": sum(bool(r.get("no_accession")) for r in t),
        "P1_rmsd": compare([r["rmsd"] for r in ok], close, "rmsd", 0, 1),
        "S4_median_rmsd_by_class": {c: float(np.median([r["rmsd"] for r in ok if r["class"] == c]))
                                    for c in ("SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL")
                                    if any(r["class"] == c for r in ok)},
    }
    okp = [r for r in ok if r["plddt"] is not None]
    closep = [r["class"] in ("SAME_ACCESSION", "SEQ95") for r in okp]
    stats_t["S2_plddt"] = compare([r["plddt"] for r in okp], closep, "plddt", 3, 4)
    if len(okp) >= 20:
        y = np.log([r["rmsd"] for r in okp])
        X = np.column_stack([np.ones(len(okp)), np.array(closep, float), [r["plddt"] for r in okp]])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        rng = np.random.default_rng(2)
        bs = []
        for _ in range(10_000):
            i = rng.integers(0, len(y), len(y))
            if X[i, 1].min() == X[i, 1].max():
                continue
            bs.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0][1])
        stats_t["S3_ols_log_rmsd"] = {"n": len(okp), "coef_close": float(beta[1]), "coef_plddt": float(beta[2]),
                                      "ci95_close": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                                      "n_boot_used": len(bs)}

    # A
    test = [l.split(",") for l in (root / "data/e428/alphaflow/atlas_test.csv").read_text().splitlines()[1:]]
    train = [l.split(",") for l in (root / "data/e428/alphaflow/atlas_train.csv").read_text().splitlines()[1:]]
    train_pdb = {n.split("_")[0].upper() for n, *_ in train}
    train_acc = {}
    for n, *_ in train:
        pdb, ch = n.split("_")
        for a in sifts_accessions(f, pdb, ch) or []:
            train_acc.setdefault(a, set()).add(n)
    a_rows = []
    for name, seq, rel, *_ in test:
        pdb, ch = name.split("_")
        accs = sifts_accessions(f, pdb, ch)
        rec = {"name": name, "release": rel, "seqlen": len(seq), "before_esmfold_cutoff": rel <= ESMFOLD_CUTOFF}
        rec.update(exposure(f, accs or [], seq))
        if accs is None:
            rec["sifts_failed"] = True
        h95 = search_ids(f, [seq_node(seq, 0.95)])
        h30 = search_ids(f, [seq_node(seq, 0.30)])
        rec["train_split_seq95"] = None if h95 is None else sorted(set(h95) & train_pdb)
        rec["train_split_seq30"] = None if h30 is None else sorted(set(h30) & train_pdb)
        rec["train_split_same_accession"] = sorted({m for a in (accs or []) for m in train_acc.get(a, ())})
        a_rows.append(rec)
    n = len(a_rows)
    stats_a = {
        "n": n,
        "class_counts": {c: sum(r["class"] == c for r in a_rows)
                         for c in ("SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL", "QUERY_FAILED")},
        "n_close_training_era": sum(r["class"] in ("SAME_ACCESSION", "SEQ95") for r in a_rows),
        "n_no_accession": sum(bool(r.get("no_accession")) for r in a_rows),
        "train_split_same_accession": sum(bool(r["train_split_same_accession"]) for r in a_rows),
        "train_split_seq95": sum(bool(r["train_split_seq95"]) for r in a_rows),
        "train_split_seq30": sum(bool(r["train_split_seq30"]) for r in a_rows),
        "train_split_search_failed": sum(r["train_split_seq95"] is None or r["train_split_seq30"] is None
                                         for r in a_rows),
        "released_on_or_before_esmfold_cutoff": sum(r["before_esmfold_cutoff"] for r in a_rows),
        "n_train_chains": len(train), "n_train_chains_with_sifts_accession":
            sum(1 for nm, *_ in train if any(nm in v for v in train_acc.values())),
    }

    (out_dir / "t_structures.json").write_text(json.dumps(t, indent=1, sort_keys=True))
    (out_dir / "a_chains.json").write_text(json.dumps(a_rows, indent=1, sort_keys=True))
    (out_dir / "stats.json").write_text(json.dumps({"registration": "docs/e428_registration.md",
                                                    "utc": datetime.now(timezone.utc).isoformat(),
                                                    "T": stats_t, "A": stats_a}, indent=1, sort_keys=True))
    (out_dir / "receipts.json").write_text(json.dumps(f.receipts, indent=0))
    print(json.dumps({"T": stats_t, "A": stats_a}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
