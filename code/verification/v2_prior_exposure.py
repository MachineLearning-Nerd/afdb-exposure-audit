"""V2/V2b: prior-PDB exposure of post-snapshot proteins (public RCSB/UniProt).

See REGISTRATION.md (V2, V2b). Every HTTP response is cached under
data/paper_verification/http/<sha256(request)>.json with a receipt.
Usage: python -I v2_prior_exposure.py <repo_root> <out_json>
"""
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

TRAIN_CUTOFF = "2018-04-30"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
GRAPHQL = "https://data.rcsb.org/graphql"
UA = {"User-Agent": "moluq-paper-verification/1 (research; contact via repo)"}
DELAY = 0.3


class Fetcher:
    def __init__(self, cache):
        self.cache = Path(cache)
        (self.cache / "http").mkdir(parents=True, exist_ok=True)
        self.receipts = []

    def _key(self, method, url, body):
        return hashlib.sha256(f"{method}\n{url}\n{body or ''}".encode()).hexdigest()

    def get(self, method, url, body=None):
        k = self._key(method, url, body)
        p = self.cache / "http" / f"{k}.json"
        if p.exists():
            rec = json.loads(p.read_text())
            self.receipts.append({k2: rec[k2] for k2 in ("url", "status", "sha256", "utc")} | {"cached": True})
            return rec["status"], rec["body"]
        status, text = None, None
        for attempt in range(3):
            try:
                data = body.encode() if body else None
                hdr = dict(UA)
                if body:
                    hdr["Content-Type"] = "application/json"
                req = urllib.request.Request(url, data=data, headers=hdr, method=method)
                with urllib.request.urlopen(req, timeout=60) as r:
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
        if status == 200 or status == 204:
            p.write_text(json.dumps(rec))
        self.receipts.append({k2: rec[k2] for k2 in ("url", "status", "sha256", "utc")} | {"cached": False})
        return status, text


def search_count(f, nodes):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "entry", "request_options": {"return_all_hits": True,
                                                       "results_verbosity": "compact"}}
    st, txt = f.get("POST", SEARCH, json.dumps(q, sort_keys=True))
    if st == 204:
        return 0, []
    if st != 200:
        return None, None
    ids = json.loads(txt).get("result_set", [])
    return len(ids), sorted(ids)


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


def exposure(f, acc, seq, date):
    out = {"date": date}
    n, ids = search_count(f, [acc_node(acc), date_node(date)])
    out["same_accession"] = n
    out["same_accession_ids"] = (ids or [])[:20]
    for lab, ident in (("seq95", 0.95), ("seq30", 0.30)):
        if seq is None:
            out[lab] = None
            continue
        n, ids = search_count(f, [seq_node(seq, ident), date_node(date)])
        out[lab] = n
        out[lab + "_ids"] = (ids or [])[:20]
    vals = [out["same_accession"], out["seq95"], out["seq30"]]
    if any(v is None for v in vals):
        cls = "QUERY_FAILED"
    elif out["same_accession"] > 0:
        cls = "SAME_ACCESSION"
    elif out["seq95"] > 0:
        cls = "SEQ95"
    elif out["seq30"] > 0:
        cls = "SEQ30"
    else:
        cls = "NOVEL"
    out["class"] = cls
    return out


def uniprot_seq(f, acc):
    st, txt = f.get("GET", UNIPROT.format(acc=acc))
    if st != 200 or not txt.startswith(">"):
        return None
    return "".join(l.strip() for l in txt.splitlines()[1:])


def release_dates(f, ids):
    out = {}
    ids = sorted(ids)
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        q = '{ entries(entry_ids: %s) { rcsb_id rcsb_accession_info { initial_release_date } } }' % json.dumps(
            [c.upper() for c in chunk])
        st, txt = f.get("POST", GRAPHQL, json.dumps({"query": q}))
        if st != 200:
            for c in chunk:
                out[c] = None
            continue
        for e in json.loads(txt)["data"]["entries"] or []:
            out[e["rcsb_id"].lower()] = e["rcsb_accession_info"]["initial_release_date"][:10]
        for c in chunk:
            out.setdefault(c, None)
    return out


def main(root, out):
    root = Path(root)
    f = Fetcher(root / "data/paper_verification")
    # universes
    census = json.loads((root / "results/e421/census_results.json").read_text())
    qual = json.loads((root / "results/e421/qualifying_entries.json").read_text())
    q_mc = {q["pdb_id"]: q["model_created_date"][:10] for q in qual}
    e421_ok = sorted({v["accession"] for v in census.values() if v["status"] == "ok"})
    e421_acc_mc = {}
    for q in qual:
        e421_acc_mc.setdefault(q["afdb_accession"], q["model_created_date"][:10])
    man = json.loads((root / "results/e422/manifest_checkpoint.json").read_text())["records"]
    e422_mc = {r["accession"]: r["model_created_date"][:10] for r in man}
    e422_acc = set()
    for b in (2, 3, 4):
        d = json.loads((root / f"results/e422/batches/e422_batch0{b}_checkpoint.json").read_text())
        e422_acc |= {r["accession"] for r in d["rows"]}
    e422_acc = sorted(e422_acc)

    res = {"schema": "paper-verification-v2", "training_cutoff": TRAIN_CUTOFF, "accessions": {}}
    universe = sorted(set(e421_ok) | set(e422_acc))
    def one(acc):
        seq = uniprot_seq(f, acc)
        mc = e422_mc.get(acc) or e421_acc_mc.get(acc)
        rec = {"in_e421": acc in e421_ok, "in_e422": acc in e422_acc, "seq_len": len(seq) if seq else None,
               "model_created_date": mc, "training": exposure(f, acc, seq, TRAIN_CUTOFF)}
        # a later date can only add hits: SAME_ACCESSION at the training cutoff
        # implies SAME_ACCESSION at any later template-bound date
        if rec["training"]["class"] == "SAME_ACCESSION" and mc and mc >= TRAIN_CUTOFF:
            rec["template_bound"] = {"date": mc, "class": "SAME_ACCESSION", "implied_by_training": True}
        else:
            rec["template_bound"] = exposure(f, acc, seq, mc) if mc else None
        return acc, rec

    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=3) as ex:
        for i, (acc, rec) in enumerate(ex.map(one, universe)):
            res["accessions"][acc] = rec
            print(i + 1, len(universe), acc, rec["training"]["class"],
                  rec["template_bound"]["class"] if rec["template_bound"] else None, flush=True)

    # V2b: e421 release dates vs model dates
    rel = release_dates(f, list(q_mc))
    rows = []
    for pid, mc in sorted(q_mc.items()):
        c = census.get(pid, {})
        n_res = len(c.get("labels", [])) if c.get("status") == "ok" else 0
        rows.append({"entry": pid, "model_created": mc, "release": rel.get(pid), "status": c.get("status"),
                     "n_residue_rows": n_res,
                     "release_after_model": (rel.get(pid) > mc) if rel.get(pid) else None})
    res["v2b_e421_entries"] = rows
    res["receipts"] = f.receipts
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
