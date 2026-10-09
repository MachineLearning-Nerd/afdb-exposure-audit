"""V2c: homology coverage extension (REGISTRATION.md Addendum A3).

Runs the UNCHANGED V2 query and class rule (imported from v2_prior_exposure:
Fetcher, uniprot_seq, exposure; same TRAIN_CUTOFF; same model-created-date
lookup order and SAME_ACCESSION shortcut) for accessions in the corrected rosters
that V2 did not cover:
  (a) e427 corrected e422 ledger (results/e427/batches, batches 1-4; accession field only);
  (b) e427 corrected e421 relabel (results/e427/e421_relabel/census_results.json, status ok).
Cache: cache/v2c/. Usage: python -I v2c_coverage_extension.py <repo_root> <out_json>
"""
import collections
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import v2_prior_exposure as v2  # noqa: E402  (no import-time side effects)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(root, out):
    root = Path(root)
    inputs = {}

    def load(rel):
        p = root / rel
        inputs[rel] = sha(p)
        return json.loads(p.read_text())

    v2out = load("paper/arxiv/verification/out/v2_prior_exposure.json")
    covered = set(v2out["accessions"])
    e427_acc = set()
    for b in (1, 2, 3, 4):
        e427_acc |= {r["accession"] for r in load(f"results/e427/batches/e422_batch0{b}_checkpoint.json")["rows"]}
    relabel = load("results/e427/e421_relabel/census_results.json")
    rel_ok = {v["accession"] for v in relabel.values() if v["status"] == "ok"}
    new_a = sorted(e427_acc - covered)
    new_b = sorted(rel_ok - covered)

    # model-created-date sources exactly as V2 (e422 manifest first, then e421 qualifying entries)
    qual = load("results/e421/qualifying_entries.json")
    e421_acc_mc = {}
    for q in qual:
        e421_acc_mc.setdefault(q["afdb_accession"], q["model_created_date"][:10])
    man = load("results/e422/manifest_checkpoint.json")["records"]
    e422_mc = {r["accession"]: r["model_created_date"][:10] for r in man}

    f = v2.Fetcher(HERE / "cache/v2c")
    TC = v2.TRAIN_CUTOFF

    def one(acc):  # verbatim V2 per-accession logic
        seq = v2.uniprot_seq(f, acc)
        mc = e422_mc.get(acc) or e421_acc_mc.get(acc)
        rec = {"in_e427_e422_ledger": acc in e427_acc, "in_e427_e421_relabel_ok": acc in rel_ok,
               "seq_len": len(seq) if seq else None, "model_created_date": mc,
               "training": v2.exposure(f, acc, seq, TC)}
        if rec["training"]["class"] == "SAME_ACCESSION" and mc and mc >= TC:
            rec["template_bound"] = {"date": mc, "class": "SAME_ACCESSION", "implied_by_training": True}
        else:
            rec["template_bound"] = v2.exposure(f, acc, seq, mc) if mc else None
        return rec

    tables = {}
    for name, lst in (("a_e427_e422_ledger_new", new_a), ("b_e427_e421_relabel_new", new_b)):
        t = {}
        for i, acc in enumerate(lst):  # sequential (V2 used 3 threads); query/class logic unchanged
            t[acc] = one(acc)
            print(name, i + 1, len(lst), acc, t[acc]["training"]["class"],
                  (t[acc]["template_bound"] or {}).get("class"), flush=True)
        tables[name] = {"n": len(lst), "accessions": t,
                        "training_class_counts": dict(collections.Counter(r["training"]["class"] for r in t.values())),
                        "template_bound_class_counts": dict(collections.Counter(
                            (r["template_bound"] or {}).get("class") for r in t.values())),
                        "query_failed": sorted(a for a, r in t.items() if r["training"]["class"] == "QUERY_FAILED"
                                               or (r["template_bound"] or {}).get("class") == "QUERY_FAILED")}

    res = {
        "schema": "paper-verification-v2c",
        "registration": "REGISTRATION.md Addendum A3",
        "utc": datetime.now(timezone.utc).isoformat(),
        "training_cutoff": TC,
        "inputs_sha256": inputs,
        "v2_code_sha256": sha(HERE / "v2_prior_exposure.py"),
        "universe_checks": {"v2_covered": len(covered), "e427_e422_ledger": len(e427_acc),
                            "e427_e421_relabel_ok": len(rel_ok), "new_a": new_a, "new_b": new_b,
                            "expected_counts": {"a": 3, "b": 22},
                            "counts_match": len(new_a) == 3 and len(new_b) == 22,
                            "v2_e421_accessions_not_in_relabel_ok": sorted(
                                a for a, r in v2out["accessions"].items() if r["in_e421"] and a not in rel_ok)},
        "results": tables,
        "receipts_summary": {"n": len(f.receipts), "n_cached": sum(r["cached"] for r in f.receipts),
                             "status_counts": dict(collections.Counter(str(r["status"]) for r in f.receipts))},
        "receipts": f.receipts,
        "deviations": [],
        "notes": ["V2 functions imported unchanged (sha256 of v2_prior_exposure.py recorded). Requests issued "
                  "sequentially instead of V2's 3 worker threads (politeness only; queries/classes unchanged).",
                  "Batch rows of results/e427/batches are read for the `accession` field only.",
                  "RCSB is queried on the run date; entries released after V2's run (2026-10-07) cannot affect "
                  "<= cutoff counts except via retroactive metadata changes."],
    }
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: {kk: v[kk] for kk in ("n", "training_class_counts", "template_bound_class_counts",
                                                "query_failed")} for k, v in tables.items()}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
