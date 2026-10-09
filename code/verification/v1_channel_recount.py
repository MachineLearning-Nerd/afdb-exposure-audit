"""V1: channel-separated recount of the e420 ATLAS temporal audit.

Reads only committed artifacts. See REGISTRATION.md (V1).
  training channel (T): every dated reference entry released <= 2018-04-30
  template channel (P): every dated reference entry released <= its AFDB
                        modelCreatedDate (the paper's original rule)
Usage: python -I v1_channel_recount.py <repo_root> <out_json>
"""
import collections
import hashlib
import json
import re
import sys
from pathlib import Path

TRAIN_CUTOFF = "2018-04-30"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(root, out):
    root = Path(root)
    rec_p = root / "results/e420/label_run_receipt.json"
    ovf_p = root / "results/e420/overflow_temporal_check.json"
    rec = json.loads(rec_p.read_text())
    ovf = json.loads(ovf_p.read_text())

    # admitted pool: one row per (accession, entry, chain); collapse to entry
    adm = {}
    undated = []
    pat = re.compile(r"modelCreatedDate=(\d{4}-\d{2}-\d{2})\S*\s*!<\s*release=(\d{4}-\d{2}-\d{2})")
    for t in rec["typed_failures"]:
        m = pat.search(t.get("detail", ""))
        if m:
            adm[(t["accession"], t["entry"])] = (m.group(1), m.group(2))
        else:
            undated.append({"accession": t["accession"], "code": t["code"]})
    ovf_rows = {}
    for r in ovf["temporal_fail"] + ovf["temporal_pass"]:
        ovf_rows[(r["accession"], r["entry"])] = (r["model_created"][:10], r["release"][:10])

    def summarise(rows):
        by_acc = collections.defaultdict(list)
        for (acc, _e), (mc, rel) in rows.items():
            by_acc[acc].append((mc, rel))
        n_ent = len(rows)
        ent_T = sum(rel <= TRAIN_CUTOFF for mc, rel in rows.values())
        ent_P = sum(rel <= mc for mc, rel in rows.values())
        acc_T = sum(all(rel <= TRAIN_CUTOFF for _, rel in v) for v in by_acc.values())
        acc_P = sum(all(rel <= mc for mc, rel in v) for v in by_acc.values())
        acc_TP = sum(all(rel <= TRAIN_CUTOFF for _, rel in v) and all(rel <= mc for mc, rel in v)
                     for v in by_acc.values())
        return {
            "entries": n_ent,
            "entries_training_channel": ent_T,
            "entries_template_channel": ent_P,
            "accessions": len(by_acc),
            "accessions_all_refs_pre_training_cutoff": acc_T,
            "accessions_all_refs_pre_model_date": acc_P,
            "accessions_both": acc_TP,
            "accessions_with_post_cutoff_ref": len(by_acc) - acc_T,
            "model_created_dates": dict(collections.Counter(mc for mc, _ in rows.values())),
        }

    union = {**adm, **ovf_rows}
    years = collections.Counter(rel[:4] for _, rel in union.values())
    res = {
        "schema": "paper-verification-v1",
        "training_cutoff": TRAIN_CUTOFF,
        "inputs": {str(rec_p.relative_to(root)): sha(rec_p), str(ovf_p.relative_to(root)): sha(ovf_p)},
        "admitted_pool": summarise(adm),
        "admitted_undated_typed_failures": undated,
        "overflow": summarise(ovf_rows),
        "union": summarise(union),
        "release_year_histogram": dict(sorted(years.items())),
    }
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    for k in ("admitted_pool", "overflow", "union"):
        s = dict(res[k]); s.pop("model_created_dates")
        print(k, s)
    print("undated:", undated)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
