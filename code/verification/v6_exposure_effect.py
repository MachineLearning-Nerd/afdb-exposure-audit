"""V6: EXPLORATORY exposure effect on the e421 corrected relabel (REGISTRATION.md Addendum A3).

Input: results/e427/e421_relabel/census_results.json (status ok entries only).
Unit = accession (residue rows pooled over that accession's ok entries).
Per accession: median lDDT; fraction lDDT < 0.60; fraction (pLDDT >= 90 and lDDT < 0.60);
median pLDDT. Groups from the V2/V2c training-cutoff class:
CLOSE = SAME_ACCESSION or SEQ95; NOT_CLOSE = SEQ30 or NOVEL (QUERY_FAILED/missing excluded, listed).
Statistic: mean(CLOSE) - mean(NOT_CLOSE); 95% percentile bootstrap over accessions,
10,000 resamples, seed 0. No e422/e427 batch ledger is read.
Usage: python -I v6_exposure_effect.py <repo_root> <out_json>
"""
import collections
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
B, SEED, ERR, CONF = 10_000, 0, 0.60, 90.0
CLOSE = {"SAME_ACCESSION", "SEQ95"}
NOT_CLOSE = {"SEQ30", "NOVEL"}
QTY = ("median_lddt", "frac_lddt_lt_060", "frac_confident_error", "median_plddt")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fin(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def main(root, out):
    root = Path(root)
    inputs = {}
    paths = {"census": root / "results/e427/e421_relabel/census_results.json",
             "v2": HERE / "out/v2_prior_exposure.json", "v2c": HERE / "out/v2c_extension.json"}
    for p in paths.values():
        inputs[str(p.relative_to(root))] = sha(p)
    census = json.loads(paths["census"].read_text())
    v2 = json.loads(paths["v2"].read_text())["accessions"]
    v2c = json.loads(paths["v2c"].read_text())["results"]

    cls, src = {}, {}
    for a, r in v2.items():
        cls[a], src[a] = r["training"]["class"], "V2"
    for tab in v2c.values():
        for a, r in tab["accessions"].items():
            if a in cls and cls[a] != r["training"]["class"]:
                raise SystemExit(f"class conflict V2 vs V2c for {a}")
            cls.setdefault(a, r["training"]["class"])
            src.setdefault(a, "V2c")

    rows = collections.defaultdict(lambda: ([], []))
    entries = collections.defaultdict(list)
    n_nonfinite = 0
    for eid, e in sorted(census.items()):
        if e["status"] != "ok":
            continue
        entries[e["accession"]].append(eid)
        for lab in e["labels"]:
            if lab["accession"] != e["accession"]:
                raise SystemExit(f"label accession mismatch in {eid}")
            if not (fin(lab["lddt"]) and fin(lab["plddt"])):
                n_nonfinite += 1
                continue
            rows[e["accession"]][0].append(lab["lddt"])
            rows[e["accession"]][1].append(lab["plddt"])

    per = {}
    for a, (l, p) in sorted(rows.items()):
        l, p = np.asarray(l, float), np.asarray(p, float)
        per[a] = {"n_entries": len(entries[a]), "n_residues": int(l.size),
                  "median_lddt": float(np.median(l)), "frac_lddt_lt_060": float(np.mean(l < ERR)),
                  "frac_confident_error": float(np.mean((p >= CONF) & (l < ERR))),
                  "median_plddt": float(np.median(p)),
                  "training_class": cls.get(a, "MISSING"), "class_source": src.get(a)}

    gC = sorted(a for a in per if per[a]["training_class"] in CLOSE)
    gN = sorted(a for a in per if per[a]["training_class"] in NOT_CLOSE)
    excluded = sorted(a for a in per if a not in gC and a not in gN)
    XC = np.array([[per[a][q] for q in QTY] for a in gC])
    XN = np.array([[per[a][q] for q in QTY] for a in gN])
    rng = np.random.default_rng(SEED)
    iC = rng.integers(0, len(gC), size=(B, len(gC)))
    iN = rng.integers(0, len(gN), size=(B, len(gN)))
    boot = XC[iC].mean(axis=1) - XN[iN].mean(axis=1)  # (B, n_qty)
    eff = {}
    for k, q in enumerate(QTY):
        lo, hi = np.percentile(boot[:, k], [2.5, 97.5])
        eff[q] = {"mean_CLOSE": float(XC[:, k].mean()), "mean_NOT_CLOSE": float(XN[:, k].mean()),
                  "diff_CLOSE_minus_NOT_CLOSE": float(XC[:, k].mean() - XN[:, k].mean()),
                  "ci95_percentile_bootstrap": [float(lo), float(hi)],
                  "median_CLOSE": float(np.median(XC[:, k])), "median_NOT_CLOSE": float(np.median(XN[:, k]))}

    four = {}
    for c in ("SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL"):
        m = [a for a in per if per[a]["training_class"] == c]
        four[c] = {"n_accessions": len(m), "n_residues": int(sum(per[a]["n_residues"] for a in m)),
                   **({q: {"mean": float(np.mean([per[a][q] for a in m])),
                           "median": float(np.median([per[a][q] for a in m]))} for q in QTY} if m else {})}

    res = {
        "schema": "paper-verification-v6-exposure-effect",
        "label": "EXPLORATORY (descriptive; no significance claim)",
        "registration": "REGISTRATION.md Addendum A3",
        "utc": datetime.now(timezone.utc).isoformat(),
        "inputs_sha256": inputs,
        "receipts_summary": {"network_requests": 0, "e422_or_e427_batch_ledgers_read": False},
        "parameters": {"bootstrap_resamples": B, "seed": SEED, "lddt_error_threshold": ERR,
                       "confident_plddt": CONF, "groups": {"CLOSE": sorted(CLOSE), "NOT_CLOSE": sorted(NOT_CLOSE)},
                       "bootstrap_scheme": "independent resampling of accessions with replacement within each group "
                                           "(group sizes fixed); one index draw shared by all four quantities"},
        "results": {
            "n_ok_entries": sum(len(v) for v in entries.values()),
            "n_accessions": len(per), "n_residues": int(sum(v["n_residues"] for v in per.values())),
            "n_nonfinite_rows_excluded": n_nonfinite,
            "group_sizes": {"CLOSE": len(gC), "NOT_CLOSE": len(gN), "excluded": len(excluded)},
            "group_residues": {"CLOSE": int(sum(per[a]["n_residues"] for a in gC)),
                               "NOT_CLOSE": int(sum(per[a]["n_residues"] for a in gN))},
            "excluded_accessions": {a: per[a]["training_class"] for a in excluded},
            "effects": eff,
            "four_class_breakdown": four,
            "class_source_counts": dict(collections.Counter(per[a]["class_source"] for a in per)),
            "per_accession": per,
        },
        "deviations": [],
        "notes": ["Registration is silent on the bootstrap resampling scheme; stratified (within-group) resampling "
                  "was chosen before computing any V6 number.",
                  "Rows with non-finite lDDT or pLDDT would be excluded (count reported); quantities are fractions "
                  "of finite rows."],
    }
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k: res["results"][k] for k in ("n_ok_entries", "n_accessions", "n_residues", "group_sizes",
                                                      "group_residues", "excluded_accessions", "effects")}, indent=1))
    print(json.dumps({c: {"n": v["n_accessions"], **{q: round(v[q]["mean"], 4) for q in QTY if q in v}}
                      for c, v in four.items()}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
