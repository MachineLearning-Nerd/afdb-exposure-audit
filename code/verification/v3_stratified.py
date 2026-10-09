"""V3: exposure-stratified re-analysis (exploratory, descriptive).

See REGISTRATION.md (V3). Independent re-implementation of the e422 pooled
split-conformal interval from the registration text (no import of
experiments/e422_protocol.py). Computes no P6/P7 quantity.
Usage: python -I v3_stratified.py <repo_root> <v2_json> <out_json>
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

SEEDS = range(200)
LEVELS = (0.90, 0.95)
BINS = [(0, 50), (50, 70), (70, 80), (80, 90), (90, 100.01)]
HIGH = {"SAME_ACCESSION", "SEQ95"}


def stratum(cls):
    if cls in HIGH:
        return "prior_relative"      # same protein or >=95% identity before cutoff
    if cls in ("SEQ30", "NOVEL"):
        return "no_close_relative"
    return "unresolved"


def e421_part(root, expo):
    rows = json.loads((root / "results/e421/labeled_residues_with_plddt.json").read_text())
    out = {}
    for name in ("all", "prior_relative", "no_close_relative", "unresolved"):
        sel = [r for r in rows if r.get("lddt") is not None and r.get("plddt") is not None
               and math.isfinite(r["lddt"]) and math.isfinite(r["plddt"]) and (name == "all" or stratum(expo.get(r["accession"], "QUERY_FAILED")) == name)]
        if not sel:
            out[name] = {"n": 0}
            continue
        p = np.array([r["plddt"] for r in sel], float)
        y = np.array([r["lddt"] for r in sel], float)
        bins = []
        for lo, hi in BINS:
            m = (p >= lo) & (p < hi)
            bins.append({"bin": [lo, hi], "n": int(m.sum()),
                         "err_rate_lddt_lt_0.60": float((y[m] < 0.60).mean()) if m.any() else None})
        rho = spearmanr(p, y).statistic if len(sel) > 2 else None
        out[name] = {"n": len(sel), "accessions": len({r["accession"] for r in sel}),
                     "spearman_plddt_lddt": float(rho) if rho is not None else None, "bins": bins}
    return out


def eligible(r):
    vals = [r.get("lddt"), r.get("plddt"), r.get("d_kabsch"), r.get("b_factor_z")]
    if any(v is None for v in vals) or not all(math.isfinite(v) for v in vals):
        return False
    return 0.0 <= r["plddt"] < 100.01


def order_stat(s, alpha):
    s = np.sort(s)
    k = math.ceil((s.size + 1) * (1 - alpha))
    assert k <= s.size
    return s[k - 1]


def e422_part(root, expo, exclude=()):
    rows = []
    for b in (1, 2, 3, 4):
        d = json.loads((root / f"results/e422/batches/e422_batch0{b}_checkpoint.json").read_text())
        rows += d["rows"]
    rows = [r for r in rows if eligible(r) and r["accession"] not in exclude]
    accs = sorted({r["accession"] for r in rows})
    a_index = {a: i for i, a in enumerate(accs)}
    ai = np.array([a_index[r["accession"]] for r in rows])
    y = np.array([r["lddt"] for r in rows], float)
    strat = np.array([stratum(expo.get(a, "QUERY_FAILED")) for a in accs])
    n = len(accs)
    res = {"n_rows": len(rows), "n_accessions": n,
           "accessions_by_stratum": {s: int((strat == s).sum()) for s in set(strat)},
           "rows_by_stratum": {s: int((strat[ai] == s).sum()) for s in set(strat)}, "levels": {}}
    for level in LEVELS:
        alpha = 1 - level
        cov_all, cov_prot, by = [], [], {s: [] for s in set(strat)}
        by_prot = {s: [] for s in set(strat)}  # NOT pre-specified (addendum A1)
        for seed in SEEDS:
            perm = np.random.default_rng(seed).permutation(n)
            t = n // 2
            c = (n - t) // 2
            role = np.empty(n, int)
            role[perm[:t]], role[perm[t:t + c]], role[perm[t + c:]] = 0, 1, 2
            r_row = role[ai]
            centre = np.median(y[r_row == 0])
            q = order_stat(np.abs(y[r_row == 1] - centre), alpha)
            te = r_row == 2
            covered = np.abs(y[te] - centre) <= q
            cov_all.append(covered.mean())
            te_acc = ai[te]
            per_prot = [covered[te_acc == a].mean() for a in np.unique(te_acc)]
            cov_prot.append(float(np.mean(per_prot)))
            for s in by:
                m = strat[te_acc] == s
                if m.any():
                    by[s].append(covered[m].mean())
                    by_prot[s].append(float(np.mean([covered[te_acc == a].mean()
                                                     for a in np.unique(te_acc[m])])))

        def agg(v):
            v = np.asarray(v, float)
            return {"mean": float(v.mean()), "mcse": float(v.std(ddof=1) / math.sqrt(v.size)), "R": int(v.size)}
        res["levels"][str(level)] = {"residue_weighted": agg(cov_all), "protein_averaged": agg(cov_prot),
                                     "by_stratum_residue_weighted": {s: agg(v) for s, v in by.items() if len(v) > 1},
                                     "by_stratum_protein_averaged_NOT_PRESPECIFIED":
                                         {s: agg(v) for s, v in by_prot.items() if len(v) > 1}}
    # per-protein lDDT summary by stratum (descriptive)
    sizes = np.bincount(ai, minlength=n)
    res["rows_per_protein_by_stratum"] = {s: {"median": float(np.median(sizes[strat == s])),
                                              "mean": float(sizes[strat == s].mean())} for s in set(strat)}
    res["median_lddt_by_stratum"] = {s: float(np.median(y[strat[ai] == s])) for s in set(strat)}
    res["frac_lddt_lt_0.60_by_stratum"] = {s: float((y[strat[ai] == s] < 0.60).mean()) for s in set(strat)}
    return res


def main(root, v2, out):
    root = Path(root)
    v2d = json.loads(Path(v2).read_text())
    expo = {a: rec["training"]["class"] for a, rec in v2d["accessions"].items()}
    res = {"schema": "paper-verification-v3", "exposure_channel": "training (<= 2018-04-30)",
           "e421": e421_part(root, expo), "e422_pooled": e422_part(root, expo),
           # post hoc (REGISTRATION.md addendum A2)
           "sensitivity_excluding_K7PQ54": e422_part(root, expo, exclude=("K7PQ54",))}
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps(res["e422_pooled"], indent=1))
    for k, v in res["e421"].items():
        print("e421", k, {kk: vv for kk, vv in v.items() if kk != "bins"})
        for b in v.get("bins", []):
            print("    ", b)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
