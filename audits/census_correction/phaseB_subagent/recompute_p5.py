"""e427 Phase B (c): independent recomputation of P5-90 / P5-95 (pooled marginal
split-conformal coverage) on the corrected e422 ledger results/e427/batches.

Written from docs/e422_registration.md section 3 / 3A.1 / section 7 (P5).
Imports NO experiments/ or moluq module.  Details the registration leaves
implicit were pinned by reading (not importing) experiments/e422_protocol.py:
  * interval coverage is closed on both ends: lo <= y <= hi with
    lo = center - q, hi = center + q (float64 ops in that order);
  * alpha_c = 1.0 - level as a float; k = ceil((n_cal+1)*(1.0-alpha_c));
  * MCSE = sd(ddof=1)/sqrt(R); threshold = level - 3.0*MCSE; HOLD iff mean >= threshold.
Computes ONLY the pooled marginal coverage per seed.  No pLDDT-binned output.

Usage: python -I recompute_p5.py <repo_root> [n_seeds]
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np

ROOT = sys.argv[1]
N_SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 200
OUT_DIR = os.path.join(ROOT, "results/e427/audit/phaseB_subagent")
LEVELS = (0.90, 0.95)

t_start = time.time()

# ---- load the corrected cumulative ledger (batches 1..4), own digest check --
rows = []
for k in range(1, 5):
    path = os.path.join(ROOT, "results/e427/batches", f"e422_batch{k:02d}_checkpoint.json")
    with open(path, "rb") as fh:
        ck = json.loads(fh.read().decode("utf-8"))
    body = {key: val for key, val in ck.items() if key != "sha256"}
    import hashlib
    dig = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert dig == ck["sha256"], f"digest mismatch batch {k}"
    rows.extend(ck["rows"])
n_rows = len(rows)

# ---- eligibility (registration section 3): all four channels finite, pLDDT in [0, 100.01)
def fin(v):
    return v is not None and math.isfinite(float(v))

by_acc = {}
n_excl = 0
for r in rows:
    y, p, d, z = r.get("lddt"), r.get("plddt"), r.get("d_kabsch"), r.get("b_factor_z")
    if not (fin(y) and fin(p) and 0.0 <= float(p) < 100.01 and fin(d) and fin(z)):
        n_excl += 1
        continue
    by_acc.setdefault(r["accession"], []).append(float(y))
accessions = sorted(by_acc)          # canonical pool order: Python sorted()
n = len(accessions)
ys = [np.asarray(by_acc[a], dtype=np.float64) for a in accessions]
counts = np.array([v.size for v in ys])
acc_digest = hashlib.sha256(json.dumps(accessions, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

# flat array + accession offsets for fast gathering
flat = np.concatenate(ys)
offs = np.concatenate([[0], np.cumsum(counts)])

def gather(idx):
    return np.concatenate([flat[offs[i]:offs[i + 1]] for i in idx])

def split(nacc, seed):
    perm = np.random.default_rng(seed).permutation(nacc)   # the only draw
    t = nacc // 2
    c = (nacc - t) // 2
    return perm[:t], perm[t:t + c], perm[t + c:]

per_seed = {str(lv): [] for lv in LEVELS}
kdiag = {str(lv): {"float_vs_exact_k_differs": 0} for lv in LEVELS}
for seed in range(N_SEEDS):
    tr, ca, te = split(n, seed)
    tr_y, ca_y, te_y = gather(tr), gather(ca), gather(te)
    center = float(np.median(tr_y))
    scores = np.sort(np.abs(ca_y - center))
    ncal = scores.size
    for lv in LEVELS:
        alpha_c = 1.0 - lv
        kf = math.ceil((ncal + 1) * (1.0 - alpha_c))
        ke = math.ceil((ncal + 1) * (1 - (1 - Fraction(str(lv)))))  # exact-rational k
        if kf != ke:
            kdiag[str(lv)]["float_vs_exact_k_differs"] += 1
        assert kf <= ncal
        q = float(scores[kf - 1])
        lo, hi = center - q, center + q
        cov = float(((lo <= te_y) & (te_y <= hi)).mean())
        per_seed[str(lv)].append({"seed": seed, "cov": cov, "center": center, "q": q,
                                  "k": kf, "n_cal": int(ncal), "n_test": int(te_y.size),
                                  "n_train": int(tr_y.size),
                                  "split_acc": [int(tr.size), int(ca.size), int(te.size)]})

summary = {}
for lv in LEVELS:
    v = np.array([e["cov"] for e in per_seed[str(lv)]])
    R = v.size
    mean = float(np.mean(v))
    sd = float(v.std(ddof=1)) if R > 1 else None
    mcse = sd / math.sqrt(R) if sd is not None else None
    thr = lv - 3.0 * mcse if mcse is not None else None
    summary[str(lv)] = {"R": R, "mean": mean, "sd": sd, "mcse": mcse, "threshold": thr,
                        "outcome": ("HOLD" if mean >= thr else "FAIL") if thr is not None else None,
                        "margin_mean_minus_threshold": (mean - thr) if thr is not None else None,
                        "margin_in_mcse": ((mean - thr) / mcse) if mcse else None,
                        "min_cov": float(v.min()), "max_cov": float(v.max())}

out = {"n_rows_ledger": n_rows, "n_rows_excluded": n_excl, "n_accessions": n,
       "canonical_accessions_sha256": acc_digest, "n_seeds": N_SEEDS,
       "k_convention_diagnostic": kdiag, "summary": summary, "per_seed": per_seed,
       "runtime_seconds": time.time() - t_start}
name = "recompute_p5.json" if N_SEEDS == 200 else f"recompute_p5_{N_SEEDS}seeds.json"
with open(os.path.join(OUT_DIR, name), "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in out.items() if k != "per_seed"}, indent=1, sort_keys=True))
