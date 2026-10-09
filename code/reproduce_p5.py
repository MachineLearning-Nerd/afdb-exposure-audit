"""Standalone recomputation of the registered P5 coverage prediction (Table 4).

Recomputes pooled, marginal residue-level split-conformal coverage over the
cumulative batch 1-4 ledger, following the registered protocol
(registrations/e422_registration.md, section 3 / 3A.1), and writes nothing
but a small JSON summary. It imports no code from code/experiments.

Protocol (as registered):
  * a row is eligible iff lddt, plddt, d_kabsch and b_factor_z are present
    and finite and 0 <= plddt < 100.01 (pLDDT is used ONLY for this
    eligibility test; no pLDDT-conditional quantity is computed);
  * accessions = sorted(set(accession of eligible rows));
  * for seed in 0..199: perm = np.random.default_rng(seed).permutation(n);
    t = n // 2, c = (n - t) // 2; train/cal/test = perm[:t] / perm[t:t+c] /
    perm[t+c:] (accession level);
  * centre = median lDDT of train residues; q = ceil((n_cal+1)(1-alpha))-th
    smallest calibration score |y - centre|; a test residue is covered iff
    lo <= y <= hi with lo = centre - q, hi = centre + q (computed in this
    floating-point form, exactly as the registered evaluator does; testing
    |y - centre| <= q instead changes a few boundary residues);
  * mean = average coverage over the 200 splits; MCSE = sd(ddof=1)/sqrt(200);
    the prediction HOLDs iff mean >= nominal - 3 * MCSE.

Usage: python -I code/reproduce_p5.py [<derived_root>] [--open-interval]
  <derived_root> defaults to data/derived (the directory containing results/).
  --open-interval counts residues on the interval boundary as uncovered
  (the boundary sensitivity reported in Section 7.5).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

SEEDS = range(200)
LEVELS = (0.90, 0.95)
K_MCSE = 3.0


def finite(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


def load_ledger(batch_dir):
    rows = []
    for b in (1, 2, 3, 4):
        rows += json.loads((batch_dir / f"e422_batch0{b}_checkpoint.json").read_text())["rows"]
    acc, y = [], []
    for r in rows:
        vals = [r.get(k) for k in ("lddt", "plddt", "d_kabsch", "b_factor_z")]
        if all(finite(v) for v in vals) and 0.0 <= r["plddt"] < 100.01:
            acc.append(r["accession"])
            y.append(float(r["lddt"]))
    return len(rows), acc, np.asarray(y, dtype=float)


def p5(acc, y, open_interval=False):
    accessions = sorted(set(acc))
    idx = {a: i for i, a in enumerate(accessions)}
    a = np.fromiter((idx[x] for x in acc), dtype=np.int64, count=len(acc))
    n = len(accessions)
    covs = {lv: [] for lv in LEVELS}
    for seed in SEEDS:
        perm = np.random.default_rng(seed).permutation(n)
        t = n // 2
        c = (n - t) // 2
        role = np.empty(n, dtype=np.int8)
        role[perm[:t]], role[perm[t:t + c]], role[perm[t + c:]] = 0, 1, 2
        r = role[a]
        centre = float(np.median(y[r == 0]))
        scores = np.sort(np.abs(y[r == 1] - centre))
        yt = y[r == 2]
        for lv in LEVELS:
            k = math.ceil((scores.size + 1) * (1.0 - (1.0 - lv)))  # registered float form, alpha_c = 1 - level
            assert k <= scores.size
            q = float(scores[k - 1])
            lo, hi = centre - q, centre + q  # registered form: lo <= y <= hi (float order matters)
            covered = ((lo < yt) & (yt < hi)) if open_interval else ((lo <= yt) & (yt <= hi))
            covs[lv].append(float(covered.mean()))
    out = {}
    for lv in LEVELS:
        v = np.asarray(covs[lv])
        mean = float(v.mean())
        mcse = float(v.std(ddof=1) / math.sqrt(v.size))
        thr = lv - K_MCSE * mcse
        out[f"P5-{int(round(lv * 100))}"] = {
            "nominal": lv, "mean_cov": mean, "mcse": mcse, "threshold": thr,
            "outcome": "HOLD" if mean >= thr else "FAIL", "r_eff": int(v.size)}
    return n, out


def main(argv):
    open_interval = "--open-interval" in argv
    args = [x for x in argv if not x.startswith("--")]
    root = Path(args[0]) if args else Path("data/derived")
    res = {"open_interval": open_interval}
    for label, sub in (("original_e422", "results/e422/batches"), ("corrected_e427", "results/e427/batches")):
        n_rows, acc, y = load_ledger(root / sub)
        n_acc, out = p5(acc, y, open_interval)
        res[label] = {"ledger_rows": n_rows, "eligible_rows": int(y.size), "accessions": n_acc, **out}
    print(json.dumps(res, indent=1))
    # Compare with the registered evaluator's stored outcomes when available.
    ref = root / "p5_batch04_outcomes.json"
    if ref.exists() and not open_interval:
        stored = json.loads(ref.read_text())
        ok = True
        for label in ("original_e422", "corrected_e427"):
            for key, pred in (("p5_90", "P5-90"), ("p5_95", "P5-95")):
                s = stored[label][key]["outcome"]
                m = res[label][pred]
                same = (abs(s["mean_cov"] - m["mean_cov"]) < 1e-12 and abs(s["threshold"] - m["threshold"]) < 1e-12
                        and s["outcome"] == m["outcome"])
                ok &= same
                print(f"{label:15s} {pred}: recomputed mean {m['mean_cov']:.15f} thr {m['threshold']:.6f} "
                      f"{m['outcome']} | registered evaluator {s['mean_cov']:.15f} {s['outcome']} -> "
                      f"{'MATCH' if same else 'DIFFERS'}")
        print("ALL MATCH" if ok else "MISMATCH")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
