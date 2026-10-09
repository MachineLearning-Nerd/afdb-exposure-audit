"""Phase B (c) cross-checks (P5 pooled marginal only; no binned output):
 1. same implementation on the ORIGINAL frozen ledger results/e422/batches
    (should reproduce the anchored original-ledger P5 0.8580/0.9152 if the
    implementation is right);
 2. corrected ledger with an OPEN interval (lo < y < hi) to show how much the
    registered closed-interval convention matters at the 0.95 margin."""
import hashlib, json, math, os, sys
import numpy as np
ROOT = sys.argv[1]
def load(d):
    rows = []
    for k in range(1, 5):
        ck = json.load(open(os.path.join(ROOT, d, f"e422_batch{k:02d}_checkpoint.json")))
        rows.extend(ck["rows"])
    by = {}
    for r in rows:
        v = [r.get(x) for x in ("lddt", "plddt", "d_kabsch", "b_factor_z")]
        if all(x is not None and math.isfinite(float(x)) for x in v) and 0 <= float(v[1]) < 100.01:
            by.setdefault(r["accession"], []).append(float(v[0]))
    acc = sorted(by)
    return len(rows), acc, [np.array(by[a]) for a in acc]
def run(ys, closed=True):
    n = len(ys); res = {}
    for lv in (0.90, 0.95):
        covs = []
        for s in range(200):
            perm = np.random.default_rng(s).permutation(n); t = n // 2; c = (n - t) // 2
            g = lambda idx: np.concatenate([ys[i] for i in idx])
            tr, ca, te = g(perm[:t]), g(perm[t:t + c]), g(perm[t + c:])
            ctr = float(np.median(tr)); sc = np.sort(np.abs(ca - ctr))
            q = float(sc[math.ceil((sc.size + 1) * (1.0 - (1.0 - lv))) - 1])
            lo, hi = ctr - q, ctr + q
            covs.append(float(((lo <= te) & (te <= hi)).mean() if closed else ((lo < te) & (te < hi)).mean()))
        v = np.array(covs); m = float(v.mean()); se = float(v.std(ddof=1) / math.sqrt(v.size)); thr = lv - 3 * se
        res[str(lv)] = {"mean": m, "mcse": se, "threshold": thr, "outcome": "HOLD" if m >= thr else "FAIL"}
    return res
out = {}
nr, acc, ys = load("results/e422/batches")
out["original_ledger_closed"] = {"n_rows": nr, "n_accessions": len(acc), **run(ys)}
nr, acc, ys = load("results/e427/batches")
out["corrected_ledger_open_interval"] = {"n_rows": nr, "n_accessions": len(acc), **run(ys, closed=False)}
json.dump(out, open(os.path.join(ROOT, "results/e427/audit/phaseB_subagent/sensitivity_p5.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
