"""V4 (post hoc, REGISTRATION.md addendum A2): K7PQ54 miss concentration,
per-chain label consistency, and calibration-quantile mechanism on the e422
batch-1..4 ledger, using the registered pooled interval at level 0.90.

Usage: python -I v4_k7pq54_diagnostics.py <repo_root> <out_json>
"""
import collections
import json
import math
import sys
from pathlib import Path

import numpy as np

ACC = "K7PQ54"


def eligible(r):
    vals = [r.get(k) for k in ("lddt", "plddt", "d_kabsch", "b_factor_z")]
    return all(v is not None and math.isfinite(v) for v in vals) and 0.0 <= r["plddt"] < 100.01


def main(root, out):
    root = Path(root)
    rows = []
    for b in (1, 2, 3, 4):
        rows += json.loads((root / f"results/e422/batches/e422_batch0{b}_checkpoint.json").read_text())["rows"]
    rows = [r for r in rows if eligible(r)]
    accs = sorted({r["accession"] for r in rows})
    ai = {a: i for i, a in enumerate(accs)}
    a = np.array([ai[r["accession"]] for r in rows])
    y = np.array([r["lddt"] for r in rows])
    n, k = len(accs), ai[ACC]

    miss = np.zeros(n)
    tested = np.zeros(n)
    by_role = {0: [], 1: [], 2: []}
    for seed in range(200):
        perm = np.random.default_rng(seed).permutation(n)
        t = n // 2
        c = (n - t) // 2
        role = np.empty(n, int)
        role[perm[:t]], role[perm[t:t + c]], role[perm[t + c:]] = 0, 1, 2
        rr = role[a]
        ctr = np.median(y[rr == 0])
        s = np.sort(np.abs(y[rr == 1] - ctr))
        q = s[math.ceil((s.size + 1) * 0.9) - 1]
        te = rr == 2
        m = np.abs(y[te] - ctr) > q
        np.add.at(miss, a[te], m)
        np.add.at(tested, a[te], 1)
        oth = te & (a != k)
        cov = np.abs(y[oth] - ctr) <= q
        pa = float(np.mean([cov[a[oth] == p].mean() for p in np.unique(a[oth])]))
        by_role[int(role[k])].append((float(q), pa))
    order = np.argsort(-miss)
    share = np.cumsum(miss[order]) / miss.sum()

    kr = [r for r in rows if r["accession"] == ACC]
    chains = collections.defaultdict(list)
    for r in kr:
        chains[(r["entry"], r["chain"])].append(r)
    chain_tab = [{"entry": e, "chain": ch, "n": len(v),
                  "pos_min": min(r["uniprot_pos"] for r in v), "pos_max": max(r["uniprot_pos"] for r in v),
                  "median_lddt": float(np.median([r["lddt"] for r in v])),
                  "median_d": float(np.median([r["d_kabsch"] for r in v])),
                  "median_plddt": float(np.median([r["plddt"] for r in v]))}
                 for (e, ch), v in sorted(chains.items())]
    res = {
        "schema": "paper-verification-v4", "level": 0.90, "n_proteins": n, "n_rows": len(rows),
        "k7pq54_rows": len(kr), "k7pq54_entries": len({e for e, _ in chains}), "k7pq54_chains": len(chains),
        "miss_share_top": {str(t): float(share[t - 1]) for t in (1, 3, 5, 10)},
        "top10": [{"accession": accs[i], "rows": int((a == i).sum()), "miss_rate_when_tested": float(miss[i] / tested[i])}
                  for i in order[:10]],
        "k7pq54_role_mechanism": {name: {"splits": len(v), "mean_q": float(np.mean([x[0] for x in v])),
                                         "protein_avg_cov_other_test": float(np.mean([x[1] for x in v]))}
                                  for name, v in (("train", by_role[0]), ("calibration", by_role[1]),
                                                  ("test", by_role[2]))},
        "k7pq54_chains_table": chain_tab,
    }
    Path(out).write_text(json.dumps(res, indent=1, sort_keys=True))
    print(json.dumps({k2: v for k2, v in res.items() if k2 != "k7pq54_chains_table"}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
