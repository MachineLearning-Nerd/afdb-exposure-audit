"""Standalone recomputation of the registered P5 coverage prediction (Table 4).

Recomputes pooled, marginal residue-level split-conformal coverage over the
cumulative batch 1-4 ledger, following the registered protocol
(registrations/e422_registration.md, section 3 / 3A.1). It imports no code
from code/experiments and writes nothing; it prints a JSON summary and a
comparison with the registered evaluator's stored outcomes.

Input integrity, checked before any computation (any failure -> exit 2):
  * every ledger checkpoint's stored canonical digest (`sha256` field) equals
    SHA-256 of the checkpoint minus that field, serialised as
    json.dumps(sort_keys=True, separators=(',', ':'), ensure_ascii=True) and
    UTF-8 encoded (the rule documented with the e422 runner; re-implemented
    here), and its `n_rows` equals the number of rows;
  * every ledger checkpoint's file SHA-256 equals its entry in the release's
    MANIFEST.sha256 (pass --no-manifest to skip, e.g. on a modified copy).

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
    |y - centre| <= q instead changes a few boundary residues of the
    original ledger and none of the corrected one);
  * mean = average coverage over the 200 splits; MCSE = sd(ddof=1)/sqrt(200);
    the prediction HOLDs iff mean >= nominal - 3 * MCSE.

Comparison: mean, MCSE and threshold must equal the registered evaluator's
stored values (data/derived/p5_batch04_outcomes.json) to 1e-12 and the
outcome must be identical; otherwise the exit status is 1.

Usage: python -I code/reproduce_p5.py [<derived_root>] [--open-interval] [--manifest PATH | --no-manifest]
  <derived_root> defaults to data/derived (the directory containing results/).
  --open-interval counts residues on the interval boundary as uncovered
  (the boundary sensitivity reported in Section 7.5); no comparison is made.
Exit status: 0 all inputs verified and all values match; 1 a recomputed
value differs from the registered evaluator; 2 an input is missing or fails
its digest check.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

SEEDS = range(200)
LEVELS = (0.90, 0.95)
K_MCSE = 3.0
TOL = 1e-12
LEDGERS = (("original_e422", "results/e422/batches"), ("corrected_e427", "results/e427/batches"))
REPO = Path(__file__).resolve().parents[1]


class InputError(Exception):
    pass


def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def canonical_sha256(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode("utf-8")).hexdigest()


def read_manifest(path):
    if not path.is_file():
        raise InputError(f"manifest not found: {path} (use --manifest PATH or --no-manifest)")
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            digest, name = line.split(None, 1)
            out[name.strip().lstrip("*")] = digest
    return out


def load_ledger(root, sub, manifest=None):
    """Rows of batches 1-4 after verifying each checkpoint; returns (n_rows, accessions, lddt array)."""
    rows = []
    for b in (1, 2, 3, 4):
        p = root / sub / f"e422_batch0{b}_checkpoint.json"
        if not p.is_file():
            raise InputError(f"missing ledger checkpoint: {p}")
        raw = p.read_bytes()
        ck = json.loads(raw.decode("utf-8"))
        stored = ck.get("sha256")
        calc = canonical_sha256({k: v for k, v in ck.items() if k != "sha256"})
        if stored != calc:
            raise InputError(f"{p}: stored canonical digest {stored} != recomputed {calc}")
        if ck.get("n_rows") != len(ck.get("rows", [])):
            raise InputError(f"{p}: n_rows {ck.get('n_rows')} != {len(ck.get('rows', []))} rows")
        if manifest is not None:
            key = f"data/derived/{sub}/e422_batch0{b}_checkpoint.json"
            want = manifest.get(key)
            got = hashlib.sha256(raw).hexdigest()
            if want is None:
                raise InputError(f"{key} has no MANIFEST.sha256 entry")
            if want != got:
                raise InputError(f"{p}: file SHA-256 {got} != MANIFEST.sha256 {want}")
        rows += ck["rows"]
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
            if k > scores.size:
                raise ValueError(f"seed {seed}: order statistic {k} exceeds {scores.size} calibration scores")
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
    manifest_path = REPO / "MANIFEST.sha256"
    use_manifest = "--no-manifest" not in argv
    args = []
    it = iter(argv)
    for x in it:
        if x == "--manifest":
            manifest_path = Path(next(it, ""))
        elif not x.startswith("--"):
            args.append(x)
    root = Path(args[0]) if args else Path("data/derived")
    try:
        manifest = read_manifest(manifest_path) if use_manifest else None
        res = {"open_interval": open_interval}
        for label, sub in LEDGERS:
            n_rows, acc, y = load_ledger(root, sub, manifest)
            n_acc, out = p5(acc, y, open_interval)
            res[label] = {"ledger_rows": n_rows, "eligible_rows": int(y.size), "accessions": n_acc, **out}
        print("inputs verified: 8 checkpoint canonical digests"
              + (f" and MANIFEST.sha256 entries ({manifest_path})" if manifest is not None else
                 " (MANIFEST.sha256 check skipped)"))
        print(json.dumps(res, indent=1))
        if open_interval:
            return 0
        ref = root / "p5_batch04_outcomes.json"
        if not ref.is_file():
            raise InputError(f"registered outcomes not found: {ref}")
        stored = json.loads(ref.read_text())
    except (InputError, OSError, ValueError, KeyError) as e:
        print(f"INPUT ERROR: {type(e).__name__}: {e}")
        return 2
    ok = True
    for label, _ in LEDGERS:
        for key, pred in (("p5_90", "P5-90"), ("p5_95", "P5-95")):
            s = stored.get(label, {}).get(key, {}).get("outcome", {})
            m = res[label][pred]
            diffs = []
            for f in ("mean_cov", "mcse", "threshold"):
                if not finite(s.get(f)) or abs(s[f] - m[f]) > TOL:
                    diffs.append(f"{f} {s.get(f)} vs {m[f]!r}")
            for f in ("outcome", "r_eff"):
                if s.get(f) != m[f]:
                    diffs.append(f"{f} {s.get(f)} vs {m[f]}")
            ok &= not diffs
            print(f"{label:15s} {pred}: mean {m['mean_cov']:.15f} MCSE {m['mcse']:.15f} thr {m['threshold']:.15f} "
                  f"{m['outcome']} | registered evaluator -> {'MATCH' if not diffs else 'DIFFERS: ' + '; '.join(diffs)}")
    print("ALL MATCH" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
