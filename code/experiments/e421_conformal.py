"""e421 AD-CP split-conformal analysis: does pLDDT give valid coverage?

Offline over the committed exploratory artifact
results/e421/labeled_residues_with_plddt.json (24,457 rows). No network,
no source access, deterministic seeds.

Question (assignment H11-14): using non-conformity score s = 1 - pLDDT/100,
does pLDDT provide valid prediction coverage at the 90% and 95% levels, and
what is the empirical-vs-nominal gap?

Two score views, both derived from s:
  - classwise split conformal (primary): confidence in class k is pLDDT for
    k=correct and 100-pLDDT for k=error, so s_correct = 1 - pLDDT/100 and
    s_error = pLDDT/100. Prediction set S(x) = {k : s_k(x) <= q_k} with q_k the
    ceil((n_k+1)(1-alpha)) order statistic of class-k calibration scores.
    Guarantee (exchangeability): P(Y in S) >= 1-alpha marginally and per class.
  - per-class recall view (assignment-literal score, reported as the per-class
    component of the same calibration): empirical P(s <= q | Y=k).

Two split protocols, because the 24,457 residues cluster within proteins:
  - residue_random: stratified 50/50 residue-level split (exchangeable within
    the pooled set, but shares proteins across calibration/test);
  - accession_level: 50/50 split by UniProt accession (protein-disjoint; the
    scientifically meaningful generalisation test for AD-CP). The accession
    split is drawn ONCE over the union of accessions and applied to all rows,
    so every accession has test probability exactly 0.5 regardless of how many
    classes it contributes rows to.

R = 200 repetitions per protocol/level; canonical split seed 0 additionally
reports per-pLDDT-bin conditional coverage (marginal validity can hide
conditional failure).
"""
import hashlib
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_PATH = os.path.join(ROOT, "results", "e421", "labeled_residues_with_plddt.json")
OUT_PATH = os.path.join(ROOT, "results", "e421", "adcp_conformal_analysis.json")

LEVELS = (0.90, 0.95)
N_REPEATS = 200
TEST_FRACTION = 0.5
CLASSES = ("correct", "error")
PLDDT_BINS = ((0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100.01))


def conformal_quantile(scores, coverage):
    """ceil((n+1)*coverage) order statistic (1-indexed); inf if it exceeds n.

    `coverage` is the target coverage 1-alpha (e.g. 0.90), so the threshold is
    the ceil((n+1)(1-alpha)) = ceil((n+1)*coverage)-th smallest score.
    """
    s = np.sort(np.asarray(scores, dtype=float))
    n = s.size
    k = math.ceil((n + 1) * coverage)
    return float(s[k - 1]) if k <= n else math.inf


def score_for(class_name, plddt):
    """Non-conformity of plddt for membership in class_name, on 0-1 scale."""
    p = np.asarray(plddt, dtype=float) / 100.0
    return (1.0 - p) if class_name == "correct" else p


def one_repeat(rng, data, level, accession_level):
    """One stratified split -> coverage metrics at one level.

    data: {class: list of (plddt, accession, bin_idx)} rows with definitive e.
    """
    cal, test = {"correct": [], "error": []}, {"correct": [], "error": []}
    if accession_level:
        # one draw over the union of accessions; rows follow their accession.
        # (A per-class draw into shared sets would double-draw dual-class
        # accessions and shrink the realized test fraction below 0.5 - the
        # v1 artifact's protocol, corrected after independent audit.)
        accs = sorted({a for c in CLASSES for _, a, _ in data[c]})
        idx = rng.permutation(len(accs))
        n_cal = max(1, int(round(len(accs) * TEST_FRACTION)))
        cal_accs = {accs[j] for j in idx[:n_cal]}
        for c in CLASSES:
            for row in data[c]:
                (cal if row[1] in cal_accs else test)[c].append(row)
    else:
        for c in CLASSES:
            idx = rng.permutation(len(data[c]))
            n_cal = int(round(len(data[c]) * TEST_FRACTION))
            for i, j in enumerate(idx):
                (cal if i < n_cal else test)[c].append(data[c][j])

    q = {}
    for c in CLASSES:
        s_cal = score_for(c, [r[0] for r in cal[c]])
        q[c] = conformal_quantile(s_cal, level)
        if not math.isfinite(q[c]):  # cannot happen at these n; fail loud, not silent
            raise RuntimeError(f"non-finite conformal quantile for class {c}")

    # per-class coverage + set composition on test
    cov = {}
    for c in CLASSES:
        s_test = score_for(c, [r[0] for r in test[c]])
        cov[c] = float(np.mean(s_test <= q[c]))
    all_test = [(c, r) for c in CLASSES for r in test[c]]
    in_set = np.zeros((len(all_test), 2), dtype=bool)
    for ki, c in enumerate(CLASSES):
        s = score_for(c, [r[0] for _, r in all_test])
        in_set[:, ki] = s <= q[c]
    y = np.array([0 if c == "correct" else 1 for c, _ in all_test])
    marginal = float(np.mean(in_set[np.arange(len(y)), y]))
    sizes = in_set.sum(axis=1)
    return {
        "marginal": marginal,
        "per_class": cov,
        "set_size_mean": float(np.mean(sizes)),
        "empty_rate": float(np.mean(sizes == 0)),
        "singleton_rate": float(np.mean(sizes == 1)),
        "doubleton_rate": float(np.mean(sizes == 2)),
    }


def aggregate(reps, nominal):
    keys = reps[0].keys()
    out = {}
    for k in keys:
        if not np.isscalar(reps[0][k]):  # per_class dict of scalars
            vals = np.array([[r[k][c] for c in ("correct", "error")] for r in reps], dtype=float)
            out[k] = {
                c: {
                    "mean": round(float(vals[:, i].mean()), 5),
                    "std": round(float(vals[:, i].std(ddof=1)), 5),
                    "mcse": round(float(vals[:, i].std(ddof=1) / math.sqrt(len(vals))), 5),
                }
                for i, c in enumerate(("correct", "error"))
            }
            continue
        vals = np.array([r[k] for r in reps], dtype=float)
        out[k] = {
                "mean": round(float(vals.mean()), 5),
                "std": round(float(vals.std(ddof=1)), 5),
                "mcse": round(float(vals.std(ddof=1) / math.sqrt(len(vals))), 5),
                "min": round(float(vals.min()), 5),
                "max": round(float(vals.max()), 5),
            }
    out["empirical_minus_nominal"] = round(out["marginal"]["mean"] - nominal, 5)
    return out


def main():
    rows = json.load(open(IN_PATH))
    assert len(rows) == 24457, f"unexpected row count {len(rows)}"

    # ---- assemble class-evaluable rows; typed accounting for the rest
    n_e_ambiguous = sum(1 for r in rows if r["e"] not in CLASSES)
    data = {c: [] for c in CLASSES}
    for r in rows:
        if r["e"] in CLASSES:
            bin_idx = next(i for i, (lo, hi) in enumerate(PLDDT_BINS) if lo <= r["plddt"] < hi)
            data[r["e"]].append((r["plddt"], r["accession"], bin_idx))
    n_by_class = {c: len(data[c]) for c in CLASSES}

    # ---- descriptive recheck (ties this analysis to adcp_analysis.json)
    plddt_all = np.array([r["plddt"] for r in rows], dtype=float)
    lddt_all = np.array([r["lddt"] if r["lddt"] is not None else np.nan for r in rows], dtype=float)
    ok = ~np.isnan(lddt_all)
    rho, pval = spearmanr(plddt_all[ok], lddt_all[ok])
    err_by_bin = {}
    for lo, hi in ((0, 50), (50, 70), (70, 80), (80, 90), (90, 100.01)):
        m = (plddt_all >= lo) & (plddt_all < hi)
        err_by_bin[f"{lo}-{hi if hi <= 100 else 100}"] = {
            "n": int(m.sum()),
            "error_rate_lddt_below_060": round(float(np.mean(lddt_all[m] < 0.60)), 4),
        }
    descriptives = {
        "spearman_plddt_vs_lddt": {
            "rho": round(float(rho), 4),
            "p_value": ("<1e-300" if pval < 1e-300 else float(pval)),
            "n_nonnull_lddt": int(ok.sum()),
            "n_null_lddt_excluded": int((~ok).sum()),
            "adcp_analysis_recorded_rho": 0.617,
            "matches_adcp_analysis": bool(abs(rho - 0.617) < 5e-4),
        },
        "error_rate_by_plddt_bin": err_by_bin,
        "error_rate_definition": "fraction of rows with lddt < 0.60 (matches "
                                 "adcp_analysis.json); strict e-channel error additionally "
                                 "requires dist > 4.0 A",
    }

    # ---- conformal protocols
    protocols = {}
    for name, acc_level in (("residue_random_split", False), ("accession_level_split", True)):
        protocols[name] = {}
        for level in LEVELS:
            reps = [one_repeat(np.random.default_rng(seed), data, level, acc_level)
                    for seed in range(N_REPEATS)]
            protocols[name][f"nominal_{int(level*100)}"] = aggregate(reps, level)

    # ---- canonical split (seed 0, residue-random): thresholds + conditional coverage
    canonical = {}
    rng0 = np.random.default_rng(0)
    for level in LEVELS:
        cal, test = {"correct": [], "error": []}, {"correct": [], "error": []}
        for c in CLASSES:
            idx = rng0.permutation(len(data[c]))
            n_cal = int(round(len(data[c]) * TEST_FRACTION))
            for i, j in enumerate(idx):
                (cal if i < n_cal else test)[c].append(data[c][j])
        entry = {"thresholds": {}, "conditional_coverage_by_plddt_bin": {}}
        for c in CLASSES:
            entry["thresholds"][c] = round(conformal_quantile(
                score_for(c, [r[0] for r in cal[c]]), level), 5)
        for c in CLASSES:
            q_c = entry["thresholds"][c]
            s_cal = score_for(c, [r[0] for r in cal[c]])
            s_test = score_for(c, [r[0] for r in test[c]])
            entry.setdefault("per_class_recall_view", {})[c] = {
                "calibration_n": int(s_cal.size),
                "test_n": int(s_test.size),
                "empirical_coverage": round(float(np.mean(s_test <= q_c)), 5),
                "nominal": level,
            }
        # conditional: per pLDDT bin of test rows, P(true class in set)
        all_test = [(c, r) for c in CLASSES for r in test[c]]
        for bi, (lo, hi) in enumerate(PLDDT_BINS):
            members = [(c, r) for c, r in all_test if r[2] == bi]
            if not members:
                continue
            hits = 0
            for c, r in members:
                s_ok = (score_for("correct", [r[0]])[0] <= entry["thresholds"]["correct"]) if c == "correct" \
                    else (score_for("error", [r[0]])[0] <= entry["thresholds"]["error"])
                hits += bool(s_ok)
            entry["conditional_coverage_by_plddt_bin"][f"{lo}-{int(hi) if hi <= 100 else 100}"] = {
                "n": len(members),
                "coverage": round(hits / len(members), 5),
                "nominal": level,
            }
        canonical[f"nominal_{int(level*100)}"] = entry

    # ---- verdicts
    verdicts, notes = {}, []
    all_valid = True
    for name in protocols:
        for lv in protocols[name]:
            m = protocols[name][lv]["marginal"]
            nominal = int(lv.split("_")[1]) / 100
            ok_v = m["mean"] >= nominal - 3 * m["mcse"]
            verdicts[f"{name}/{lv}"] = "COVERAGE_VALID" if ok_v else "COVERAGE_INVALID"
            all_valid &= ok_v
    # conditional-coverage honesty check on the canonical split
    for lv in canonical:
        for b, d in canonical[lv]["conditional_coverage_by_plddt_bin"].items():
            if d["coverage"] < d["nominal"] - 0.05 and d["n"] >= 30:
                notes.append(f"conditional coverage below nominal in pLDDT bin {b} "
                             f"({lv}): {d['coverage']} vs {d['nominal']} (n={d['n']})")
    accession_ok = all(
        protocols["accession_level_split"][lv]["marginal"]["mean"] >= int(lv.split("_")[1]) / 100 - 3 * protocols["accession_level_split"][lv]["marginal"]["mcse"]
        for lv in protocols["accession_level_split"])

    verdict = ("COVERAGE_VALID_MARGINAL" if all_valid else "COVERAGE_INVALID") + \
              ("_PROTEIN_DISJOINT_OK" if accession_ok else "_PROTEIN_DISJOINT_FAILS")

    out = {
        "schema": "e421-adcp-conformal-v2",
        "status": "EXPLORATORY_NOT_FULL_REGISTERED",
        "v2_correction": "accession_level_split redrawn as a single 50/50 draw over "
                         "the union of accessions applied to all rows (v1 drew "
                         "per-class into shared sets, giving dual-class accessions "
                         "test probability 0.25 and a realized test fraction of "
                         "~0.33; found by the detached e421 final audit, 2026-09-08). "
                         "residue_random_split and canonical_split_seed0 are "
                         "byte-identical to v1 (same seeds, same draws).",
        "exploratory_note": "Computed offline from the exploratory e421 ledger "
                            "(provenance envelopes incomplete); not admissible for the "
                            "AD-CP positive claim without a prospective re-run.",
        "question": "Does pLDDT provide valid conformal coverage for e-channel class "
                    "membership at 90%/95%, and what is the empirical-vs-nominal gap?",
        "score_definition": "s = 1 - pLDDT/100; classwise view: s_correct = 1 - pLDDT/100, "
                            "s_error = pLDDT/100 (confidence in class k is pLDDT for correct, "
                            "100-pLDDT for error)",
        "calibration_rule": "q_k = ceil((n_k+1)(1-alpha)) order statistic of class-k "
                            "calibration scores; prediction set {k : s_k(x) <= q_k}",
        "n_rows_total": len(rows),
        "n_class_evaluable": n_by_class,
        "n_e_ambiguous_excluded": n_e_ambiguous,
        "n_accessions": len({r["accession"] for r in rows}),
        "levels": list(LEVELS),
        "n_repeats": N_REPEATS,
        "test_fraction": TEST_FRACTION,
        "descriptive_recheck": descriptives,
        "protocols": protocols,
        "canonical_split_seed0": canonical,
        "verdicts_per_protocol_level": verdicts,
        "verdict": verdict,
        "findings": notes,
        "digest_rule": "analysis_digest = SHA256(UTF8(json.dumps(payload, sort_keys=True, "
                       "separators=(',',':')))) over all fields except analysis_digest",
    }
    out["analysis_digest"] = hashlib.sha256(
        json.dumps({k: v for k, v in out.items() if k != "analysis_digest"},
                   sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()

    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"classes: {n_by_class} (excluded e-ambiguous: {n_e_ambiguous})")
    print(f"spearman recheck: rho={rho:.4f} (adcp recorded 0.617)")
    for name in protocols:
        for lv in protocols[name]:
            agg = protocols[name][lv]
            m, pc = agg["marginal"], agg["per_class"]
            print(f"{name} {lv}: marginal={m['mean']:.4f} (nominal {lv.split('_')[1]}/100, "
                  f"gap {agg['empirical_minus_nominal']:+.4f}, "
                  f"mcse {m['mcse']:.4f}) set_size={agg['set_size_mean']['mean']:.3f} "
                  f"empty={agg['empty_rate']['mean']:.3f} | cov_correct={pc['correct']['mean']:.4f} "
                  f"cov_error={pc['error']['mean']:.4f}")
    for lv in canonical:
        print(f"canonical {lv} conditional:", {b: d["coverage"] for b, d in
                                               canonical[lv]["conditional_coverage_by_plddt_bin"].items()})
    print("VERDICT:", verdict)
    for n in notes:
        print("FINDING:", n)
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
