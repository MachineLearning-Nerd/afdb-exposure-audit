#!/usr/bin/env python3
"""Fresh, standalone reproduction of registered e422 P5 and post-hoc checks.

Implements only the pooled interval and split algorithm written in
docs/e422_registration.md §3. It imports neither e422_protocol.py nor the
paper's verification scripts. Run from outside the repository with
``python -I``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_eligible(row: dict) -> tuple[bool, str | None]:
    for key, reason in (("lddt", "null_label"), ("plddt", "null_feature"),
                        ("d_kabsch", "null_descriptive_channel"),
                        ("b_factor_z", "null_descriptive_channel")):
        if row.get(key) is None:
            return False, reason
        try:
            value = float(row[key])
        except (TypeError, ValueError):
            return False, reason
        if not math.isfinite(value):
            return False, "nonfinite_label" if key == "lddt" else (
                "nonfinite_feature" if key == "plddt" else "null_descriptive_channel")
    p = float(row["plddt"])
    if not 0.0 <= p < 100.01:
        return False, "feature_out_of_grid"
    return True, None


def aggregate(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    return {
        "mean": float(arr.mean()),
        "mcse": float(arr.std(ddof=1) / math.sqrt(arr.size)) if arr.size > 1 else None,
        "R": int(arr.size),
    }


def classify(raw_class: str | None) -> str:
    if raw_class in ("SAME_ACCESSION", "SEQ95"):
        return "prior_relative"
    if raw_class in ("SEQ30", "NOVEL"):
        return "no_close_relative"
    return "unresolved"


def reproduce(rows: list[dict], alpha_levels: tuple[float, ...],
              strata: dict[str, str] | None = None,
              mechanism_accession: str | None = None) -> dict:
    accessions = sorted({r["accession"] for r in rows})
    index = {acc: i for i, acc in enumerate(accessions)}
    row_acc = np.asarray([index[r["accession"]] for r in rows], dtype=int)
    y = np.asarray([float(r["lddt"]) for r in rows], dtype=float)
    acc_stratum = np.asarray([
        classify((strata or {}).get(acc)) for acc in accessions
    ], dtype=object)
    n = len(accessions)
    result = {
        "n_rows": int(len(rows)),
        "n_accessions": n,
        "accessions_by_stratum": {s: int(np.sum(acc_stratum == s)) for s in
                                  ("prior_relative", "no_close_relative", "unresolved")},
        "rows_by_stratum": {
            s: int(np.sum(acc_stratum[row_acc] == s)) for s in
            ("prior_relative", "no_close_relative", "unresolved")
        },
        "levels": {},
    }
    mechanism: dict[str, list] = {
        role: [] for role in ("train", "calibration", "test")
    }

    for level in alpha_levels:
        alpha = 1.0 - level
        overall_residue, overall_protein = [], []
        per_stratum_residue = defaultdict(list)
        per_stratum_protein = defaultdict(list)
        misses_by_accession = Counter()
        test_rows_evaluated = 0
        for seed in range(200):
            permutation = np.random.default_rng(seed).permutation(n)
            n_train = n // 2
            n_cal = (n - n_train) // 2
            train = permutation[:n_train]
            cal = permutation[n_train:n_train + n_cal]
            test = permutation[n_train + n_cal:]
            role = np.full(n, 2, dtype=np.int8)
            role[train] = 0
            role[cal] = 1
            row_role = role[row_acc]
            center = float(np.median(y[row_role == 0]))
            cal_scores = np.sort(np.abs(y[row_role == 1] - center))
            k = math.ceil((len(cal_scores) + 1) * (1 - alpha))
            if not 1 <= k <= len(cal_scores):
                raise ValueError(f"registered precondition fails: seed={seed}, k={k}, ncal={len(cal_scores)}")
            half_width = float(cal_scores[k - 1])
            test_mask = row_role == 2
            test_y = y[test_mask]
            test_acc = row_acc[test_mask]
            covered = np.abs(test_y - center) <= half_width
            overall_residue.append(float(covered.mean()))
            test_rows_evaluated += int(len(covered))
            for acc_idx in np.unique(test_acc):
                acc_name = accessions[int(acc_idx)]
                misses_by_accession[acc_name] += int(np.sum(~covered[test_acc == acc_idx]))
            protein_rates = [float(covered[test_acc == a].mean()) for a in np.unique(test_acc)]
            overall_protein.append(float(np.mean(protein_rates)))

            for s in ("prior_relative", "no_close_relative", "unresolved"):
                selected = acc_stratum[test_acc] == s
                if selected.any():
                    per_stratum_residue[s].append(float(covered[selected].mean()))
                    stratum_accs = np.unique(test_acc[selected])
                    per_stratum_protein[s].append(float(np.mean([
                        covered[test_acc == acc].mean() for acc in stratum_accs
                    ])))

            if mechanism_accession is not None and level == 0.90:
                k7_idx = index.get(mechanism_accession)
                if k7_idx is not None:
                    k7_role = ("train", "calibration", "test")[int(role[k7_idx])]
                    other_rates = [float(covered[test_acc == acc].mean())
                                   for acc in np.unique(test_acc) if acc != k7_idx]
                    mechanism[k7_role].append({
                        "half_width": half_width,
                        "other_test_protein_coverage": float(np.mean(other_rates)) if other_rates else None,
                    })

        level_result = {
            "residue_weighted": aggregate(overall_residue),
            "protein_averaged": aggregate(overall_protein),
            "by_stratum_residue_weighted": {
                s: aggregate(values) for s, values in sorted(per_stratum_residue.items())
            },
            "by_stratum_protein_averaged": {
                s: aggregate(values) for s, values in sorted(per_stratum_protein.items())
            },
            "miss_concentration_across_test_predictions": {
                "n_test_predictions": test_rows_evaluated,
                "n_misses": int(sum(misses_by_accession.values())),
                "K7PQ54_misses": int(misses_by_accession.get("K7PQ54", 0)),
                "K7PQ54_share": float(misses_by_accession.get("K7PQ54", 0) / sum(misses_by_accession.values()))
                    if sum(misses_by_accession.values()) else None,
                "top_10_accessions_by_miss_count": [
                    {"accession": acc, "misses": count,
                     "share_of_misses": count / sum(misses_by_accession.values())}
                    for acc, count in misses_by_accession.most_common(10)
                ],
            },
        }
        result["levels"][str(level)] = level_result

    if mechanism_accession is not None:
        result["mechanism_K7PQ54"] = {
            role: {
                "n_repeats": len(vals),
                "mean_half_width": float(np.mean([r["half_width"] for r in vals])) if vals else None,
                "mean_other_test_protein_coverage": float(np.mean([
                    r["other_test_protein_coverage"] for r in vals
                    if r["other_test_protein_coverage"] is not None
                ])) if any(r["other_test_protein_coverage"] is not None for r in vals) else None,
            }
            for role, vals in mechanism.items()
        }
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoints", type=Path, nargs="+", required=True)
    ap.add_argument("--v2", type=Path, required=True,
                    help="V2 classes, independently sample-checked against cached raw HTTP responses")
    ap.add_argument("--registered-evaluation", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    all_rows = []
    exclusion_counts = Counter()
    inputs = []
    for path in args.checkpoints:
        raw = path.read_bytes()
        checkpoint = json.loads(raw)
        inputs.append({"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                       "batch_index": checkpoint.get("batch_index"),
                       "n_rows": len(checkpoint.get("rows", []))})
        for row in checkpoint.get("rows", []):
            eligible, reason = finite_eligible(row)
            if eligible:
                all_rows.append(row)
            else:
                exclusion_counts[reason or "unknown"] += 1

    v2 = json.loads(args.v2.read_text())
    training = {acc: rec["training"]["class"] for acc, rec in v2["accessions"].items()}
    no_k7_rows = [r for r in all_rows if r["accession"] != "K7PQ54"]
    full = reproduce(all_rows, (0.90, 0.95), training, "K7PQ54")
    no_k7 = reproduce(no_k7_rows, (0.90, 0.95), training, None)
    registered = json.loads(args.registered_evaluation.read_text())
    registered_vals = {}
    for lev in ("0.9", "0.95"):
        pooled = registered["evaluation_six_bin"]["aggregates"][lev]["marginal"]["pooled"]
        registered_vals[lev] = {
            "registered_mean_coverage": pooled["mean"],
            "registered_mcse": pooled["mcse"],
            "registered_n_repeats": pooled["r_eff"],
        }
    output = {
        "method": {
            "registration": "docs/e422_registration.md §3A.1",
            "split": "sorted eligible accessions; default_rng(seed).permutation(n); floor-half train, floor-half-of-remainder calibration, remainder test; seeds 0..199",
            "interval": "training median lDDT; absolute calibration residual; ceil((n_cal_rows+1)*(1-alpha)) order statistic; inclusive coverage",
            "strata": "prior_relative = SAME_ACCESSION or SEQ95; no_close_relative = SEQ30 or NOVEL; strata are descriptive",
            "mcse": "sample standard deviation across 200 repeat statistics divided by sqrt(R)",
            "numpy_version": np.__version__,
            "no_p6_p7": "only the pooled interval was computed; no per-bin conditional coverage and no Mondrian/CQR comparison",
        },
        "inputs": inputs + [{"path": str(args.v2), "sha256": sha(args.v2)},
                            {"path": str(args.registered_evaluation), "sha256": sha(args.registered_evaluation)}],
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
        "registered_evaluation_comparison": registered_vals,
        "including_K7PQ54": full,
        "excluding_K7PQ54_posthoc": no_k7,
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "inputs": inputs,
        "exclusions": dict(exclusion_counts),
        "including_K7PQ54": {
            lev: val["residue_weighted"] for lev, val in full["levels"].items()
        },
        "excluding_K7PQ54": {
            lev: {"residue_weighted": val["residue_weighted"], "protein_averaged": val["protein_averaged"]}
            for lev, val in no_k7["levels"].items()
        },
        "strata_including_K7PQ54": {
            lev: val["by_stratum_residue_weighted"] for lev, val in full["levels"].items()
        },
        "strata_excluding_K7PQ54": {
            lev: val["by_stratum_residue_weighted"] for lev, val in no_k7["levels"].items()
        },
        "mechanism": full["mechanism_K7PQ54"],
        "registered": registered_vals,
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
