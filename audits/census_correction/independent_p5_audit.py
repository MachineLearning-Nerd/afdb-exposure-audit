#!/usr/bin/env python3
"""Independent pooled P5 recomputation; no e42x/moluq imports.

The evaluation artifact is parsed for its integrity digest, corrected-ledger
references, and only the P5 fields requested by the audit.  No conditional or
per-bin result is inspected or emitted.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LEDGER = ROOT / "results/e427/batches"
EVAL = ROOT / "results/e427/evals/e422_batch04_evaluation.json"
RECEIPT = ROOT / "results/e427/evals/e422_batch04_execution_receipt.json"
SIDECAR = LEDGER / "e422_binding_state.json"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_sha256(obj) -> str:
    # Independently reproduce the frozen canonical-JSON rule.
    blob = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256_bytes(blob)


def finite_number(row, key):
    value = row.get(key)
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def eligible_rows(checkpoints):
    by_accession = {}
    n_input = 0
    for checkpoint in checkpoints:
        for row in checkpoint["rows"]:
            n_input += 1
            y = finite_number(row, "lddt")
            p = finite_number(row, "plddt")
            d = finite_number(row, "d_kabsch")
            z = finite_number(row, "b_factor_z")
            # The registered grid is contiguous from 0 through 100.01;
            # check only membership, never form/report pLDDT strata.
            if y is None or p is None or not (0.0 <= p < 100.01) or d is None or z is None:
                continue
            by_accession.setdefault(row["accession"], []).append(y)
    accessions = sorted(by_accession)
    return accessions, by_accession, n_input


def pooled_repeat(by_accession, accessions, seed: int, level: float) -> float:
    n = len(accessions)
    if n < 4:
        raise ValueError(f"registered split requires >=4 accessions, got {n}")
    permutation = np.random.default_rng(seed).permutation(n)
    n_train = n // 2
    n_cal = (n - n_train) // 2
    train_indices = permutation[:n_train]
    cal_indices = permutation[n_train:n_train + n_cal]
    test_indices = permutation[n_train + n_cal:]

    train = np.asarray([y for i in train_indices for y in by_accession[accessions[int(i)]]], dtype=float)
    calibration = np.asarray([y for i in cal_indices for y in by_accession[accessions[int(i)]]], dtype=float)
    test = np.asarray([y for i in test_indices for y in by_accession[accessions[int(i)]]], dtype=float)
    if not train.size or not calibration.size or not test.size:
        raise ValueError(f"empty registered split at seed {seed}")
    center = float(np.median(train))
    scores = np.sort(np.abs(calibration - center))
    # alpha_c = 1-level; conformal order statistic is ceil((n+1)*level).
    k = math.ceil((scores.size + 1) * level)
    if k > scores.size:
        raise ValueError(f"registered order statistic exceeds calibration n at seed {seed}")
    radius = float(scores[k - 1])
    low, high = center - radius, center + radius
    return float(np.mean((low <= test) & (test <= high)))


def summarize_repeat_coverages(values, level):
    values = np.asarray(values, dtype=np.float64)
    if values.size != 200:
        raise ValueError(f"expected 200 repeats, got {values.size}")
    mean = float(np.mean(values))
    mcse = float(np.std(values, ddof=1) / math.sqrt(values.size))
    threshold = level - 3.0 * mcse
    margin = mean - threshold
    return {"nominal": level, "r_eff": int(values.size), "mean_cov": mean,
            "mcse": mcse, "threshold": threshold, "margin": margin,
            "margin_in_mcse": margin / mcse if mcse else None,
            "outcome": "HOLD" if mean >= threshold else "FAIL"}


def verified_document(path: Path):
    doc = read_json(path)
    stored = doc.get("sha256")
    actual = canonical_sha256({k: v for k, v in doc.items() if k != "sha256"})
    return doc, {"stored_canonical_sha256": stored, "computed_canonical_sha256": actual,
                 "canonical_digest_matches": stored == actual,
                 "file_sha256": sha256_bytes(path.read_bytes())}


def main():
    evaluation, eval_digest = verified_document(EVAL)
    receipt, receipt_digest = verified_document(RECEIPT)
    sidecar, sidecar_digest = verified_document(SIDECAR)

    checkpoints, checkpoint_digest_rows = [], []
    for index in range(1, 5):
        path = LEDGER / f"e422_batch{index:02d}_checkpoint.json"
        checkpoint = read_json(path)
        stored = checkpoint.get("sha256")
        computed = canonical_sha256({k: v for k, v in checkpoint.items() if k != "sha256"})
        checkpoint_digest_rows.append({"batch_index": index, "stored_sha256": stored,
                                       "computed_sha256": computed,
                                       "digest_matches": stored == computed,
                                       "n_rows": len(checkpoint.get("rows", []))})
        checkpoints.append(checkpoint)

    accessions, by_accession, n_input = eligible_rows(checkpoints)
    recomputed = {}
    artifact_p5 = evaluation["evaluation_six_bin"]["predictions"]["p5"]
    aggregate_p5 = evaluation["evaluation_six_bin"]["aggregates"]
    p5_binding = {key: evaluation["binding"][key] for key in ("p5_90", "p5_95")}
    sidecar_binding = sidecar.get("binding", {})

    for level, p5_key in ((0.90, "0.9"), (0.95, "0.95")):
        repeat_values = [pooled_repeat(by_accession, accessions, seed, level)
                         for seed in range(200)]
        computed = summarize_repeat_coverages(repeat_values, level)
        claimed = artifact_p5[p5_key]
        agg = aggregate_p5[p5_key]["marginal"]["pooled"]
        bind_key = "p5_90" if level == 0.90 else "p5_95"
        bound = p5_binding[bind_key]
        side = sidecar_binding.get(bind_key, {})
        computed["artifact_match"] = (
            claimed.get("outcome") == computed["outcome"]
            and claimed.get("r_eff") == computed["r_eff"]
            and math.isclose(float(claimed.get("mean_cov")), computed["mean_cov"], rel_tol=0, abs_tol=1e-15)
            and math.isclose(float(claimed.get("mcse")), computed["mcse"], rel_tol=0, abs_tol=1e-15)
            and math.isclose(float(claimed.get("threshold")), computed["threshold"], rel_tol=0, abs_tol=1e-15)
            and math.isclose(float(agg.get("mean")), computed["mean_cov"], rel_tol=0, abs_tol=1e-15)
            and math.isclose(float(agg.get("mcse")), computed["mcse"], rel_tol=0, abs_tol=1e-15)
            and bound.get("bound_at_batch") == 4
            and bound.get("outcome", {}).get("outcome") == computed["outcome"]
            and math.isclose(float(bound.get("outcome", {}).get("mean_cov")), computed["mean_cov"], rel_tol=0, abs_tol=1e-15)
            and side.get("bound_at_batch") == 4
            and side.get("outcome", {}).get("outcome") == computed["outcome"]
            and math.isclose(float(side.get("outcome", {}).get("mean_cov")), computed["mean_cov"], rel_tol=0, abs_tol=1e-15)
        )
        recomputed[f"P5-{int(level * 100)}"] = computed

    expected_ledger = []
    for i, checkpoint in enumerate(checkpoints, start=1):
        expected_ledger.append({"batch_index": i, "sha256": checkpoint.get("sha256"),
                                "n_rows": len(checkpoint.get("rows", [])),
                                "manifest_checkpoint_sha256": checkpoint.get("manifest_checkpoint_sha256")})
    ledger_matches = evaluation.get("ledger") == expected_ledger
    sidecar_keys = sorted(sidecar_binding)
    p6p7_unbound = not any(key in sidecar_binding for key in ("p6", "p7a", "p7b"))
    actual_seed_span = evaluation["evaluation_six_bin"].get("seeds")

    receipt_links = {
        "same_batch": receipt.get("batch_index") == evaluation.get("batch_index") == 4,
        "same_t0": receipt.get("t0_utc") == evaluation.get("t0_utc"),
        "same_role": receipt.get("role") == evaluation.get("role"),
        "same_evaluated_utc": receipt.get("evaluated_utc") == evaluation.get("evaluated_utc"),
        "within_registered_budget": receipt.get("within_budget") is True,
    }
    output = {
        "method": "independent NumPy/stdlib recomputation; no e42x or moluq imports",
        "evaluation": {"schema": evaluation.get("schema"), "e_number": evaluation.get("e_number"),
                       "batch_index": evaluation.get("batch_index"),
                       "n_rows_cumulative": evaluation.get("n_rows_cumulative"),
                       "seed_span": actual_seed_span,
                       "registered_seed_span_matches": actual_seed_span == [0, 199],
                       "prior_binding_path": evaluation.get("prior_binding_path"),
                       "canonical_digest": eval_digest},
        "execution_receipt": {"schema": receipt.get("schema"), "e_number": receipt.get("e_number"),
                              "batch_index": receipt.get("batch_index"),
                              "wall_seconds": receipt.get("wall_seconds"),
                              "registered_budget_hours": receipt.get("registered_budget_hours"),
                              "within_budget": receipt.get("within_budget"),
                              "links_to_evaluation": receipt_links,
                              "canonical_digest": receipt_digest},
        "binding_sidecar": {"schema": sidecar.get("schema"), "batch_index": sidecar.get("batch_index"),
                            "binding_keys": sidecar_keys, "p6_p7_unbound": p6p7_unbound,
                            "canonical_digest": sidecar_digest},
        "corrected_checkpoint_digests": checkpoint_digest_rows,
        "ledger_references_corrected_checkpoints": ledger_matches,
        "ledger_rows_sum_to_evaluation": sum(r["n_rows"] for r in expected_ledger) == evaluation.get("n_rows_cumulative"),
        "eligible_rows_for_p5": sum(len(v) for v in by_accession.values()),
        "eligible_accessions_for_p5": len(accessions),
        "split_sizes_accessions": {"train": len(accessions) // 2,
                                   "cal": (len(accessions) - len(accessions) // 2) // 2,
                                   "test": len(accessions) - len(accessions) // 2
                                           - (len(accessions) - len(accessions) // 2) // 2},
        "input_rows": n_input,
        "p5_recomputed": recomputed,
    }
    output["all_internal_digests_match"] = (
        eval_digest["canonical_digest_matches"]
        and receipt_digest["canonical_digest_matches"]
        and sidecar_digest["canonical_digest_matches"]
        and all(r["digest_matches"] for r in checkpoint_digest_rows))
    output["all_p5_match"] = all(v["artifact_match"] for v in recomputed.values())
    (HERE / "independent_p5_audit.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
