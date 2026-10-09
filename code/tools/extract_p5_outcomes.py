"""Extract ONLY the registered P5-90/P5-95 outcome blocks from the batch-4
evaluation artifacts (original e422 run and corrected e427 rerun).

The full evaluation artifacts also contain pLDDT-binned / six-bin /
decided-grid sections that feed predictions P6/P7. Those predictions are
registered to stay blind until the census's batch-8 evaluation, so the full
artifacts are NOT released here (they will be released after the census
closes). This script reads each artifact, selects the P5 keys by name,
recomputes the artifact's canonical self-digest (SHA-256 of canonical JSON
without the `sha256` field) and the file SHA-256, and writes nothing else.
It never prints or stores any other section.

Usage (run once by the authors against the unreleased artifacts):
  python -I code/tools/extract_p5_outcomes.py <orig_eval.json> <corr_eval.json> <out.json>
"""
import hashlib
import json
import sys

P5_KEYS = ("prediction", "outcome", "mean_cov", "mcse", "r_eff", "threshold", "nominal")


def canonical_sha256(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def extract(path, label):
    raw = open(path, "rb").read()
    art = json.loads(raw)
    body = {k: v for k, v in art.items() if k != "sha256"}
    out = {
        "run": label,
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "stored_canonical_sha256": art["sha256"],
        "recomputed_canonical_sha256_matches": canonical_sha256(body) == art["sha256"],
        "schema": art.get("schema"),
        "batch_index": art.get("batch_index"),
        "evaluated_utc": art.get("evaluated_utc"),
        "n_rows_cumulative": art.get("n_rows_cumulative"),
    }
    for key in ("p5_90", "p5_95"):
        blk = art["binding"][key]
        out[key] = {"bound_at_batch": blk["bound_at_batch"],
                    "outcome": {k: blk["outcome"][k] for k in P5_KEYS}}
    return out


def main(orig, corr, dest):
    res = {
        "schema": "p5-batch04-outcomes-extract-v1",
        "note": ("Key-selected extract of the registered P5 outcome blocks only. The full batch-4 "
                 "evaluation artifacts contain sections bound to predictions P6/P7, which stay blind "
                 "until the batch-8 evaluation; they will be released after the census closes."),
        "original_e422": extract(orig, "original labels (results/e422/evals/e422_batch04_evaluation.json)"),
        "corrected_e427": extract(corr, "corrected labels (results/e427/evals/e422_batch04_evaluation.json)"),
    }
    with open(dest, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
        fh.write("\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
