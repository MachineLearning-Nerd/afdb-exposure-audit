"""Phase B (a)/(d): independent digest verification + allowed-field extraction.
Blinded: prints only digests, ledger/provenance fields, P5 fields, and binding
status strings. Never prints per-bin / floor-table / Mondrian / CQR / P6/P7 values.
No experiments/ imports; canonical rule re-implemented from its documented form
(json.dumps sort_keys=True, separators=(',',':'), UTF-8; digest over the object
minus its 'sha256' key)."""
import hashlib, json, os, sys
ROOT = sys.argv[1]
def canon(o):
    return hashlib.sha256(json.dumps(o, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode("utf-8")).hexdigest()
def check(path):
    raw = open(path, "rb").read()
    o = json.loads(raw.decode("utf-8"))
    body = {k: v for k, v in o.items() if k != "sha256"}
    calc = canon(body)
    return o, calc, o.get("sha256"), hashlib.sha256(raw).hexdigest()
out = {}
# checkpoints
ck = {}
total = 0
for k in range(1, 5):
    p = os.path.join(ROOT, "results/e427/batches", f"e422_batch{k:02d}_checkpoint.json")
    o, calc, rec, filesha = check(p)
    ck[k] = {"recorded": rec, "recomputed": calc, "match": calc == rec,
             "n_rows": len(o["rows"]), "file_sha256": filesha,
             "e427_relabel_present": "e427_relabel" in o,
             "manifest_checkpoint_sha256": o.get("manifest_checkpoint_sha256")}
    total += len(o["rows"])
out["checkpoints"] = ck
out["checkpoint_rows_total"] = total
# evaluation artifact
o, calc, rec, filesha = check(os.path.join(ROOT, "results/e427/evals/e422_batch04_evaluation.json"))
ev = o["evaluation_six_bin"]
p5 = ev.get("predictions", {}).get("p5")
pooled = {lv: ev["aggregates"][lv]["marginal"]["pooled"] for lv in ev["aggregates"]}
out["evaluation"] = {
    "recorded": rec, "recomputed": calc, "match": calc == rec, "file_sha256": filesha,
    "schema": o.get("schema"), "batch_index": o.get("batch_index"), "role": o.get("role"),
    "t0_utc": o.get("t0_utc"), "evaluated_utc": o.get("evaluated_utc"),
    "ledger": o.get("ledger"), "n_rows_cumulative": o.get("n_rows_cumulative"),
    "prior_binding_path": o.get("prior_binding_path"),
    "floor_decision_is_null": o.get("floor_decision") is None,
    "evaluation_on_decided_grid_is_null": o.get("evaluation_on_decided_grid") is None,
    "inner_scope": ev.get("scope"), "inner_status": ev.get("status"),
    "inner_registration": ev.get("registration"), "inner_seeds": ev.get("seeds"),
    "inner_levels": ev.get("levels"), "inner_n_rows_input": ev.get("n_rows_input"),
    "inner_n_accessions": ev.get("n_accessions"),
    "inner_exclusions": ev.get("exclusions"),
    "canonical_accessions_sha256": canon(ev.get("canonical_accessions")),
    "inner_keys": sorted(ev.keys()),
    "inner_prediction_keys": sorted(ev.get("predictions", {}).keys()),
    "aggregate_keys_per_level": {lv: sorted(ev["aggregates"][lv].keys()) for lv in ev["aggregates"]},
    "p5": p5, "marginal_pooled": pooled,
    "binding_view_status_strings": {k: v.get("outcome") if k.startswith("p5") is False else v
                                    for k, v in o.get("binding", {}).items()},
    "top_keys": sorted(o.keys()),
}
# receipt
o, calc, rec, filesha = check(os.path.join(ROOT, "results/e427/evals/e422_batch04_execution_receipt.json"))
out["receipt"] = {"recorded": rec, "recomputed": calc, "match": calc == rec, "file_sha256": filesha,
                  **{k: o.get(k) for k in ("schema", "batch_index", "role", "t0_utc", "evaluated_utc",
                                           "started_utc", "completed_utc", "wall_seconds",
                                           "registered_budget_hours", "within_budget", "deviation_note")}}
# binding state
o, calc, rec, filesha = check(os.path.join(ROOT, "results/e427/batches/e422_binding_state.json"))
out["binding_state"] = {"recorded": rec, "recomputed": calc, "match": calc == rec, "file_sha256": filesha,
                        "schema": o.get("schema"), "batch_index": o.get("batch_index"),
                        "binding_keys": sorted(o.get("binding", {}).keys()),
                        "binding": {k: v for k, v in o.get("binding", {}).items() if k in ("p5_90", "p5_95")},
                        "non_p5_keys_present": [k for k in o.get("binding", {}) if k not in ("p5_90", "p5_95")]}
json.dump(out, open(os.path.join(ROOT, "results/e427/audit/phaseB_subagent/digests.json"), "w"), indent=1, sort_keys=True)
print(json.dumps(out, indent=1, sort_keys=True))
