"""e422 section 6 evaluation-schedule driver.

Orchestrates the three registered cumulative evaluations (batch-4 day-28,
batch-8 day-56, batch-12 day-84 close) over the batch checkpoints
produced by experiments/e422_runner.py:

  1. load the cumulative label ledger (batches 1..k), verifying each
     batch checkpoint's canonical-JSON SHA-256 (typed tamper failure);
  2. record the six-bin floor trajectory (per-bin eligible residues and
     proteins) — DESCRIPTIVE at batch 4;
  3. at batch >= 8 (the binding point): make the registered section 4
     floor decision — SIX_BIN_GRID, MERGED_STEP1, MERGED_STEP2, or
     TOP3_RESCOPE — and write the decision receipt BEFORE any coverage
     evaluation exists (decided-before-unblinding fence);
  4. run the registered six-bin evaluation (protocol.run_evaluation);
     on a merged/restricted decision, run the conditional evaluation on
     the decided grid via protocol.use_grid (pool and splits unchanged;
     only the conditional strata change);
  5. apply the section 6 binding map: batch 4 binds P5-90/P5-95 only;
     batch >= 8 binds P6/P7a/P7b (on the decided grid) and any still
     unbound P5; batch 12 is the close — unbound predictions become
     DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE (a registered outcome);
  6. emit the driver artifact (canonical-JSON SHA-256, atomic write).

Offline: deterministic given the checkpoint files; no network.
"""

from __future__ import annotations

import json
import os
import sys
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e422_protocol as protocol  # noqa: E402


DRIVER_SCHEMA = "e422-evaldriver-v1"
DECISION_SCHEMA = "e422-floor-decision-v1"

SCHEDULE = {4: "batch4_day28", 8: "batch8_day56", 12: "batch12_day84_close"}

RESIDUE_FLOOR = 1000
PROTEIN_FLOOR = 20

SIX_BIN_GRID = protocol.PLDDT_BINS
MERGED_STEP1 = ((0.0, 50.0), (50.0, 70.0), (70.0, 80.0), (80.0, 90.0),
                (90.0, 100.01))
MERGED_STEP2 = ((0.0, 80.0), (80.0, 90.0), (90.0, 100.01))
TOP3_RESCOPE = ((70.0, 80.0), (80.0, 90.0), (90.0, 100.01))

SIX_BIN = "SIX_BIN_GRID"
STEP1 = "MERGED_STEP1"
STEP2 = "MERGED_STEP2"
TOP3 = "TOP3_RESCOPE"

DECISION_GRIDS = {SIX_BIN: SIX_BIN_GRID, STEP1: MERGED_STEP1,
                  STEP2: MERGED_STEP2, TOP3: TOP3_RESCOPE}

AT_CLOSE = "DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE"
NOT_YET_EVALUABLE = "DESCRIPTIVE_AT_BATCH4_FIRST_EVALUABLE_BATCH8"


class DriverError(RuntimeError):
    """Typed driver failure (never silent)."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


# --------------------------------------------------------------------------- #
# Ledger loading + digest verification.
# --------------------------------------------------------------------------- #
def load_cumulative_ledger(out_dir: str, through_batch: int) -> tuple[list[dict], list[dict]]:
    """Load batch checkpoints 1..k in order, verifying each digest."""
    rows: list[dict] = []
    ledger: list[dict] = []
    for k in range(1, through_batch + 1):
        path = os.path.join(out_dir, f"e422_batch{k:02d}_checkpoint.json")
        if not os.path.exists(path):
            raise DriverError("E422_CHECKPOINT_MISSING", path)
        with open(path, "rb") as handle:
            checkpoint = json.loads(handle.read().decode("utf-8"))
        if "sha256" not in checkpoint:
            raise DriverError("E422_CHECKPOINT_DIGEST_MISSING", path)
        body = {key: val for key, val in checkpoint.items() if key != "sha256"}
        if manifest_mod.canonical_sha256(body) != checkpoint["sha256"]:
            raise DriverError("E422_CHECKPOINT_DIGEST_MISMATCH", path)
        rows.extend(checkpoint["rows"])
        ledger.append({"batch_index": k,
                       "sha256": checkpoint["sha256"],
                       "n_rows": len(checkpoint["rows"]),
                       "manifest_checkpoint_sha256":
                           checkpoint.get("manifest_checkpoint_sha256")})
    return rows, ledger


# --------------------------------------------------------------------------- #
# Floor table + section 4 decision.
# --------------------------------------------------------------------------- #
def floor_table(rows: list[dict]) -> tuple[list[dict], dict]:
    """Six-bin floor trajectory over ELIGIBLE rows (residues + proteins),
    plus the exclusion tally.  Protein membership is carried exactly (per
    bin accession sets) so merged-bin floors fold without approximation;
    the artifact form carries counts only."""
    if protocol._ANALYSIS_GRID is not protocol.PLDDT_BINS:
        # the registered driver order always decides floors BEFORE any
        # grid rebind; a restricted grid here would silently fold
        # sub-low-edge rows (bin -1) into the last indexed bin
        raise DriverError("E422_FLOOR_TABLE_ON_RESTRICTED_GRID",
                          "floor_table requires the registered six-bin grid")
    accessions, by_acc, excluded = protocol._prep_rows(rows)
    residues = [0] * len(protocol.PLDDT_BINS)
    proteins: list[set] = [set() for _ in protocol.PLDDT_BINS]
    for accession in accessions:
        ps, _ys, _zs = by_acc[accession]
        for bi in protocol._to_bins(ps):
            residues[int(bi)] += 1
            proteins[int(bi)].add(accession)
    table = [{"bin_lo": lo, "bin_hi": hi, "residues": residues[i],
              "proteins": len(proteins[i])}
             for i, (lo, hi) in enumerate(protocol.PLDDT_BINS)]
    return table, excluded, proteins


def _floors_pass(table: list[dict]) -> bool:
    return all(b["residues"] >= RESIDUE_FLOOR and b["proteins"] >= PROTEIN_FLOOR
               for b in table)


def floor_decision(six_table: list[dict], six_protein_sets: list[set]) -> dict:
    """Registered section 4 decision: six-bin -> step 1 -> step 2 -> top-3.

    Floors on merged/restricted grids fold EXACTLY: a merged bin's
    residues are the constituent sum; its protein count is the size of
    the union of constituent per-bin accession sets (a protein with rows
    in several constituents counts once).
    """
    def fold(grid):
        merged = []
        for lo, hi in grid:
            idx = [i for i, (blo, bhi) in enumerate(protocol.PLDDT_BINS)
                   if blo >= lo and bhi <= hi]
            merged.append({"bin_lo": lo, "bin_hi": hi,
                           "residues": sum(six_table[i]["residues"] for i in idx),
                           "proteins": len(set().union(*(six_protein_sets[i]
                                                         for i in idx)))
                           if idx else 0})
        return merged

    for name in (SIX_BIN, STEP1, STEP2, TOP3):
        grid = DECISION_GRIDS[name]
        table = six_table if name == SIX_BIN else fold(grid)
        if _floors_pass(table):
            return {"decision": name, "grid": [list(b) for b in grid],
                    "floor_table_on_decided_grid": table, "status": "FLOORS_MET"}
    # no candidate passes: registered last fallback is the top-3 re-scope
    # with low bins descriptive — the conditional analysis is scoped there
    # regardless of its own floors (registered text, section 4).
    grid = DECISION_GRIDS[TOP3]
    return {"decision": TOP3, "grid": [list(b) for b in grid],
            "floor_table_on_decided_grid": fold(grid),
            "status": "FLOORS_UNMET_TOP3_RESCOPE"}


# --------------------------------------------------------------------------- #
# Binding map.
# --------------------------------------------------------------------------- #
def _outcome_or_none(outcome: Any) -> Any:
    if isinstance(outcome, str) and outcome == protocol.DEFERRED:
        return None
    return outcome


def _scope_batch4(evaluation: dict) -> dict:
    """Registered batch-4 scope (section 6 checkpoint roles): the batch-4
    artifact carries P5-90/P5-95 outcomes + the descriptive floor table
    ONLY — all conditional-coverage outputs (per-bin tables, P6/P7
    predictions) are SUPPRESSED until the batch-8 checkpoint, so no
    conditional table exists before the count-based grid decision
    (decided-before-unblinding fence, section 4)."""
    import copy
    scoped = copy.deepcopy(evaluation)
    aggregates = scoped.get("aggregates", {})
    for lv in aggregates:
        aggregates[lv].pop("per_bin", None)
    preds = scoped.get("predictions", {})
    for key in ("p6", "p7a", "p7b"):
        preds.pop(key, None)
    # the §3B joint grid is likewise a conditional-coverage output:
    # suppressed at batch-4 with the per-bin tables
    scoped.pop("joint_grid", None)
    scoped["scope"] = "BATCH4_P5_ONLY_REGISTERED_SUPPRESSION"
    return scoped


def binding_map(batch_index: int, predictions: dict, prior: dict) -> dict:
    """Apply the section 6 binding schedule; returns the updated state.

    prior: {prediction_key: {"outcome": ..., "bound_at_batch": k}}
    """
    state = dict(prior)
    evaluable = batch_index >= 8
    for key in ("p5_90", "p5_95"):
        if key in state:
            continue
        raw = predictions["p5"]["0.9" if key == "p5_90" else "0.95"]
        outcome = _outcome_or_none(raw)
        if outcome is not None:
            state[key] = {"outcome": outcome, "bound_at_batch": batch_index}
    if evaluable:
        for key in ("p6", "p7a", "p7b"):
            if key in state:
                continue
            outcome = _outcome_or_none(predictions[key])
            if outcome is not None:
                state[key] = {"outcome": outcome, "bound_at_batch": batch_index}
    if batch_index >= 12:
        for key in ("p5_90", "p5_95", "p6", "p7a", "p7b"):
            if key not in state:
                state[key] = {"outcome": AT_CLOSE, "bound_at_batch": batch_index}
    return state


def _binding_view(state: dict, batch_index: int) -> dict:
    """Per-prediction view for the artifact at this checkpoint."""
    view = {}
    for key in ("p5_90", "p5_95", "p6", "p7a", "p7b"):
        if key in state:
            view[key] = dict(state[key])
        elif batch_index < 8 and key in ("p6", "p7a", "p7b"):
            view[key] = {"outcome": NOT_YET_EVALUABLE}
        else:
            view[key] = {"outcome": "NON_OUTCOME_AT_THIS_CHECKPOINT"}
    return view


# --------------------------------------------------------------------------- #
# Scheduled evaluation.
# --------------------------------------------------------------------------- #
def run_scheduled_evaluation(
    t0_utc: str,
    batch_index: int,
    out_dir: str,
    result_dir: str,
    *,
    seeds=None,
    prior_binding_path: str | None = None,
    now_iso: str | None = None,
) -> dict:
    """Run the registered evaluation for batch_index in {4, 8, 12}."""
    if batch_index not in SCHEDULE:
        raise DriverError("E422_UNSCHEDULED_EVALUATION", f"batch {batch_index}")
    seeds = seeds if seeds is not None else range(200)
    os.makedirs(result_dir, exist_ok=True)
    now = now_iso or census._utc_now()

    rows, ledger = load_cumulative_ledger(out_dir, batch_index)
    six_table, excluded, six_sets = floor_table(rows)

    # ---- floor decision at the binding point (BEFORE any evaluation) -----
    decision = None
    if batch_index >= 8:
        decision = floor_decision(six_table, six_sets)
        receipt = {
            "schema": DECISION_SCHEMA,
            "e_number": "e422",
            "t0_utc": t0_utc,
            "batch_index": batch_index,
            "decided_utc": now,
            "six_bin_floor_table": six_table,
            "residue_floor": RESIDUE_FLOOR,
            "protein_floor": PROTEIN_FLOOR,
            **decision,
            "note": ("decided before any coverage evaluation exists "
                     "(section 4 fence); binds P6/P7 grids from here on"),
        }
        receipt["sha256"] = manifest_mod.canonical_sha256(
            {k: v for k, v in receipt.items() if k != "sha256"})
        census._atomic_write_json(
            os.path.join(result_dir,
                         f"e422_batch{batch_index:02d}_floor_decision.json"),
            receipt)

    # ---- registered six-bin evaluation (always reported) -----------------
    evaluation = protocol.run_evaluation(rows, seeds=seeds)
    if batch_index == 4:
        evaluation = _scope_batch4(evaluation)

    # ---- conditional evaluation on the decided grid (merged/restricted) --
    merged_evaluation = None
    if decision is not None and decision["decision"] != SIX_BIN:
        restore = protocol.use_grid(
            tuple(tuple(b) for b in decision["grid"]))
        try:
            merged_evaluation = protocol.run_evaluation(rows, seeds=seeds)
        finally:
            restore()
        merged_evaluation["grid_decision"] = decision["decision"]

    # ---- binding state ----------------------------------------------------
    prior = {}
    if prior_binding_path and os.path.exists(prior_binding_path):
        with open(prior_binding_path, "rb") as handle:
            prior = json.loads(handle.read().decode("utf-8")).get("binding", {})
    source = merged_evaluation or evaluation
    predictions = source.get("predictions") or {
        "p5": {"0.9": protocol.DEFERRED, "0.95": protocol.DEFERRED},
        "p6": protocol.DEFERRED, "p7a": protocol.DEFERRED,
        "p7b": protocol.DEFERRED}
    state = binding_map(batch_index, predictions, prior)

    artifact = {
        "schema": DRIVER_SCHEMA,
        "e_number": "e422",
        "t0_utc": t0_utc,
        "batch_index": batch_index,
        "role": SCHEDULE[batch_index],
        "evaluated_utc": now,
        "ledger": ledger,
        "n_rows_cumulative": len(rows),
        "eligibility_exclusions": excluded,
        "six_bin_floor_table": six_table,
        "floor_decision": decision,
        "evaluation_six_bin": evaluation,
        "evaluation_on_decided_grid": merged_evaluation,
        "binding": _binding_view(state, batch_index),
        "prior_binding_path": prior_binding_path,
    }
    artifact["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in artifact.items() if k != "sha256"})
    census._atomic_write_json(
        os.path.join(result_dir,
                     f"e422_batch{batch_index:02d}_evaluation.json"),
        artifact)

    # ---- binding-state sidecar (input to the next scheduled evaluation) --
    if prior_binding_path:
        state_artifact = {
            "schema": "e422-binding-state-v1",
            "batch_index": batch_index,
            "binding": state,
        }
        state_artifact["sha256"] = manifest_mod.canonical_sha256(
            {k: v for k, v in state_artifact.items() if k != "sha256"})
        census._atomic_write_json(prior_binding_path, state_artifact)
    return artifact


__all__ = [
    "DRIVER_SCHEMA",
    "DECISION_SCHEMA",
    "SCHEDULE",
    "RESIDUE_FLOOR",
    "PROTEIN_FLOOR",
    "DECISION_GRIDS",
    "DriverError",
    "load_cumulative_ledger",
    "floor_table",
    "floor_decision",
    "binding_map",
    "run_scheduled_evaluation",
    "AT_CLOSE",
    "NOT_YET_EVALUABLE",
]
