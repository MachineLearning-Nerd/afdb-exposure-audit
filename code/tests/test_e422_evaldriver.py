"""Offline tests for the e422 section 6 evaluation-schedule driver."""

from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_evaldriver as driver  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e422_protocol as protocol  # noqa: E402

T0 = "2026-09-10T00:00:00Z"


# --------------------------------------------------------------------------- #
# Helpers.
# --------------------------------------------------------------------------- #
def write_checkpoint(out_dir, batch, rows, t0=T0, extra=None):
    checkpoint = {
        "schema": "e422-batch-checkpoint-v1",
        "e_number": "e422",
        "t0_utc": t0,
        "batch_index": batch,
        "window": {"start": f"2026-09-{9 + 7 * (batch - 1):02d}T00:00:00Z",
                   "close": f"2026-09-{16 + 7 * (batch - 1):02d}T00:00:00Z"},
        "completed_utc": f"2026-10-{9 + batch:02d}T00:00:00Z",
        "census_entries": 1,
        "entries_with_mappings": 1,
        "manifest_checkpoint_sha256": "0" * 64,
        "manifest_admitted_new": [],
        "freeze_receipt_sha256": "1" * 64,
        "n_rows": len(rows),
        "rows": rows,
        "exclusions": {"mapping_failures": [], "manifest_failures": [],
                       "triple_failures": []},
        "verified_accessions": [],
        "receipts": [],
    }
    if extra:
        checkpoint.update(extra)
    checkpoint["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in checkpoint.items() if k != "sha256"})
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"e422_batch{batch:02d}_checkpoint.json")
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(checkpoint) + b"\n")
    return checkpoint


def synth_rows(n_accessions=40, per_acc=40, seed=7, starve=None):
    """Rows spread across all six bins (10/10/12/16/22/30) so an 8-batch
    cumulative ledger clears the registered floors on the six-bin grid."""
    rng = __import__("numpy").random.default_rng(seed)
    rows = []
    for a in range(n_accessions):
        accession = f"P{a:05d}"
        for r in range(per_acc):
            u = rng.uniform()
            if u < 0.10:
                p = rng.uniform(0, 50)
            elif u < 0.20:
                p = rng.uniform(50, 60)
            elif u < 0.32:
                p = rng.uniform(60, 70)
            elif u < 0.48:
                p = rng.uniform(70, 80)
            elif u < 0.70:
                p = rng.uniform(80, 90)
            else:
                p = rng.uniform(90, 100)
            y = min(max(1.0 - p / 150.0 + rng.normal(0, 0.05), 0.0), 1.0)
            rows.append({"entry": "aaaa", "chain": "A", "accession": accession,
                         "uniprot_pos": r + 1, "plddt": round(p, 4),
                         "lddt": round(y, 4), "d_kabsch": 1.0,
                         "b_factor_z": 0.5})
    # (rows carry the descriptive channels: eligibility requires them)
    return rows


# --------------------------------------------------------------------------- #
# Ledger loading.
# --------------------------------------------------------------------------- #
def test_load_ledger_verifies_digests(tmp_path):
    rows = synth_rows(n_accessions=4, per_acc=5)
    write_checkpoint(str(tmp_path), 1, rows)
    loaded, ledger = driver.load_cumulative_ledger(str(tmp_path), 1)
    assert len(loaded) == len(rows)
    assert ledger[0]["batch_index"] == 1
    with pytest.raises(driver.DriverError) as excinfo:
        driver.load_cumulative_ledger(str(tmp_path), 2)
    assert excinfo.value.code == "E422_CHECKPOINT_MISSING"


def test_load_ledger_detects_tamper(tmp_path):
    checkpoint = write_checkpoint(str(tmp_path), 1, synth_rows(n_accessions=4,
                                                               per_acc=5))
    path = os.path.join(str(tmp_path), "e422_batch01_checkpoint.json")
    stored = json.loads(open(path, "rb").read().decode())
    stored["rows"][0]["plddt"] = 42.0
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(stored) + b"\n")
    with pytest.raises(driver.DriverError) as excinfo:
        driver.load_cumulative_ledger(str(tmp_path), 1)
    assert excinfo.value.code == "E422_CHECKPOINT_DIGEST_MISMATCH"
    assert checkpoint  # the original dict is intact in memory


# --------------------------------------------------------------------------- #
# Floor table + decision.
# --------------------------------------------------------------------------- #
def test_floor_table_counts_and_exclusions():
    rows = [
        {"accession": "P00001", "plddt": 45.0, "lddt": 0.7,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00001", "plddt": 45.0, "lddt": 0.8,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00002", "plddt": 55.0, "lddt": 0.6,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00003", "plddt": 95.0, "lddt": 0.9,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00003", "plddt": 95.0, "lddt": None,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00004", "plddt": 500.0, "lddt": 0.5,
         "d_kabsch": 1.0, "b_factor_z": 0.0},
        {"accession": "P00004", "plddt": 75.0, "lddt": 0.5,
         "d_kabsch": None, "b_factor_z": 0.0},
    ]
    table, excluded, _sets = driver.floor_table(rows)
    assert [b["residues"] for b in table] == [2, 1, 0, 0, 0, 1]
    assert [b["proteins"] for b in table] == [1, 1, 0, 0, 0, 1]
    assert excluded == {"null_label": 1, "nonfinite_label": 0,
                        "null_feature": 0, "nonfinite_feature": 0,
                        "feature_out_of_grid": 1,
                        "null_descriptive_channel": 1}


def _table_from_counts(counts):
    """counts: list of (residues, proteins) per six-bin."""
    return [{"bin_lo": lo, "bin_hi": hi, "residues": r, "proteins": p}
            for (lo, hi), (r, p) in zip(protocol.PLDDT_BINS, counts)]


def test_floor_decision_six_bin_when_all_pass():
    counts = [(1200, 25)] * 6
    decision = driver.floor_decision(_table_from_counts(counts),
                                     [set(range(25))] * 6)
    assert decision["decision"] == driver.SIX_BIN
    assert decision["status"] == "FLOORS_MET"


def test_floor_decision_step1_on_mid_bin_starvation():
    # 50-60 starved on residues; 0-50 fine -> step 1 merges 50-60 + 60-70
    counts = [(1200, 25), (999, 25), (500, 25), (1200, 25), (2000, 30),
              (5000, 40)]
    sets = [set(range(25)), set(range(25)), set(range(25, 45)),
            set(range(25)), set(range(30)), set(range(40))]
    decision = driver.floor_decision(_table_from_counts(counts), sets)
    assert decision["decision"] == driver.STEP1
    merged = [b for b in decision["floor_table_on_decided_grid"]
              if b["bin_lo"] == 50.0][0]
    assert merged["residues"] == 999 + 500
    assert merged["proteins"] == 45  # union of {0..24} and {25..44}


def test_floor_decision_step2_on_low_bin_starvation():
    # 0-50 starved on proteins -> step 1's [0,50) bin still fails -> step 2
    counts = [(1200, 19), (1200, 25), (500, 25), (1200, 25), (2000, 30),
              (5000, 40)]
    decision = driver.floor_decision(_table_from_counts(counts),
                                     [set(range(19))] + [set(range(25))] * 2
                                     + [set(range(25))] + [set(range(30))]
                                     + [set(range(40))])
    assert decision["decision"] == driver.STEP2
    low = decision["floor_table_on_decided_grid"][0]
    assert low["bin_lo"] == 0.0 and low["bin_hi"] == 80.0
    assert low["residues"] == 1200 + 1200 + 500 + 1200


def test_floor_decision_top3_when_nothing_passes():
    counts = [(10, 2)] * 6
    decision = driver.floor_decision(_table_from_counts(counts),
                                     [{f"P{i:05d}"} for i in range(2)] * 6)
    assert decision["decision"] == driver.TOP3
    assert decision["status"] == "FLOORS_UNMET_TOP3_RESCOPE"
    assert decision["grid"] == [[70.0, 80.0], [80.0, 90.0], [90.0, 100.01]]


# --------------------------------------------------------------------------- #
# use_grid (protocol merge fail-safe).
# --------------------------------------------------------------------------- #
def test_use_grid_roundtrip_and_binning():
    restore = protocol.use_grid(driver.MERGED_STEP1)
    try:
        assert protocol._to_bins(__import__("numpy").array([55.0]))[0] == 1
        assert protocol._to_bins(__import__("numpy").array([45.0]))[0] == 0
    finally:
        restore()
    assert protocol._to_bins(__import__("numpy").array([55.0]))[0] == 1  # six-bin
    assert protocol._to_bins(__import__("numpy").array([65.0]))[0] == 2


def test_use_grid_validation():
    with pytest.raises(protocol.E422SpecViolation):
        protocol.use_grid(((0.0, 50.0), (60.0, 100.01)))  # not contiguous
    with pytest.raises(protocol.E422SpecViolation):
        protocol.use_grid(((0.0, 50.0), (50.0, 99.0)))  # wrong top edge
    with pytest.raises(protocol.E422SpecViolation):
        protocol.use_grid(())


def test_use_grid_restricted_top3_unstratifies_low_rows():
    import numpy as np
    restore = protocol.use_grid(driver.TOP3_RESCOPE)
    try:
        idx = protocol._to_bins(np.array([45.0, 75.0, 95.0]))
        assert list(idx) == [-1, 0, 2]
        stats = protocol._bin_stats(np.array([45.0, 75.0, 95.0]),
                                    np.array([0.5, 0.6, 0.7]),
                                    np.array([0.0, 0.0, 0.0]),
                                    np.array([1.0, 1.0, 1.0]))
        assert set(stats) == {0, 1, 2}
        assert stats[0]["n_test"] == 1
    finally:
        restore()
    # eligibility authority is UNCHANGED by grid rebinds
    assert protocol.bin_index(45.0) == 0


# --------------------------------------------------------------------------- #
# Scheduled evaluations end-to-end (offline, small).
# --------------------------------------------------------------------------- #
def _run_batch4(tmp_path, rows):
    write_checkpoint(str(tmp_path), 1, synth_rows())
    write_checkpoint(str(tmp_path), 2, synth_rows(seed=8))
    write_checkpoint(str(tmp_path), 3, synth_rows(seed=9))
    write_checkpoint(str(tmp_path), 4, rows)
    return driver.run_scheduled_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-10-08T00:00:00Z")


def test_binding_map_prior_bindings_immutable_under_conflict():
    # Dated test amendment 2026-09-10 (w1:pG): FIRST-BINDING-WINS is the
    # registered semantics — a prior-bound key survives any conflicting
    # fresh evaluation, at every batch index.
    prior = {"p5_90": {"outcome": "BOUND_AT_DAY_28", "bound_at_batch": 4}}
    predictions = {"p5": {"0.9": "CONFLICTING_FRESH", "0.95": "NEW"},
                   "p6": "X", "p7a": "X", "p7b": "X"}
    state = driver.binding_map(8, predictions, prior)
    assert state["p5_90"] == {"outcome": "BOUND_AT_DAY_28",
                              "bound_at_batch": 4}
    assert state["p5_95"] == {"outcome": "NEW", "bound_at_batch": 8}


def test_binding_map_p6_p7_cannot_bind_before_batch8():
    # before batch 8 the P6/P7a/P7b channels are NOT EVALUABLE — their
    # prediction values are ignored even when non-null (no early bind)
    predictions = {"p5": {"0.9": None, "0.95": None},
                   "p6": "EARLY?", "p7a": "EARLY?", "p7b": "EARLY?"}
    state = driver.binding_map(4, predictions, {})
    assert set(state) == set()


def test_binding_map_at_close_fills_only_unbound():
    # the batch-12 close must never overwrite a real binding with
    # DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE — gaps only
    prior = {"p5_90": {"outcome": "B", "bound_at_batch": 4}}
    predictions = {"p5": {"0.9": "B", "0.95": "B"},
                   "p6": "C", "p7a": "C", "p7b": None}
    state = driver.binding_map(12, predictions, prior)
    assert state["p5_90"] == {"outcome": "B", "bound_at_batch": 4}
    assert state["p6"] == {"outcome": "C", "bound_at_batch": 12}
    assert state["p7b"]["outcome"] == driver.AT_CLOSE
    assert state["p7b"]["bound_at_batch"] == 12


def test_batch4_descriptive_p5_binding(tmp_path):
    artifact = _run_batch4(tmp_path, synth_rows(seed=10))
    assert artifact["role"] == "batch4_day28"
    assert artifact["floor_decision"] is None
    assert artifact["evaluation_on_decided_grid"] is None
    assert not os.path.exists(os.path.join(str(tmp_path / "results"),
                                           "e422_batch04_floor_decision.json"))
    assert artifact["sha256"] == manifest_mod.canonical_sha256(
        {k: v for k, v in artifact.items() if k != "sha256"})
    binding = artifact["binding"]
    assert binding["p6"]["outcome"] == driver.NOT_YET_EVALUABLE
    assert binding["p7a"]["outcome"] == driver.NOT_YET_EVALUABLE
    # P5 binds at batch 4 if the decision rule yields an outcome
    state = json.loads(open(str(tmp_path / "binding.json"), "rb").read())
    for key in ("p5_90", "p5_95"):
        if key in state["binding"]:
            assert state["binding"][key]["bound_at_batch"] == 4
            assert binding[key]["outcome"] == state["binding"][key]["outcome"]
    # six-bin floor table recorded descriptively
    assert len(artifact["six_bin_floor_table"]) == 6
    assert artifact["n_rows_cumulative"] > 0


def test_batch8_decision_receipt_and_binding(tmp_path):
    _run_batch4(tmp_path, synth_rows(seed=10))
    for k in range(5, 9):
        write_checkpoint(str(tmp_path), k, synth_rows(seed=10 + k))
    artifact = driver.run_scheduled_evaluation(
        T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-11-05T00:00:00Z")
    assert artifact["role"] == "batch8_day56"
    decision = artifact["floor_decision"]
    assert decision is not None
    receipt_path = os.path.join(str(tmp_path / "results"),
                                "e422_batch08_floor_decision.json")
    assert os.path.exists(receipt_path)
    receipt = json.loads(open(receipt_path, "rb").read())
    assert receipt["decision"] == decision["decision"]
    assert receipt["sha256"] == manifest_mod.canonical_sha256(
        {k: v for k, v in receipt.items() if k != "sha256"})
    binding = artifact["binding"]
    # with a healthy synthetic ledger all six floors pass -> six-bin binds
    assert decision["decision"] == driver.SIX_BIN
    assert artifact["evaluation_on_decided_grid"] is None
    # p5 bindings carried from batch 4 unchanged
    state = json.loads(open(str(tmp_path / "binding.json"), "rb").read())
    for key in ("p5_90", "p5_95"):
        if key in state["binding"]:
            assert state["binding"][key]["bound_at_batch"] == 4


def test_batch8_receipt_survives_evaluation_failure(tmp_path, monkeypatch):
    # Dated test amendment 2026-09-10 (w1:pG): pins the section 4 fence as
    # MECHANICS, not just happy-path ordering — the floor-decision receipt
    # is written BEFORE protocol.run_evaluation is ever called, so an
    # evaluation that crashes mid-flight must still leave a valid,
    # digest-verifying receipt on disk (decided-before-unblinding).
    _run_batch4(tmp_path, synth_rows(seed=10))
    for k in range(5, 9):
        write_checkpoint(str(tmp_path), k, synth_rows(seed=10 + k))
    receipt_path = os.path.join(str(tmp_path / "results"),
                                "e422_batch08_floor_decision.json")

    def _boom(*_args, **_kwargs):
        raise RuntimeError("evaluation crashed")

    monkeypatch.setattr(protocol, "run_evaluation", _boom)
    with pytest.raises(RuntimeError, match="evaluation crashed"):
        driver.run_scheduled_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=str(tmp_path / "binding.json"),
            now_iso="2026-11-05T00:00:00Z")
    assert os.path.exists(receipt_path)  # the fence held
    receipt = json.loads(open(receipt_path, "rb").read())
    assert receipt["decision"] == driver.SIX_BIN
    assert receipt["sha256"] == manifest_mod.canonical_sha256(
        {k: v for k, v in receipt.items() if k != "sha256"})
    assert receipt["note"].startswith("decided before any coverage "
                                      "evaluation exists")


def test_batch8_merged_grid_conditional_evaluation(tmp_path, monkeypatch):
    # force a merged decision by starving bins 1-2 in the ledger
    def starved_rows(**kw):
        rows = synth_rows(**kw)
        return [r for r in rows
                if not (50.0 <= r["plddt"] < 70.0)]
    write_checkpoint(str(tmp_path), 1, synth_rows())
    write_checkpoint(str(tmp_path), 2, synth_rows(seed=8))
    write_checkpoint(str(tmp_path), 3, synth_rows(seed=9))
    write_checkpoint(str(tmp_path), 4, synth_rows(seed=10))
    for k in range(5, 9):
        write_checkpoint(str(tmp_path), k, starved_rows(seed=10 + k))
    artifact = driver.run_scheduled_evaluation(
        T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-11-05T00:00:00Z")
    decision = artifact["floor_decision"]
    assert decision["decision"] in (driver.STEP1, driver.STEP2, driver.TOP3)
    assert artifact["evaluation_on_decided_grid"] is not None
    assert artifact["evaluation_on_decided_grid"]["grid_decision"] == \
        decision["decision"]
    # the decided-grid evaluation ran on the DECIDED grid
    assert artifact["evaluation_on_decided_grid"]["grid"] == decision["grid"]


def test_batch12_close_sets_at_close(tmp_path):
    # tiny ledger -> every evaluation DEFERRED -> close records AT_CLOSE
    for k in range(1, 13):
        write_checkpoint(str(tmp_path), k, synth_rows(n_accessions=2,
                                                      per_acc=5, seed=k))
    artifact = driver.run_scheduled_evaluation(
        T0, 12, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-12-03T00:00:00Z")
    assert artifact["role"] == "batch12_day84_close"
    for key in ("p5_90", "p5_95", "p6", "p7a", "p7b"):
        assert artifact["binding"][key]["outcome"] == driver.AT_CLOSE
        assert artifact["binding"][key]["bound_at_batch"] == 12


def test_unscheduled_batch_rejected(tmp_path):
    with pytest.raises(driver.DriverError) as excinfo:
        driver.run_scheduled_evaluation(T0, 5, str(tmp_path),
                                        str(tmp_path / "results"))
    assert excinfo.value.code == "E422_UNSCHEDULED_EVALUATION"


def test_batch4_conditional_outputs_suppressed(tmp_path):
    # registered section 6 scope: the batch-4 artifact carries P5 outcomes
    # + the descriptive floor table ONLY (decided-before-unblinding fence)
    artifact = _run_batch4(tmp_path, synth_rows(seed=10))
    ev = artifact["evaluation_six_bin"]
    assert ev["scope"] == "BATCH4_P5_ONLY_REGISTERED_SUPPRESSION"
    assert "p6" not in ev["predictions"] and "p7a" not in ev["predictions"]
    for lv in ev["aggregates"].values():
        assert "per_bin" not in lv
    assert "p6" not in ev and "p7a" not in ev and "p7b" not in ev
    assert "joint_grid" not in ev
    # P5 outcomes survive the scoping
    assert "p5" in ev["predictions"]
    assert ev["predictions"]["p5"]["0.9"]["prediction"] == "P5-90"


def test_floor_table_rejects_restricted_grid(tmp_path):
    # the floor trajectory is a PRE-decision artifact: on the registered
    # driver order it always runs before any use_grid; a restricted grid
    # bound at call time is a typed failure, never a silent misfold
    rows = synth_rows(n_accessions=4, per_acc=5)
    restore = protocol.use_grid(driver.MERGED_STEP1)
    try:
        with pytest.raises(driver.DriverError) as excinfo:
            driver.floor_table(rows)
        assert excinfo.value.code == "E422_FLOOR_TABLE_ON_RESTRICTED_GRID"
    finally:
        restore()
    # on the registered grid it works as before
    table, _excl, _sets = driver.floor_table(rows)
    assert len(table) == 6
