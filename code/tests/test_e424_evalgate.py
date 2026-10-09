"""Offline tests for the e424 binding-state integrity gate (2026-09-10).

The gate wraps the FROZEN e422_evaldriver; every refusal here must
precede the driver call (asserted via a sentinel monkeypatch).  Also
pins the batch-8 binding outcomes (P6/P7a/P7b bind on the decided grid)
that the frozen driver's own suite left unpinned.
"""

from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_evaldriver as driver  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e424_evalgate as gate  # noqa: E402

# the registered T0 bind: the gate pins t0 == REGISTERED_T0 (detached
# review 2026-09-10), so the offline ledgers must be built at the frozen
# gate-(a) instant, not a dev convenience T0
T0 = "2026-09-11T00:00:00Z"


# --------------------------------------------------------------------------- #
# Helpers (mirroring tests/test_e422_evaldriver.py fixtures).
# --------------------------------------------------------------------------- #
def write_checkpoint(out_dir, batch, rows, t0=T0):
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
    checkpoint["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in checkpoint.items() if k != "sha256"})
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"e422_batch{batch:02d}_checkpoint.json")
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(checkpoint) + b"\n")
    return checkpoint


def synth_rows(n_accessions=40, per_acc=40, seed=7):
    import numpy as np
    rng = np.random.default_rng(seed)
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
    return rows


def make_sidecar(path, binding=None, *, batch_index=4):
    doc = {"schema": "e422-binding-state-v1", "batch_index": batch_index,
           "binding": binding if binding is not None else {}}
    doc["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in doc.items() if k != "sha256"})
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(doc) + b"\n")
    return doc


@pytest.fixture()
def forbid_driver(monkeypatch):
    """Sentinel: the frozen driver must NEVER be reached on a refusal."""
    def _boom(*_a, **_k):
        raise AssertionError("frozen driver reached despite gate refusal")
    monkeypatch.setattr(driver, "run_scheduled_evaluation", _boom)
    return monkeypatch


# --------------------------------------------------------------------------- #
# Verification function.
# --------------------------------------------------------------------------- #
def test_verified_load_absent_optional(tmp_path):
    assert gate.load_verified_binding_state(
        str(tmp_path / "none.json"), required=False,
        max_batch_index=4) is None


def test_verified_load_absent_required_typed(tmp_path):
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_MISSING"):
        gate.load_verified_binding_state(
            str(tmp_path / "none.json"), required=True, max_batch_index=8)


def test_verified_load_valid(tmp_path):
    path = str(tmp_path / "sidecar.json")
    doc = make_sidecar(path, binding={"p5_90": {"bound_at_batch": 4}})
    assert gate.load_verified_binding_state(
        path, required=True, max_batch_index=8) == doc


def test_verified_load_tampered_binding_typed(tmp_path):
    path = str(tmp_path / "sidecar.json")
    doc = make_sidecar(path)
    doc["binding"]["p5_90"] = {"outcome": "FORGED"}
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(doc) + b"\n")
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_DIGEST_MISMATCH"):
        gate.load_verified_binding_state(path, required=True,
                                         max_batch_index=8)


def test_verified_load_malformed_json_typed(tmp_path):
    path = str(tmp_path / "sidecar.json")
    os.makedirs(tmp_path, exist_ok=True)
    with open(path, "wb") as handle:
        handle.write(b'{"schema": "e422-binding-state-v1", "binding": ')
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_SCHEMA"):
        gate.load_verified_binding_state(path, required=True,
                                         max_batch_index=8)


def test_verified_load_wrong_schema_or_scalar_typed(tmp_path):
    scalar = str(tmp_path / "scalar.json")
    with open(scalar, "wb") as handle:
        handle.write(b"17")
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_SCHEMA"):
        gate.load_verified_binding_state(scalar, required=True,
                                         max_batch_index=8)
    path = str(tmp_path / "wrong.json")
    make_sidecar(path)
    doc = json.loads(open(path, "rb").read().decode())
    doc["schema"] = "something-else-v9"
    doc["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in doc.items() if k != "sha256"})
    with open(path, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(doc) + b"\n")
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_SCHEMA"):
        gate.load_verified_binding_state(path, required=True,
                                         max_batch_index=8)


# --------------------------------------------------------------------------- #
# Gated evaluation end-to-end (offline, real frozen driver).
# --------------------------------------------------------------------------- #
def _build_ledger(tmp_path, batches=(1, 2, 3, 4)):
    for k in batches:
        write_checkpoint(str(tmp_path), k, synth_rows(seed=10 + k))


def test_gated_batch4_creates_sidecar(tmp_path):
    _build_ledger(tmp_path)
    sidecar = str(tmp_path / "binding.json")
    artifact = gate.run_gated_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=sidecar, now_iso="2026-10-08T00:00:00Z")
    assert artifact["role"] == "batch4_day28"
    assert os.path.exists(sidecar)  # created + post-write verified
    state = json.loads(open(sidecar, "rb").read())
    assert state["schema"] == "e422-binding-state-v1"


def test_gated_batch4_refuses_preexisting_tampered_sidecar(tmp_path,
                                                           forbid_driver):
    _build_ledger(tmp_path)
    sidecar = str(tmp_path / "binding.json")
    make_sidecar(sidecar)
    doc = json.loads(open(sidecar, "rb").read())
    doc["binding"]["p5_90"] = {"outcome": "FORGED"}
    with open(sidecar, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(doc) + b"\n")
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_DIGEST_MISMATCH"):
        gate.run_gated_evaluation(
            T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=sidecar, now_iso="2026-10-08T00:00:00Z")


def test_gated_batch8_requires_existing_sidecar(tmp_path, forbid_driver):
    _build_ledger(tmp_path, batches=range(1, 9))
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_MISSING"):
        gate.run_gated_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=str(tmp_path / "absent.json"),
            now_iso="2026-11-05T00:00:00Z")


def test_gated_batch8_refuses_tampered_sidecar(tmp_path, forbid_driver):
    _build_ledger(tmp_path, batches=range(1, 9))
    sidecar = str(tmp_path / "binding.json")
    make_sidecar(sidecar)
    doc = json.loads(open(sidecar, "rb").read())
    doc["binding"]["p6"] = {"outcome": "FORGED", "bound_at_batch": 99}
    with open(sidecar, "wb") as handle:
        handle.write(manifest_mod.canonical_json_bytes(doc) + b"\n")
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_DIGEST_MISMATCH"):
        gate.run_gated_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=sidecar, now_iso="2026-11-05T00:00:00Z")


def test_gated_batch8_happy_chain_and_binding_outcomes(tmp_path):
    # the healthy-ledger chain: batch-4 creates, batch-8 chains; AND the
    # previously unpinned pin — P6/P7a/P7b bind at batch 8 on the
    # decided (six-bin) grid, P5 carried, sidecar re-verified post-write
    _build_ledger(tmp_path, batches=range(1, 9))
    sidecar = str(tmp_path / "binding.json")
    gate.run_gated_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=sidecar, now_iso="2026-10-08T00:00:00Z")
    b4_state = json.loads(open(sidecar, "rb").read())
    artifact = gate.run_gated_evaluation(
        T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=sidecar, now_iso="2026-11-05T00:00:00Z")
    assert artifact["floor_decision"]["decision"] == driver.SIX_BIN
    binding = artifact["binding"]
    for key in ("p6", "p7a", "p7b"):
        entry = binding[key]
        if entry["outcome"] != driver.NOT_YET_EVALUABLE:
            assert entry["bound_at_batch"] == 8
    # prior P5 bindings are immutable across the chain
    for key in ("p5_90", "p5_95"):
        if key in b4_state["binding"]:
            assert binding[key]["outcome"] == \
                b4_state["binding"][key]["outcome"]
            assert binding[key]["bound_at_batch"] == 4


# --------------------------------------------------------------------------- #
# Detached-review refusals (2026-09-10, wf_d8bff144): registered-T0 pin,
# None-path chain break, chain-position hole, unregistered keys,
# double-evaluation guard.
# --------------------------------------------------------------------------- #
DEV_T0 = "2026-09-10T00:00:00Z"  # deliberately NOT the registered bind


def test_gate_pins_registered_t0(tmp_path, forbid_driver):
    _build_ledger(tmp_path)
    with pytest.raises(gate.GateError, match="E424_T0_MISMATCH"):
        gate.run_gated_evaluation(
            DEV_T0, 4, str(tmp_path), str(tmp_path / "results"),
            seeds=range(3), prior_binding_path=str(tmp_path / "b.json"),
            now_iso="2026-10-08T00:00:00Z")


def test_gate_none_path_is_chain_break_at_every_batch(tmp_path,
                                                      forbid_driver):
    # prior_binding_path=None must refuse at EVERY scheduled batch: at
    # batch >= 8 the frozen driver would silently start a fresh study;
    # at batch 4 it would silently never WRITE the sidecar (the driver
    # only writes when a path is supplied).  Batch 4's registered fresh
    # study is an absent file with a supplied path.
    _build_ledger(tmp_path)
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_MISSING"):
        gate.run_gated_evaluation(
            T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=None, now_iso="2026-10-08T00:00:00Z")
    _build_ledger(tmp_path, batches=range(1, 9))
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_MISSING"):
        gate.run_gated_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=None, now_iso="2026-11-05T00:00:00Z")


def test_gate_refuses_self_and_future_chain_position(tmp_path,
                                                     forbid_driver):
    _build_ledger(tmp_path, batches=range(1, 9))
    # the review's live-reproduced hole: a digest-valid sidecar written
    # AT (or after) the evaluated batch is not a registered prior state
    for stale in (4, 8):
        sidecar = str(tmp_path / f"stale{stale}.json")
        make_sidecar(sidecar, batch_index=stale)
        with pytest.raises(gate.GateError,
                           match="E424_BINDING_STATE_CHAIN_POSITION"):
            gate.run_gated_evaluation(
                T0, 4, str(tmp_path), str(tmp_path / "results"),
                seeds=range(3), prior_binding_path=sidecar,
                now_iso="2026-10-08T00:00:00Z")
    # non-integer and zero positions refuse typed as well
    for bad in ("4", 0):
        sidecar = str(tmp_path / f"pos{bad}.json")
        make_sidecar(sidecar, batch_index=bad)
        with pytest.raises(gate.GateError,
                           match="E424_BINDING_STATE_CHAIN_POSITION"):
            gate.load_verified_binding_state(sidecar, required=True,
                                             max_batch_index=4)


def test_gate_refuses_unregistered_binding_keys(tmp_path, forbid_driver):
    _build_ledger(tmp_path, batches=range(1, 9))
    sidecar = str(tmp_path / "sidecar.json")
    make_sidecar(sidecar, batch_index=4,
                 binding={"p5_90": {"bound_at_batch": 4},
                          "rogue_prediction": {"outcome": "BOUND"}})
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_CHAIN_POSITION"):
        gate.run_gated_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=sidecar, now_iso="2026-11-05T00:00:00Z")


def test_gate_non_dict_binding_with_valid_digest_typed(tmp_path):
    # digest-valid but structurally wrong (binding is a scalar): the
    # schema check must fire AFTER the digest verifies, not before
    path = str(tmp_path / "scalar_binding.json")
    make_sidecar(path, binding=5)
    with pytest.raises(gate.GateError, match="E424_BINDING_STATE_SCHEMA"):
        gate.load_verified_binding_state(path, required=True,
                                         max_batch_index=8)


def test_gate_refuses_rerun_when_outputs_exist(tmp_path, forbid_driver):
    # any anchored output (evaluation artifact, batch-8 floor receipt,
    # e425 execution receipt) existing means the evaluation already ran
    _build_ledger(tmp_path, batches=(1, 2, 3, 4))
    results = str(tmp_path / "results")
    os.makedirs(results, exist_ok=True)
    variants = [
        "e422_batch04_evaluation.json",
        "e422_batch04_execution_receipt.json",
    ]
    for name in variants:
        with open(os.path.join(results, name), "wb") as handle:
            handle.write(b"{}\n")
        with pytest.raises(gate.GateError, match="E424_EVALUATION_EXISTS"):
            gate.run_gated_evaluation(
                T0, 4, str(tmp_path), results, seeds=range(3),
                prior_binding_path=str(tmp_path / "b.json"),
                now_iso="2026-10-08T00:00:00Z")
        os.remove(os.path.join(results, name))
    # batch 8: the floor-decision receipt also anchors
    _build_ledger(tmp_path, batches=range(1, 9))
    with open(os.path.join(results,
                           "e422_batch08_floor_decision.json"), "wb") as h:
        h.write(b"{}\n")
    with pytest.raises(gate.GateError, match="E424_EVALUATION_EXISTS"):
        gate.run_gated_evaluation(
            T0, 8, str(tmp_path), results, seeds=range(3),
            prior_binding_path=str(tmp_path / "b.json"),
            now_iso="2026-11-05T00:00:00Z")


def test_gate_exists_guard_precedes_sidecar_checks(tmp_path, forbid_driver):
    # guard order: an existing output refuses even when the sidecar
    # argument would refuse first — the re-run is the primary signal
    _build_ledger(tmp_path, batches=(1, 2, 3, 4))
    results = str(tmp_path / "results")
    os.makedirs(results, exist_ok=True)
    with open(os.path.join(results,
                           "e422_batch04_evaluation.json"), "wb") as h:
        h.write(b"{}\n")
    with pytest.raises(gate.GateError, match="E424_EVALUATION_EXISTS"):
        gate.run_gated_evaluation(
            T0, 4, str(tmp_path), results, seeds=range(3),
            prior_binding_path=str(tmp_path / "absent.json"),
            now_iso="2026-10-08T00:00:00Z")


def test_gate_batch4_absent_file_supplied_path_creates_sidecar(tmp_path):
    # batch 4's registered fresh study: the sidecar FILE is absent but
    # the PATH is supplied — the driver creates it, post-write verified
    _build_ledger(tmp_path)
    sidecar = str(tmp_path / "binding.json")
    artifact = gate.run_gated_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=sidecar, now_iso="2026-10-08T00:00:00Z")
    assert artifact["role"] == "batch4_day28"
    assert os.path.exists(sidecar)


def test_gate_post_write_failure_is_typed(tmp_path, monkeypatch):
    # negative test for the write-path check (detached review MINOR):
    # if the (re)written sidecar fails re-verification, the gate raises
    # WRITE_UNVERIFIED (wrapping the underlying code) — the evaluation
    # artifact exists but is NOT silently trusted as chained
    _build_ledger(tmp_path)
    sidecar = str(tmp_path / "binding.json")
    real_load = gate.load_verified_binding_state
    calls = []

    def _flaky(path, *, required, max_batch_index):
        calls.append(max_batch_index)
        if len(calls) == 1:
            return real_load(path, required=required,
                             max_batch_index=max_batch_index)  # pre-verify OK
        raise gate.GateError(
            gate.GateFailureCode.BINDING_STATE_DIGEST_MISMATCH.value,
            "synthetic post-write tamper")

    monkeypatch.setattr(gate, "load_verified_binding_state", _flaky)
    with pytest.raises(gate.GateError,
                       match="E424_BINDING_STATE_WRITE_UNVERIFIED"):
        gate.run_gated_evaluation(
            T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=sidecar, now_iso="2026-10-08T00:00:00Z")
    assert len(calls) == 2  # pre-driver + post-write both ran
