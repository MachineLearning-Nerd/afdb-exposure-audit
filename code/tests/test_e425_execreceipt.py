"""Offline tests for the e425 execution-receipt layer (2026-09-10).

Layering under test: run_logged_evaluation = timing + receipt + e424
gate + FROZEN driver.  The receipt is a per-run measurement (real
clock); the evaluation artifacts stay deterministic.
"""

from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_evaldriver as driver  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e425_execreceipt as receipt_mod  # noqa: E402

from tests.test_e424_evalgate import (make_sidecar, synth_rows,  # noqa: E402
                                      write_checkpoint)

# the e424 gate pins t0 == the frozen registration bind (detached review
# 2026-09-10), so the receipt-layer tests run at the registered T0
T0 = "2026-09-11T00:00:00Z"


def _build_ledger(tmp_path, batches=(1, 2, 3, 4)):
    for k in batches:
        write_checkpoint(str(tmp_path), k, synth_rows(seed=10 + k))


def _receipt_path(tmp_path):
    return os.path.join(str(tmp_path / "results"),
                        "e422_batch04_execution_receipt.json")


def test_receipt_written_and_digest_verifies(tmp_path):
    _build_ledger(tmp_path)
    artifact = receipt_mod.run_logged_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-10-08T00:00:00Z")
    assert artifact["role"] == "batch4_day28"  # artifact returned unchanged
    path = _receipt_path(tmp_path)
    assert os.path.exists(path)
    receipt = json.loads(open(path, "rb").read())
    assert receipt["schema"] == receipt_mod.RECEIPT_SCHEMA
    assert receipt["role"] == "batch4_day28"
    # evaluated_utc is copied from the artifact (detached review 2026-09-10:
    # the receipt names the evaluation instant it measured around)
    assert receipt["evaluated_utc"] == artifact["evaluated_utc"]
    assert receipt["evaluated_utc"] == "2026-10-08T00:00:00Z"
    assert receipt["registered_budget_hours"] == 24.0
    assert receipt["within_budget"] is True
    assert receipt["deviation_note"] is None
    assert receipt["wall_seconds"] >= 0.0
    assert receipt["sha256"] == manifest_mod.canonical_sha256(
        {k: v for k, v in receipt.items() if k != "sha256"})


def test_receipt_reports_deviation_when_budget_exceeded(tmp_path,
                                                        monkeypatch):
    monkeypatch.setitem(receipt_mod.BUDGET_HOURS, 4, -1.0)  # force exceed
    _build_ledger(tmp_path)
    receipt_mod.run_logged_evaluation(
        T0, 4, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
        prior_binding_path=str(tmp_path / "binding.json"),
        now_iso="2026-10-08T00:00:00Z")
    receipt = json.loads(open(_receipt_path(tmp_path), "rb").read())
    assert receipt["within_budget"] is False
    note = receipt["deviation_note"]
    assert note is not None
    assert "budget exceeded" in note
    assert "never truncated" in note
    assert "no repeat is dropped" in note


def test_no_receipt_when_gate_refuses(tmp_path, monkeypatch):
    # a refused run produced NO execution — nothing to report; the
    # receipt must not exist (and the frozen driver is never reached)
    def _boom(*_a, **_k):
        raise AssertionError("driver reached despite gate refusal")
    monkeypatch.setattr(driver, "run_scheduled_evaluation", _boom)
    _build_ledger(tmp_path, batches=range(1, 9))
    sidecar = str(tmp_path / "binding.json")
    make_sidecar(sidecar)  # valid, but batch-8 chain requires prior state:
    # make the gate refuse via a MISSING sidecar instead
    missing = str(tmp_path / "absent.json")
    with pytest.raises(Exception, match="E424_BINDING_STATE_MISSING"):
        receipt_mod.run_logged_evaluation(
            T0, 8, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=missing, now_iso="2026-11-05T00:00:00Z")
    assert not os.path.exists(os.path.join(
        str(tmp_path / "results"), "e422_batch08_execution_receipt.json"))


def test_budgets_match_registration(tmp_path):
    # registered section 6 wall-clock budgets (frozen text dc7a3627,
    # docs/e422_registration.md lines 548, 566-570: <= 24 h at batches
    # 4 and 8; <= 72 h at batch 12)
    assert receipt_mod.BUDGET_HOURS == {4: 24.0, 8: 24.0, 12: 72.0}
    assert set(receipt_mod.BUDGET_HOURS) == set(driver.SCHEDULE)


def test_unscheduled_batch_is_typed_not_keyerror(tmp_path):
    # an operator typo (batch outside SCHEDULE) must surface as the
    # frozen driver's typed refusal, not a receipt-layer KeyError —
    # the budget lookup happens only after the gated call (2026-09-10)
    _build_ledger(tmp_path, batches=(1,))
    with pytest.raises(driver.DriverError, match="E422_UNSCHEDULED_EVALUATION"):
        receipt_mod.run_logged_evaluation(
            T0, 5, str(tmp_path), str(tmp_path / "results"), seeds=range(3),
            prior_binding_path=str(tmp_path / "absent.json"),
            now_iso="2026-10-08T00:00:00Z")
    assert not os.path.exists(os.path.join(
        str(tmp_path / "results"), "e422_batch05_execution_receipt.json"))


# --------------------------------------------------------------------------- #
# Registered CLI entry point + typed-failure envelope (detached review
# 2026-09-10): the runbook's scheduled-evaluation command is this main(),
# and a gate/driver/census refusal exits 2 with TYPED_FAILURE <code> —
# never a bare traceback, never a silent drop.
# --------------------------------------------------------------------------- #
def test_main_passes_registered_arguments(tmp_path, monkeypatch, capsys):
    seen = {}

    def _recorder(t0, batch, out_dir, result_dir, *, prior_binding_path,
                  seeds=None, now_iso=None):
        seen.update(t0=t0, batch=batch, out_dir=out_dir,
                    result_dir=result_dir, prior_binding_path=prior_binding_path)
        return {"role": "stub"}

    monkeypatch.setattr(receipt_mod, "run_logged_evaluation", _recorder)
    code = receipt_mod.main([
        "--batch", "4",
        "--out-dir", str(tmp_path / "batches"),
        "--result-dir", str(tmp_path / "evals"),
        "--prior-binding-path", str(tmp_path / "binding.json"),
    ])
    assert code == 0
    assert seen == {"t0": receipt_mod.REGISTERED_T0, "batch": 4,
                    "out_dir": str(tmp_path / "batches"),
                    "result_dir": str(tmp_path / "evals"),
                    "prior_binding_path": str(tmp_path / "binding.json")}
    out = capsys.readouterr()
    assert out.err == ""


def test_main_typed_envelope_on_gate_refusal(tmp_path, monkeypatch, capsys):
    # a batch-8 chain break (missing sidecar) through the CLI envelope
    _build_ledger(tmp_path, batches=range(1, 9))
    argv = ["--batch", "8",
            "--out-dir", str(tmp_path),
            "--result-dir", str(tmp_path / "results"),
            "--prior-binding-path", str(tmp_path / "absent.json")]
    assert receipt_mod.main(argv) == 2
    err = capsys.readouterr().err
    assert err.startswith("TYPED_FAILURE E424_BINDING_STATE_MISSING")
    # and no execution receipt was written for the refused run
    assert not os.path.exists(os.path.join(
        str(tmp_path / "results"),
        "e422_batch08_execution_receipt.json"))


def test_main_typed_envelope_on_driver_and_census_refusal(tmp_path,
                                                          monkeypatch,
                                                          capsys):
    monkeypatch.setattr(
        receipt_mod, "run_logged_evaluation",
        lambda *a, **k: (_ for _ in ()).throw(
            driver.DriverError("E422_CHECKPOINT_MISSING", "gone.json")))
    assert receipt_mod.main(["--batch", "4"]) == 2
    assert capsys.readouterr().err.startswith(
        "TYPED_FAILURE E422_CHECKPOINT_MISSING")

    monkeypatch.setattr(
        receipt_mod, "run_logged_evaluation",
        lambda *a, **k: (_ for _ in ()).throw(
            receipt_mod.census.CensusError("E420_INTERNAL_ERROR", "x")))
    assert receipt_mod.main(["--batch", "4"]) == 2
    assert capsys.readouterr().err.startswith(
        "TYPED_FAILURE E420_INTERNAL_ERROR")
