"""Offline tests for the e426 attested-empty dispatch (no network)."""

from __future__ import annotations

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import e426_attested_empty as e426  # noqa: E402
import e423_runner as runner  # noqa: E402
import e420_census as census  # noqa: E402


def _census_body(total: int) -> bytes:
    return json.dumps({"total_count": total, "result_set": []},
                      sort_keys=True).encode("utf-8")


def _fake_poster(counts: dict[str, int]) -> object:
    """Monkeypatch e426._post_query with a window-keyed fake.

    keys: 'hist' (historical week), 'event' (09-16 window),
    'b1' (batch-1 window).  An HTTP-204 window is simulated with the
    sentinel string '204'."""
    def post(query: dict) -> tuple[int, bytes]:
        start = query["query"]["nodes"][2]["parameters"]["value"]
        close = query["query"]["nodes"][3]["parameters"]["value"]
        if (start, close) == e426.HISTORICAL_WINDOW:
            key = "hist"
        elif (start, close) == e426.EVENT_WINDOW:
            key = "event"
        else:
            key = "b1"
        val = counts[key]
        if val == "204":
            return 204, b""
        return 200, _census_body(val)
    return post


def test_attested_empty_body_is_canonical_and_valid():
    doc = json.loads(e426.ATTESTED_EMPTY_BODY.decode("utf-8"))
    assert doc == {"result_set": [], "total_count": 0}
    assert e426.ATTESTED_EMPTY_BODY == json.dumps(
        doc, sort_keys=True).encode("utf-8")


def test_attestation_passes_on_healthy_recovered_lane(monkeypatch):
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 140, "event": 156, "b1": "204"}))
    probes = e426.attest_empty_batch1()
    assert [p["total_count"] for p in probes] == [140, 156, None]
    assert probes[2]["http_status"] == 204


def test_attestation_a1_lane_drift_fails(monkeypatch):
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 139, "event": 156, "b1": "204"}))
    with pytest.raises(e426.E426Error) as ei:
        e426.attest_empty_batch1()
    assert ei.value.code == "E426_LANE_HEALTH_FAILED"


def test_attestation_a2_still_frozen_fails(monkeypatch):
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 140, "event": "204", "b1": "204"}))
    with pytest.raises(e426.E426Error) as ei:
        e426.attest_empty_batch1()
    assert ei.value.code == "E426_NO_POST_FREEZE_RELEASES"


def test_attestation_a3_window_not_empty_fails(monkeypatch):
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 140, "event": 156, "b1": 7}))
    with pytest.raises(e426.E426Error) as ei:
        e426.attest_empty_batch1()
    assert ei.value.code == "E426_WINDOW_NOT_EMPTY"


def test_transport_get_refuses_and_post_is_stable(tmp_path):
    t = e426.AttestedEmptyTransport(e426.ATTESTED_EMPTY_BODY, "deadbeef")
    resp = t.post_json("https://search.rcsb.org/rcsbsearch/v2/query",
                       b"{}", 30.0)
    assert resp.status == 200
    assert resp.content == e426.ATTESTED_EMPTY_BODY
    assert resp.headers["X-E426-Attested-Empty"] == "deadbeef"
    with pytest.raises(e426.E426Error) as ei:
        t.get("https://www.ebi.ac.uk/x", 30.0)
    assert ei.value.code == "E426_TRANSPORT_GET_UNEXPECTED"


def test_zero_row_checkpoint_chains_through_frozen_path(tmp_path, monkeypatch):
    """The full dispatch path: attested-empty transport through the
    FROZEN run_batch produces a quiet-window checkpoint that the frozen
    eval-side ledger loader accepts digest-clean."""
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 140, "event": 156, "b1": "204"}))
    probes = e426.attest_empty_batch1()
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "batches"
    manifest_path = tmp_path / "manifest_checkpoint.json"
    transport = e426.AttestedEmptyTransport(
        e426.ATTESTED_EMPTY_BODY,
        e426.census._sha256_hex(e426.ATTESTED_EMPTY_BODY))
    ckpt = runner.run_batch(
        e426.REGISTERED_T0, 1, transport=transport,
        raw_dir=str(raw_dir), out_dir=str(out_dir),
        manifest_path=str(manifest_path))
    assert transport.post_calls == 1
    assert ckpt["census_entries"] == 0
    assert ckpt["n_rows"] == 0
    assert ckpt["max_entries"] is None
    assert ckpt["window"] == {"start": "2026-09-11T00:00:00Z",
                              "close": "2026-09-13T00:00:00Z"}
    assert ckpt["exclusions"] == {"mapping_failures": [],
                                  "manifest_failures": [],
                                  "triple_failures": []}
    # frozen eval-side loader accepts the zero-row chain digest-clean
    from e422_evaldriver import load_cumulative_ledger
    rows, ledger = load_cumulative_ledger(str(out_dir), 1)
    assert len(ledger) == 1
    assert ledger[0]["sha256"] == ckpt["sha256"]
    assert ledger[0]["n_rows"] == 0
    assert rows == []
    assert manifest_path.exists()


def test_main_execute_writes_receipt_and_checkpoint(tmp_path, monkeypatch,
                                                    capsys):
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "batches"
    manifest_path = tmp_path / "manifest_checkpoint.json"
    attest_path = tmp_path / "e426_batch1_empty_attestation.json"
    monkeypatch.setattr(e426, "RAW_DIR", str(raw_dir))
    monkeypatch.setattr(e426, "OUT_DIR", str(out_dir))
    monkeypatch.setattr(e426, "MANIFEST_PATH", str(manifest_path))
    monkeypatch.setattr(e426, "ATTEST_PATH", str(attest_path))
    monkeypatch.setattr(e426, "_post_query",
                        _fake_poster({"hist": 140, "event": 156, "b1": "204"}))
    rc = e426.main(["--execute"])
    assert rc == 0
    receipt = json.loads(attest_path.read_text())
    assert receipt["schema"] == "e426-empty-attestation-v1"
    assert len(receipt["probes"]) == 3
    assert receipt["synthesized_census_body_utf8"] == '{"result_set": [], "total_count": 0}'
    assert receipt["frozen_bytes_changed"] is False
    ckpt = json.loads((out_dir / "e422_batch01_checkpoint.json").read_text())
    assert ckpt["n_rows"] == 0
    out = capsys.readouterr().out
    assert "attestation PASSED" in out


def test_main_refuses_existing_checkpoint(tmp_path, monkeypatch):
    out_dir = tmp_path / "batches"
    out_dir.mkdir()
    (out_dir / "e422_batch01_checkpoint.json").write_text("{}")
    monkeypatch.setattr(e426, "OUT_DIR", str(out_dir))
    with pytest.raises(census.CensusError) as ei:
        e426.main(["--execute"])
    assert ei.value.code == "E426_CHECKPOINT_EXISTS"


def test_plan_mode_touches_nothing(tmp_path, monkeypatch, capsys):
    raw_dir = tmp_path / "raw"
    out_dir = tmp_path / "batches"
    monkeypatch.setattr(e426, "OUT_DIR", str(out_dir))
    monkeypatch.setattr(e426, "RAW_DIR", str(raw_dir))
    rc = e426.main([])
    assert rc == 0
    assert not raw_dir.exists() and not out_dir.exists()
    assert "plan" in capsys.readouterr().out


def test_cli_typed_failure_exit_2_on_refusal(tmp_path, monkeypatch):
    """The __main__ envelope: a typed refusal exits 2 with
    TYPED_FAILURE on stderr (runpy, mirroring the e423 dispatch tests).

    runpy re-executes the module FRESH from disk, so monkeypatched
    module constants do NOT apply — the module copy is placed in a
    scratch directory instead, which re-derives its registered paths
    from its own __file__ (REPO = the scratch dir), keeping the real
    results/e422 paths untouched.  This is the exact defect that made
    the first draft of this test produce a real (unanchored) dispatch
    against the registered paths — see DECISIONS 2026-09-18."""
    import io
    import contextlib
    import shutil
    import runpy
    scratch = tmp_path / "scratch"
    exp_dir = scratch / "experiments"
    exp_dir.mkdir(parents=True)
    module_src = os.path.join(os.path.dirname(e426.__file__),
                              "e426_attested_empty.py")
    module_copy = exp_dir / "e426_attested_empty.py"
    shutil.copyfile(module_src, module_copy)
    out_dir = scratch / "results" / "e422" / "batches"
    out_dir.mkdir(parents=True)
    (out_dir / "e422_batch01_checkpoint.json").write_text("{}")
    monkeypatch.setattr("sys.argv", [str(module_copy), "--execute"])
    stderr = io.StringIO()
    with pytest.raises(SystemExit) as ei:
        with contextlib.redirect_stderr(stderr):
            runpy.run_path(str(module_copy), run_name="__main__")
    assert ei.value.code == 2
    assert "TYPED_FAILURE E426_CHECKPOINT_EXISTS" in stderr.getvalue()
    # the copy derives its registered paths from the scratch dir: nothing
    # was written next to the real module
    assert not os.path.exists(os.path.join(
        os.path.dirname(e426.__file__), os.pardir, "results", "e422",
        "e426_batch01_probe_marker"))
