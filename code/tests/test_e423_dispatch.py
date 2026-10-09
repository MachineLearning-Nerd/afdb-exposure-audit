"""Offline tests for the e423 dispatch CLI (no network)."""

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e423_dispatch as dispatch  # noqa: E402
import e423_runner as runner  # noqa: E402
from e423_runner import batch_window_utc  # noqa: E402

T0 = "2026-09-11T00:00:00Z"


def test_default_t0_is_frozen_bind():
    # the CLI default must equal the T0 bound in the frozen registration
    import re
    text = open(os.path.join(os.path.dirname(__file__), "..",
                             "docs", "e422_registration.md")).read()
    match = re.search(r"T0\s*=\s*(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)", text)
    assert match, "frozen text carries no T0 bind"
    assert dispatch.DEFAULT_T0 == match.group(1)


def test_dry_run_prints_registered_window_and_query(capsys):
    rc = dispatch.main(["--t0", T0, "--batch", "3", "--raw-dir", "r",
                        "--out-dir", "o", "--dry-run"])
    out = capsys.readouterr().out
    assert rc == 0
    start, close = batch_window_utc(T0, 3)
    assert start in out and close in out
    assert "X-RAY DIFFRACTION" in out
    assert '"rows": 500' in out
    assert "dry-run" in out


def test_batch_window_arithmetic():
    # registered T0 under the 48-h batch-1 acceleration amendment
    # (DECISIONS 2026-09-10): batch 1 = (T0, T0+2d]; batches 2-12 stay
    # 7-day and shift by batch 1's early close (start(k) = close(k-1))
    assert batch_window_utc(T0, 1) == ("2026-09-11T00:00:00Z",
                                       "2026-09-13T00:00:00Z")
    assert batch_window_utc(T0, 4) == ("2026-09-27T00:00:00Z",
                                       "2026-10-04T00:00:00Z")


def test_amendment_registered_schedule_closes():
    # close(k) = T0 + 2 + 7(k-1) days for k >= 2: the nominal day-28 /
    # 56 / 84 evaluation labels land at ACTUAL days 23 / 51 / 79
    # (disclosed in the runbook + DECISIONS; driver roles are keyed to
    # batch index, not dates)
    closes = {k: batch_window_utc(T0, k)[1] for k in (2, 4, 8, 12)}
    assert closes == {
        2: "2026-09-20T00:00:00Z",
        4: "2026-10-04T00:00:00Z",   # "day-28" role, actual day 23
        8: "2026-11-01T00:00:00Z",   # "day-56" role, actual day 51
        12: "2026-11-29T00:00:00Z",  # "day-84" role, actual day 79
    }


def test_amendment_windows_are_contiguous():
    # the amended schedule stays gap-free and overlap-free: the batch-1
    # 48-h window abuts batch 2 exactly
    for k in range(1, 12):
        assert batch_window_utc(T0, k + 1)[0] == batch_window_utc(T0, k)[1]


def test_amendment_is_registered_t0_only():
    # any NON-registered t0 (dev shakedowns, historical rehearsals)
    # keeps the pure 7-day arithmetic — the amendment touches only the
    # registered cadence
    dev = "2026-09-03T00:00:00Z"
    assert batch_window_utc(dev, 1) == ("2026-09-03T00:00:00Z",
                                        "2026-09-10T00:00:00Z")
    assert batch_window_utc(dev, 4) == ("2026-09-24T00:00:00Z",
                                        "2026-10-01T00:00:00Z")


def test_registered_t0_constant_matches_dispatch_pin():
    # single source of truth for the registered bind across the runner
    # (window arithmetic) and the dispatch CLI (T0 pin)
    assert runner.REGISTERED_T0 == dispatch.DEFAULT_T0


def test_post_transport_exists_and_wraps_requests(monkeypatch):
    import e420_census as census
    captured = {}

    class FakeResponse:
        status_code = 200
        headers = {"Content-Type": "application/json"}
        content = b"{}"

    def fake_post(url, data=None, headers=None, timeout=None):
        captured.update(url=url, data=data, headers=headers,
                        timeout=timeout)
        return FakeResponse()

    import requests
    monkeypatch.setattr(requests, "post", fake_post)
    transport = dispatch.PostTransport(user_agent="ua-test")
    resp = transport.post_json("https://search.rcsb.org/rcsbsearch/v2/query",
                               b'{"query":1}', 30.0)
    assert isinstance(resp, census.HttpResponse)
    assert resp.status == 200
    assert captured["data"] == b'{"query":1}'
    assert captured["headers"]["Content-Type"] == "application/json"
    assert captured["timeout"] == 30.0


def test_dispatch_pins_registered_t0(monkeypatch, tmp_path, capsys):
    # REWRITTEN 2026-09-10 (w1:pG, detached review wf_d8bff144 MAJOR):
    # --t0 is procedurally pinned now — the registered dispatch entry
    # point refuses ANY t0 other than the frozen gate-(a) bind with
    # typed E423_T0_MISMATCH, exit 2, BEFORE window arithmetic or
    # network (any t0 yields self-consistently "closed" windows, so a
    # dev convenience t0 on the registered path was a real hole).  The
    # census-failure path this test used to exercise over a historical
    # window is covered runner-side
    # (tests/test_e423_runner.py::test_census_400_is_typed_not_silent_empty).
    import runpy
    import types

    class FakeResponse:  # must never be reached
        status_code = 200
        headers = {}
        content = b"{}"

    fake_requests = types.ModuleType("requests")
    fake_requests.post = lambda *a, **k: FakeResponse()
    monkeypatch.setitem(sys.modules, "requests", fake_requests)
    monkeypatch.setattr(sys, "argv", [
        "e423_dispatch.py", "--t0", "2026-08-21T00:00:00Z", "--batch", "1",
        "--raw-dir", str(tmp_path / "raw"), "--out-dir", str(tmp_path / "out")])
    script = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "experiments", "e423_dispatch.py"))
    with pytest.raises(SystemExit) as excinfo:
        runpy.run_path(script, run_name="__main__")
    assert excinfo.value.code == 2
    err = capsys.readouterr().err
    assert "TYPED_FAILURE E423_T0_MISMATCH" in err
    assert "2026-08-21T00:00:00Z" in err
    assert not (tmp_path / "raw").exists()   # no side effects at all
    assert not (tmp_path / "out").exists()


def test_typed_manifest_error_exits_nonzero(monkeypatch, tmp_path, capsys):
    # widened envelope 2026-09-10: the runner's pre-network manifest
    # verification (detached review) raises ManifestError on a tampered
    # checkpoint — the dispatch envelope must type it, not the residual
    # clause.  (The pre-network verification itself is proven runner-side
    # in tests/test_e423_runner.py; here the registered-T0 batch-1 window
    # is still unclosed in real time, so run_batch is stubbed with the
    # exact ManifestError the tamper path raises — the envelope is the
    # thing under test.)
    import types

    import e422_manifest as manifest_mod

    def _tampered(*_a, **_k):
        raise manifest_mod.ManifestError(
            manifest_mod.ManifestFailureCode.MANIFEST_DIGEST_MISMATCH,
            "synthetic tampered manifest checkpoint")

    class FakeResponse:
        status_code = 200
        headers = {}
        content = b"{}"

    fake_requests = types.ModuleType("requests")
    fake_requests.post = lambda *a, **k: FakeResponse()
    monkeypatch.setitem(sys.modules, "requests", fake_requests)
    monkeypatch.setattr(dispatch, "run_batch", _tampered)
    monkeypatch.setattr(sys, "argv", [
        "e423_dispatch.py", "--t0", T0, "--batch", "1",
        "--raw-dir", str(tmp_path / "raw"), "--out-dir", str(tmp_path / "out")])
    assert dispatch.main() == 2
    err = capsys.readouterr().err
    assert "TYPED_FAILURE E422_MANIFEST_DIGEST_MISMATCH" in err
    assert "E423_UNEXPECTED" not in err


def test_main_guard_widened_types_unexpected(monkeypatch, tmp_path, capsys):
    # e423 (detached review MINOR-5): the __main__ guard's second clause
    # types ANY residual non-CensusError exception — class name on stderr,
    # exit 2, never a bare traceback.
    import runpy
    import types

    class FakeResponse:
        status_code = 200
        headers = {}
        content = b"{}"

    fake_requests = types.ModuleType("requests")
    fake_requests.post = lambda *a, **k: FakeResponse()

    def boom(*a, **k):
        raise RuntimeError("synthetic residual failure")

    monkeypatch.setitem(sys.modules, "requests", fake_requests)
    import e423_runner as runner_mod
    monkeypatch.setattr(runner_mod, "run_batch", boom)
    monkeypatch.setattr(sys, "argv", [
        "e423_dispatch.py", "--t0", T0, "--batch", "1",
        "--raw-dir", str(tmp_path / "raw"), "--out-dir", str(tmp_path / "out")])
    script = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "experiments", "e423_dispatch.py"))
    with pytest.raises(SystemExit) as excinfo:
        runpy.run_path(script, run_name="__main__")
    assert excinfo.value.code == 2
    assert "TYPED_FAILURE E423_UNEXPECTED_RuntimeError" in capsys.readouterr().err
