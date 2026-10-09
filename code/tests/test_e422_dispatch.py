"""Offline tests for the e422 dispatch CLI (no network).

AMENDED 2026-09-10 under e423 (DECISIONS 2026-09-10): e422_dispatch.py is
now a refuse-to-run tombstone (detached review MAJOR-1 / F14 — its frozen
wiring still carries the invalid census query and the silent-ok POST lane,
one command away from reproducing the empty-checkpoint defect), so the CLI
tests assert the TOMBSTONE REFUSES.  The superseded DEFAULT_T0 /
PostTransport tests were removed with the dead wiring; the batch-window
test is retained against the FROZEN e422_runner module, whose bytes are
unchanged.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_dispatch as dispatch  # noqa: E402
from e422_runner import batch_window_utc  # noqa: E402

T0 = "2026-09-11T00:00:00Z"


def test_tombstone_refuses_unconditionally(capsys):
    # AMENDED 2026-09-10 under e423: main() is a refuse-to-run tombstone
    # (typed E423_DISPATCH_SUPERSEDED, exit 2); even --dry-run refuses,
    # because the query it prints is the invalid one
    rc = dispatch.main(["--t0", T0, "--batch", "3", "--raw-dir", "r",
                        "--out-dir", "o", "--dry-run"])
    err = capsys.readouterr().err
    assert rc == 2
    assert "E423_DISPATCH_SUPERSEDED" in err
    assert "e423_dispatch.py" in err


def test_tombstone_refuses_without_any_args(capsys):
    assert dispatch.main([]) == 2
    assert "SUPERSEDED" in capsys.readouterr().err


def test_batch_window_arithmetic():
    # retained: the FROZEN e422_runner window arithmetic is unchanged
    assert batch_window_utc(T0, 1) == ("2026-09-11T00:00:00Z",
                                       "2026-09-18T00:00:00Z")
    assert batch_window_utc(T0, 4) == ("2026-10-02T00:00:00Z",
                                       "2026-10-09T00:00:00Z")
