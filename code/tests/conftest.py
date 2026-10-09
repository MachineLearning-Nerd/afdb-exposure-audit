"""Shared test fixtures."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))


@pytest.fixture(autouse=True)
def _isolate_e422_pause(tmp_path, monkeypatch):
    """Tests never see the live results/e422/PAUSE.json; pause tests point
    e427_pause.PAUSE_PATH at their own file."""
    import e427_pause
    monkeypatch.setattr(e427_pause, "PAUSE_PATH", str(tmp_path / "no_pause.json"))
