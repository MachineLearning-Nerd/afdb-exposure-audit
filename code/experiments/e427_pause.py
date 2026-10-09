"""e427 program-pause guard for the e422 dispatch and evaluation CLIs.

While ``results/e422/PAUSE.json`` exists, the registered dispatch
(experiments/e423_dispatch.py) and scheduled-evaluation
(experiments/e425_execreceipt.py) entry points refuse with the typed failure
E427_PROGRAM_PAUSED (exit 2) BEFORE any network, filesystem, or ledger side
effect.  Scheduled cron prompts run these commands verbatim, so a rebuilt or
forgotten cron cannot dispatch while the pause holds.  Lifting the pause is a
human decision recorded in docs/DECISIONS.md; it is done by removing the file.
"""
from __future__ import annotations

import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAUSE_PATH = os.path.join(ROOT, "results", "e422", "PAUSE.json")
CODE = "E427_PROGRAM_PAUSED"


def pause_reason(path: str | None = None) -> str | None:
    """Return a one-line reason when paused, else None.  An unreadable or
    malformed pause file still pauses (fail-closed)."""
    path = path or PAUSE_PATH  # module global, so tests can monkeypatch it
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            doc = json.load(fh)
        return f"{doc.get('reason', 'paused')} (since {doc.get('paused_utc', '?')})"
    except Exception as exc:  # noqa: BLE001 - fail closed
        return f"pause file present but unreadable ({type(exc).__name__})"
