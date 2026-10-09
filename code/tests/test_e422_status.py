"""Tests for the e422 finalized analysis-status artifact (F-10)."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_status as st  # noqa: E402


def test_every_status_has_a_recording_location():
    for status in st.E422Status:
        assert status.name in st.RECORDING_LOCATION, status.name


def test_no_orphan_locations():
    names = {s.name for s in st.E422Status}
    assert set(st.RECORDING_LOCATION) == names


def test_registration_statuses_present():
    # e422 §3/§6 analysis statuses
    for name in ("FIT_FAILURE", "UNDER_SUPPORTED_STRATUM",
                 "DEFERRED_INSUFFICIENT_SUPPORT",
                 "DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE"):
        assert hasattr(st.E422Status, name)
    # e422 §2 manifest/source failures
    for name in ("PRE_SNAPSHOT_RELEASE", "MODEL_CREATED_DATE_MISSING",
                 "MODEL_CREATED_DATE_MISMATCH", "AFDB_VERSION_MISMATCH",
                 "CHILD_URL_NOT_V6", "ABSENT_FROM_SNAPSHOT",
                 "SOURCE_CONTRACT_FAILURE"):
        assert hasattr(st.E422Status, name)
    # e422 §7 prediction outcomes
    for name in ("HOLD", "FAIL", "REPLICATE", "FALSIFY",
                 "WIN", "MIRRORED_WIN", "NO_DECISION"):
        assert hasattr(st.E422Status, name)


def test_manifest_codes_agree_with_e422_manifest_module():
    import e422_manifest as em

    manifest_values = {c.value for c in em.ManifestFailureCode}
    registry_values = {s.value for s in st.E422Status}
    assert manifest_values <= registry_values


def test_registry_digest_is_stable():
    assert st.registry_digest() == st.registry_digest()
    assert len(st.registry_digest()) == 64
