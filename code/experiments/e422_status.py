"""e422 finalized analysis-status artifact (enumerate-before-anchor).

The e421-P §2.1 discipline, carried by e422 §3A: every typed status the
registration can emit is enumerated here BEFORE anchor, with its
recording location.  This module is the registry of record; the gate-(a)
freeze records its SHA-256.  Adding a status post-freeze is an amendment
(via DECISIONS), never a silent extension.
"""

from __future__ import annotations

from enum import Enum


class E422Status(str, Enum):
    """Every typed status/failure class in the e422 registration."""

    # -- analysis statuses (e422 §3/§6; recorded in evaluation artifacts) ----
    FIT_FAILURE = "FIT_FAILURE"  # §3A.3 - dropped repeat, R_eff reported
    UNDER_SUPPORTED_STRATUM = "UNDER_SUPPORTED_STRATUM"  # §3B joint cell
    DEFERRED_INSUFFICIENT_SUPPORT = "DEFERRED_INSUFFICIENT_SUPPORT"  # §3/§6
    DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE = (
        "DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE"  # §6 close; F-11: also analysis C
    )

    # -- source/manifest failures (e422 §2 amendment; manifest + runner) -----
    PRE_SNAPSHOT_RELEASE = "E422_PRE_SNAPSHOT_RELEASE"  # fail-closed temporal
    MODEL_CREATED_DATE_MISSING = "E422_MODEL_CREATED_DATE_MISSING"
    MODEL_CREATED_DATE_MISMATCH = "E422_MODEL_CREATED_DATE_MISMATCH"
    AFDB_VERSION_MISMATCH = "E422_AFDB_VERSION_MISMATCH"
    CHILD_URL_NOT_V6 = "E422_CHILD_URL_NOT_V6"
    ABSENT_FROM_SNAPSHOT = "E422_ABSENT_FROM_SNAPSHOT"
    SOURCE_CONTRACT_FAILURE = "E422_SOURCE_CONTRACT_FAILURE"
    MANIFEST_DIGEST_MISMATCH = "E422_MANIFEST_DIGEST_MISMATCH"
    MANIFEST_SCHEMA_INVALID = "E422_MANIFEST_SCHEMA_INVALID"

    # -- prediction outcomes (e422 §6/§7; recorded per checkpoint) -----------
    HOLD = "HOLD"  # P5: coverage >= nominal - 3 MCSE
    FAIL = "FAIL"  # P5: below
    REPLICATE = "REPLICATE"  # P6
    FALSIFY = "FALSIFY"  # P6 null replication (registered outcome)
    WIN = "WIN"  # P7a/P7b
    MIRRORED_WIN = "MIRRORED_WIN"  # P7a/P7b
    NO_DECISION = "NO_DECISION"  # P7a/P7b


#: Recording location for every status (checked against the enum above).
RECORDING_LOCATION = {
    "FIT_FAILURE": "evaluation artifact, per (repeat, level), plus R_eff",
    "UNDER_SUPPORTED_STRATUM": "evaluation artifact, per joint cell, with n",
    "DEFERRED_INSUFFICIENT_SUPPORT": "evaluation artifact, per statistic/comparison",
    "DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE": "close report (§6), per prediction",
    "PRE_SNAPSHOT_RELEASE": "batch checkpoint exclusion tally (fail-closed)",
    "MODEL_CREATED_DATE_MISSING": "manifest build / child-fetch receipt trail",
    "MODEL_CREATED_DATE_MISMATCH": "child-fetch receipt trail (typed failure)",
    "AFDB_VERSION_MISMATCH": "child-fetch receipt trail (source-contract class)",
    "CHILD_URL_NOT_V6": "manifest build receipt trail",
    "ABSENT_FROM_SNAPSHOT": "batch checkpoint exclusion tally (fail-closed)",
    "SOURCE_CONTRACT_FAILURE": "cross-fetch sha256 comparison in receipts",
    "MANIFEST_DIGEST_MISMATCH": "manifest load (tamper check)",
    "MANIFEST_SCHEMA_INVALID": "manifest build/load validation",
    "HOLD": "evaluation artifact, per level (P5)",
    "FAIL": "evaluation artifact, per level (P5)",
    "REPLICATE": "evaluation artifact, 0.90 only (P6)",
    "FALSIFY": "evaluation artifact, 0.90 only (P6)",
    "WIN": "evaluation artifact, per comparison (P7a/P7b)",
    "MIRRORED_WIN": "evaluation artifact, per comparison (P7a/P7b)",
    "NO_DECISION": "evaluation artifact, per comparison (P7a/P7b)",
}


def registry_digest() -> str:
    """SHA-256 over the canonical (sorted) value->location mapping."""
    import hashlib
    import json

    payload = {s.value: RECORDING_LOCATION[s.name] for s in E422Status}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


__all__ = ["E422Status", "RECORDING_LOCATION", "registry_digest"]
