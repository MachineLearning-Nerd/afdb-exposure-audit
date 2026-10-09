"""e422 snapshot presence-manifest (incremental per-batch checkpoints).

Registered basis: e421-P §2.2 as amended for e422 at the gate-(b)
disposition (2026-09-09).  The complete AFDB-v6 accession snapshot is not
enumerable over the registered allowlist lanes (no bulk index endpoint is
allowlisted; the per-accession metadata lane cannot enumerate ~2e8
entries within the registered rate limits), so the manifest is built
INCREMENTALLY: each census accession enters the manifest at its first
appearance, with a FAIL-CLOSED temporal check that preserves the §2.1
T0-snapshot semantics:

* the AFDB metadata response's ``modelCreatedDate`` is date-only ``D``;
  it is mapped to the instant ``D T00:00:00Z``;
* ``D >= T0.date()`` (UTC) — INCLUDING the same UTC day as T0 — is the
  typed failure ``PRE_SNAPSHOT_RELEASE`` (fail-closed exclusion): a model
  created any time during T0's UTC day may postdate the T0 instant, so
  same-day is rejected by construction;
* an accession whose metadata fetch fails or returns no v6 model is
  ``ABSENT_FROM_SNAPSHOT`` (fail-closed; no substitution).

Manifest records are ``{accession, v6_cif_url, model_created_date,
afdb_version, model_entity_id}``.  Multiple v6 model records in one
response bind the registered selection rule: earliest
``modelCreatedDate``, then lexicographic ``modelEntityId`` (the e420
grammar amendment).  Manifest digest = SHA-256 over canonical JSON
(``sort_keys=True, separators=(',',':')``), checkpointed per batch and
OTS-anchored with the batch checkpoint (e422 §6).

Typed failure vocabulary (this module):
``PRE_SNAPSHOT_RELEASE``, ``MODEL_CREATED_DATE_MISSING``,
``MODEL_CREATED_DATE_MISMATCH``, ``AFDB_VERSION_MISMATCH`` (a
source-contract disagreement: the registered v6 release is uniform, so a
manifest/fetch version disagreement is ``SOURCE_CONTRACT_FAILURE``
class), ``CHILD_URL_NOT_V6``, ``ABSENT_FROM_SNAPSHOT``,
``MANIFEST_DIGEST_MISMATCH`` (tamper/load check).
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from enum import Enum
from typing import Any, Iterable, Mapping

import e420_census as census


# --------------------------------------------------------------------------- #
# Typed failures.
# --------------------------------------------------------------------------- #
class ManifestFailureCode(str, Enum):
    """Typed manifest failure codes (e422 gate-(b) disposition)."""

    PRE_SNAPSHOT_RELEASE = "E422_PRE_SNAPSHOT_RELEASE"
    MODEL_CREATED_DATE_MISSING = "E422_MODEL_CREATED_DATE_MISSING"
    MODEL_CREATED_DATE_MISMATCH = "E422_MODEL_CREATED_DATE_MISMATCH"
    AFDB_VERSION_MISMATCH = "E422_AFDB_VERSION_MISMATCH"
    CHILD_URL_NOT_V6 = "E422_CHILD_URL_NOT_V6"
    ABSENT_FROM_SNAPSHOT = "E422_ABSENT_FROM_SNAPSHOT"
    MANIFEST_DIGEST_MISMATCH = "E422_MANIFEST_DIGEST_MISMATCH"
    MANIFEST_SCHEMA_INVALID = "E422_MANIFEST_SCHEMA_INVALID"


class ManifestError(RuntimeError):
    """Typed manifest error (fail-closed; never silently dropped)."""

    def __init__(self, code: ManifestFailureCode | str, detail: str = "") -> None:
        self.code = code.value if isinstance(code, ManifestFailureCode) else str(code)
        self.detail = detail
        super().__init__(f"{self.code}: {detail}")


# --------------------------------------------------------------------------- #
# Canonical JSON + digest.
# --------------------------------------------------------------------------- #
def canonical_json_bytes(obj: Any) -> bytes:
    """Canonical JSON: sort_keys=True, separators=(',',':'), UTF-8."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")


def canonical_sha256(obj: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(obj)).hexdigest()


# --------------------------------------------------------------------------- #
# Temporal rule.
# --------------------------------------------------------------------------- #
def parse_created_date_utc(model_created_date: str) -> datetime.datetime:
    """Date-only ``D`` -> aware UTC instant ``D T00:00:00Z``.

    Raises ``MODEL_CREATED_DATE_MISSING`` on an unparseable value.
    """
    text = (model_created_date or "").strip()
    if not text:
        raise ManifestError(
            ManifestFailureCode.MODEL_CREATED_DATE_MISSING, "empty modelCreatedDate"
        )
    try:
        day = datetime.date.fromisoformat(text[:10])
    except ValueError as exc:
        raise ManifestError(
            ManifestFailureCode.MODEL_CREATED_DATE_MISSING,
            f"unparseable modelCreatedDate {model_created_date!r}",
        ) from exc
    return datetime.datetime(day.year, day.month, day.day, tzinfo=datetime.timezone.utc)


def check_temporal(model_created_date: str, t0_utc: str) -> datetime.datetime:
    """Fail-closed §2.1 temporal check.

    Accepts only ``D < T0.date()`` (UTC day strictly before T0's UTC day);
    ``D >= T0.date()`` — including the same UTC day — raises
    ``PRE_SNAPSHOT_RELEASE``.  Returns the mapped instant on success.
    """
    t0 = parse_t0(t0_utc)
    release = parse_created_date_utc(model_created_date)
    if release.date() >= t0.date():
        raise ManifestError(
            ManifestFailureCode.PRE_SNAPSHOT_RELEASE,
            f"modelCreatedDate {model_created_date} not strictly before T0 day "
            f"{t0.date().isoformat()}",
        )
    return release


def parse_t0(t0_utc: str) -> datetime.datetime:
    """Parse the frozen second-level UTC T0 instant (``...Z`` suffix)."""
    text = (t0_utc or "").strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    parsed = datetime.datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=datetime.timezone.utc)
    return parsed.astimezone(datetime.timezone.utc)


# --------------------------------------------------------------------------- #
# Record selection + validation.
# --------------------------------------------------------------------------- #
MANIFEST_RECORD_KEYS = frozenset(
    {"accession", "v6_cif_url", "model_created_date", "afdb_version", "model_entity_id"}
)


def select_record(models: Iterable[Mapping[str, Any]]) -> dict:
    """Registered selection rule over v6 model records.

    Earliest ``modelCreatedDate``, then lexicographic ``modelEntityId``
    (the e420 grammar amendment).  Records lacking a parseable
    ``modelCreatedDate`` are skipped; an empty/invalid set raises
    ``MODEL_CREATED_DATE_MISSING``.
    """
    parsed: list[tuple[str, str, Mapping[str, Any]]] = []
    for model in models:
        created = (model.get("modelCreatedDate") or "").strip()
        try:
            parse_created_date_utc(created)
        except ManifestError:
            continue
        entity = str(model.get("modelEntityId") or "")
        parsed.append((created, entity, model))
    if not parsed:
        raise ManifestError(
            ManifestFailureCode.MODEL_CREATED_DATE_MISSING,
            "no v6 model record with a parseable modelCreatedDate",
        )
    parsed.sort(key=lambda item: (item[0], item[1]))
    return dict(parsed[0][2])


def _expected_child_url(accession: str) -> str:
    return f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif"


def build_record(models: Iterable[Mapping[str, Any]], t0_utc: str) -> dict:
    """Validate + select + temporally check one metadata response -> record.

    Raises, in order: ``MODEL_CREATED_DATE_MISSING`` (no parseable record),
    ``PRE_SNAPSHOT_RELEASE`` (temporal check), ``CHILD_URL_NOT_V6`` (the
    observed ``cifUrl`` is not the registered v6 child shape).
    """
    chosen = select_record(models)
    check_temporal(chosen["modelCreatedDate"], t0_utc)
    accession = str(chosen.get("uniprotAccession") or chosen.get("accession") or "").strip()
    cif_url = str(chosen.get("cifUrl") or "").strip()
    if not cif_url or cif_url != _expected_child_url(accession):
        raise ManifestError(
            ManifestFailureCode.CHILD_URL_NOT_V6,
            f"cifUrl {cif_url!r} is not the registered v6 child for {accession!r}",
        )
    return {
        "accession": accession,
        "v6_cif_url": cif_url,
        "model_created_date": str(chosen["modelCreatedDate"]),
        "afdb_version": int(chosen.get("latestVersion")),
        "model_entity_id": str(chosen.get("modelEntityId") or ""),
    }


def verify_against_manifest(
    manifest_record: Mapping[str, Any], fetched: Mapping[str, Any]
) -> None:
    """Batch-time verification (e421-P §2.2 point 3).

    ``fetched`` carries the metadata values observed at a verification
    fetch.  A parseable-created-date disagreement raises
    ``MODEL_CREATED_DATE_MISMATCH``; a missing created date raises
    ``MODEL_CREATED_DATE_MISSING``; an ``afdb_version`` disagreement
    raises ``AFDB_VERSION_MISMATCH`` (source-contract class).
    """
    created = (fetched.get("modelCreatedDate") or "").strip()
    try:
        parse_created_date_utc(created)
    except ManifestError as exc:
        raise ManifestError(
            ManifestFailureCode.MODEL_CREATED_DATE_MISSING, exc.detail
        ) from exc
    if created != str(manifest_record["model_created_date"]):
        raise ManifestError(
            ManifestFailureCode.MODEL_CREATED_DATE_MISMATCH,
            f"manifest {manifest_record['model_created_date']} != fetched {created}",
        )
    version = fetched.get("latestVersion")
    try:
        version = int(version)
    except (TypeError, ValueError) as exc:
        raise ManifestError(
            ManifestFailureCode.AFDB_VERSION_MISMATCH,
            f"unparseable latestVersion {version!r}",
        ) from exc
    if version != int(manifest_record["afdb_version"]):
        raise ManifestError(
            ManifestFailureCode.AFDB_VERSION_MISMATCH,
            f"manifest v{manifest_record['afdb_version']} != fetched v{version}",
        )


# --------------------------------------------------------------------------- #
# Manifest object.
# --------------------------------------------------------------------------- #
MANIFEST_SCHEMA = "e422-presence-manifest-v1"


class ManifestBuilder:
    """Incremental presence manifest with canonical-JSON SHA-256 digests."""

    def __init__(self, t0_utc: str) -> None:
        self.t0_utc = t0_utc
        parse_t0(t0_utc)  # fail fast on a malformed frozen T0
        self.records: dict[str, dict] = {}

    def admit(self, models: Iterable[Mapping[str, Any]]) -> dict:
        """Validate + append one accession; returns the stored record."""
        record = build_record(models, self.t0_utc)
        accession = record["accession"]
        prior = self.records.get(accession)
        if prior is not None and prior != record:
            raise ManifestError(
                ManifestFailureCode.MANIFEST_SCHEMA_INVALID,
                f"conflicting re-admission for {accession}",
            )
        self.records[accession] = record
        return record

    def admit_absent(self, accession: str, detail: str = "") -> str:
        """Record a fail-closed ``ABSENT_FROM_SNAPSHOT`` exclusion."""
        if accession in self.records:
            raise ManifestError(
                ManifestFailureCode.MANIFEST_SCHEMA_INVALID,
                f"ABSENT_FROM_SNAPSHOT for already-present {accession}",
            )
        return ManifestFailureCode.ABSENT_FROM_SNAPSHOT.value + (
            f" ({detail})" if detail else ""
        )

    def digest_obj(self) -> dict:
        """The digested object: schema, T0, and sorted records."""
        return {
            "schema": MANIFEST_SCHEMA,
            "t0_utc": self.t0_utc,
            "records": [self.records[k] for k in sorted(self.records)],
        }

    def checkpoint(self) -> str:
        """Canonical-JSON SHA-256 of the current manifest content."""
        return canonical_sha256(self.digest_obj())

    def save(self, path: str) -> str:
        """Atomic write of {schema, t0_utc, records, sha256}."""
        payload = dict(self.digest_obj())
        payload["sha256"] = self.checkpoint()
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=directory, suffix=".tmp")
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(canonical_json_bytes(payload))
            os.replace(tmp, path)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)
        return payload["sha256"]

    @classmethod
    def load(cls, path: str) -> "ManifestBuilder":
        """Load + digest-verify (``MANIFEST_DIGEST_MISMATCH`` on tamper)."""
        with open(path, "rb") as handle:
            payload = json.loads(handle.read().decode("utf-8"))
        for key in ("schema", "t0_utc", "records", "sha256"):
            if key not in payload:
                raise ManifestError(
                    ManifestFailureCode.MANIFEST_SCHEMA_INVALID, f"missing key {key}"
                )
        if payload["schema"] != MANIFEST_SCHEMA:
            raise ManifestError(
                ManifestFailureCode.MANIFEST_SCHEMA_INVALID,
                f"schema {payload['schema']!r}",
            )
        builder = cls(payload["t0_utc"])
        builder.records = {r["accession"]: r for r in payload["records"]}
        if builder.checkpoint() != payload["sha256"]:
            raise ManifestError(
                ManifestFailureCode.MANIFEST_DIGEST_MISMATCH,
                f"digest mismatch for {path}",
            )
        return builder


__all__ = [
    "ManifestFailureCode",
    "ManifestError",
    "ManifestBuilder",
    "MANIFEST_SCHEMA",
    "MANIFEST_RECORD_KEYS",
    "canonical_json_bytes",
    "canonical_sha256",
    "parse_t0",
    "parse_created_date_utc",
    "check_temporal",
    "select_record",
    "build_record",
    "verify_against_manifest",
    "census",
]
