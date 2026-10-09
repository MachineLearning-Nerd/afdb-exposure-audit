"""Tests for the e422 incremental presence manifest (e421-P §2.2 amended)."""

from __future__ import annotations

import copy
import datetime
import json
import os

import pytest

import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e422_manifest as em  # noqa: E402


T0 = "2026-09-10T14:23:05Z"
YESTERDAY = "2026-09-09"
TWO_DAYS_BEFORE = "2026-09-08"
SAME_DAY = "2026-09-10"
LATER = "2026-09-11"


def _model(accession="P12345", created=YESTERDAY, entity="1", version=6,
           cif_url=None):
    return {
        "uniprotAccession": accession,
        "modelCreatedDate": created,
        "modelEntityId": entity,
        "latestVersion": version,
        "cifUrl": cif_url
        or f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif",
    }


# --------------------------------------------------------------------------- #
# Canonical JSON.
# --------------------------------------------------------------------------- #
def test_canonical_json_is_sorted_and_compact():
    assert em.canonical_json_bytes({"b": 1, "a": 2}) == b'{"a":2,"b":1}'


def test_canonical_sha256_deterministic():
    assert em.canonical_sha256({"x": [1, 2]}) == em.canonical_sha256(
        {"x": [1, 2]}
    )


# --------------------------------------------------------------------------- #
# T0 + temporal rule.
# --------------------------------------------------------------------------- #
def test_parse_t0_handles_z_and_offsets():
    assert em.parse_t0("2026-09-10T14:23:05Z") == datetime.datetime(
        2026, 9, 10, 14, 23, 5, tzinfo=datetime.timezone.utc
    )
    assert em.parse_t0("2026-09-10T19:53:05+05:30") == em.parse_t0(
        "2026-09-10T14:23:05Z"
    )


def test_created_date_maps_to_midnight_utc():
    parsed = em.parse_created_date_utc("2026-09-08")
    assert parsed == datetime.datetime(2026, 9, 8, tzinfo=datetime.timezone.utc)


@pytest.mark.parametrize("bad", ["", "not-a-date", "2026-13-99"])
def test_unparseable_created_date_is_typed(bad):
    with pytest.raises(em.ManifestError) as excinfo:
        em.parse_created_date_utc(bad)
    assert excinfo.value.code == "E422_MODEL_CREATED_DATE_MISSING"


def test_temporal_accepts_day_strictly_before_t0():
    assert em.check_temporal(TWO_DAYS_BEFORE, T0) is not None


@pytest.mark.parametrize("created", [SAME_DAY, LATER])
def test_temporal_rejects_same_and_later_days_fail_closed(created):
    with pytest.raises(em.ManifestError) as excinfo:
        em.check_temporal(created, T0)
    assert excinfo.value.code == "E422_PRE_SNAPSHOT_RELEASE"


# --------------------------------------------------------------------------- #
# Record selection (registered rule).
# --------------------------------------------------------------------------- #
def test_select_record_earliest_created_then_lexicographic_entity():
    models = [
        _model(created=YESTERDAY, entity="b"),
        _model(created=TWO_DAYS_BEFORE, entity="z"),
        _model(created=TWO_DAYS_BEFORE, entity="a"),
    ]
    chosen = em.select_record(models)
    assert chosen["modelCreatedDate"] == TWO_DAYS_BEFORE
    assert chosen["modelEntityId"] == "a"


def test_select_record_skips_unparseable_and_raises_when_none():
    assert em.select_record([_model(created="garbage"), _model(created=YESTERDAY)])[
        "modelCreatedDate"
    ] == YESTERDAY
    with pytest.raises(em.ManifestError) as excinfo:
        em.select_record([_model(created="garbage")])
    assert excinfo.value.code == "E422_MODEL_CREATED_DATE_MISSING"


# --------------------------------------------------------------------------- #
# Record build + child-URL fence.
# --------------------------------------------------------------------------- #
def test_build_record_happy_path():
    record = em.build_record([_model()], T0)
    assert record == {
        "accession": "P12345",
        "v6_cif_url": "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif",
        "model_created_date": YESTERDAY,
        "afdb_version": 6,
        "model_entity_id": "1",
    }


def test_build_record_rejects_non_v6_child_url():
    with pytest.raises(em.ManifestError) as excinfo:
        em.build_record(
            [_model(cif_url="https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v7.cif")],
            T0,
        )
    assert excinfo.value.code == "E422_CHILD_URL_NOT_V6"


def test_build_record_enforces_temporal_check():
    with pytest.raises(em.ManifestError) as excinfo:
        em.build_record([_model(created=SAME_DAY)], T0)
    assert excinfo.value.code == "E422_PRE_SNAPSHOT_RELEASE"


# --------------------------------------------------------------------------- #
# Batch-time verification.
# --------------------------------------------------------------------------- #
def test_verify_against_manifest_accepts_match():
    record = em.build_record([_model()], T0)
    em.verify_against_manifest(record, {"modelCreatedDate": YESTERDAY,
                                        "latestVersion": 6})


def test_verify_against_manifest_flags_created_mismatch():
    record = em.build_record([_model()], T0)
    with pytest.raises(em.ManifestError) as excinfo:
        em.verify_against_manifest(record, {"modelCreatedDate": TWO_DAYS_BEFORE,
                                            "latestVersion": 6})
    assert excinfo.value.code == "E422_MODEL_CREATED_DATE_MISMATCH"


def test_verify_against_manifest_flags_missing_created_date():
    record = em.build_record([_model()], T0)
    with pytest.raises(em.ManifestError) as excinfo:
        em.verify_against_manifest(record, {"modelCreatedDate": "",
                                            "latestVersion": 6})
    assert excinfo.value.code == "E422_MODEL_CREATED_DATE_MISSING"


def test_verify_against_manifest_flags_version_mismatch():
    record = em.build_record([_model()], T0)
    with pytest.raises(em.ManifestError) as excinfo:
        em.verify_against_manifest(record, {"modelCreatedDate": YESTERDAY,
                                            "latestVersion": 7})
    assert excinfo.value.code == "E422_AFDB_VERSION_MISMATCH"


# --------------------------------------------------------------------------- #
# Builder + digest + persistence.
# --------------------------------------------------------------------------- #
def test_builder_admit_and_checkpoint_determinism(tmp_path):
    builder = em.ManifestBuilder(T0)
    first = builder.admit([_model()])
    assert first["accession"] == "P12345"
    builder.admit([_model(accession="Q99999", created=TWO_DAYS_BEFORE)])
    digest = builder.checkpoint()
    rebuilt = em.ManifestBuilder(T0)
    rebuilt.admit([_model(accession="Q99999", created=TWO_DAYS_BEFORE)])
    rebuilt.admit([_model()])
    assert rebuilt.checkpoint() == digest  # order-independent


def test_builder_conflicting_readmission_is_typed():
    builder = em.ManifestBuilder(T0)
    builder.admit([_model()])
    with pytest.raises(em.ManifestError) as excinfo:
        builder.admit([_model(created=TWO_DAYS_BEFORE)])
    assert excinfo.value.code == "E422_MANIFEST_SCHEMA_INVALID"


def test_builder_identical_readmission_is_idempotent():
    builder = em.ManifestBuilder(T0)
    builder.admit([_model()])
    builder.admit([_model()])
    assert len(builder.records) == 1


def test_builder_absent_accession_recorded_fail_closed():
    builder = em.ManifestBuilder(T0)
    note = builder.admit_absent("P00000", "metadata 404")
    assert "E422_ABSENT_FROM_SNAPSHOT" in note
    assert "P00000" not in builder.records  # exclusion, never a manifest row
    builder.admit([_model()])
    with pytest.raises(em.ManifestError) as excinfo:
        builder.admit_absent("P12345")  # already present -> typed error
    assert excinfo.value.code == "E422_MANIFEST_SCHEMA_INVALID"


def test_save_load_roundtrip_and_tamper_detection(tmp_path):
    builder = em.ManifestBuilder(T0)
    builder.admit([_model()])
    builder.admit([_model(accession="Q99999", created=TWO_DAYS_BEFORE)])
    path = str(tmp_path / "manifest.json")
    digest = builder.save(path)
    loaded = em.ManifestBuilder.load(path)
    assert loaded.checkpoint() == digest == builder.checkpoint()

    payload = json.loads(open(path, "rb").read().decode("utf-8"))
    payload["records"][0]["model_created_date"] = TWO_DAYS_BEFORE
    tampered = str(tmp_path / "tampered.json")
    open(tampered, "wb").write(em.canonical_json_bytes(payload))
    with pytest.raises(em.ManifestError) as excinfo:
        em.ManifestBuilder.load(tampered)
    assert excinfo.value.code == "E422_MANIFEST_DIGEST_MISMATCH"


def test_builder_rejects_malformed_t0_at_construction():
    with pytest.raises(ValueError):
        em.ManifestBuilder("not-a-time")


def test_admit_propagates_pre_snapshot_release():
    builder = em.ManifestBuilder(T0)
    with pytest.raises(em.ManifestError) as excinfo:
        builder.admit([_model(created=SAME_DAY)])
    assert excinfo.value.code == "E422_PRE_SNAPSHOT_RELEASE"
    assert builder.records == {}  # fail-closed: nothing stored
