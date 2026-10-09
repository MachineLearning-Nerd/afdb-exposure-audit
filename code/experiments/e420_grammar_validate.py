"""E420 gate-(e) phase (i) — offline grammar-validation harness.

Strict envelope parsers/validators for the eight registered e420 source roles
(``docs/e420_registration.md`` §3.5), bound to retained contract bytes only:

- PDBe v2 primary uniprot_mapping envelope   (E420_PDBE_MAPPING_PRIMARY_V1)
- PDBe v2 segment mappings envelope          (E420_PDBE_MAPPING_SEGMENT_V1)
- AFDB prediction NewEntrySummary            (E420_AFDB_METADATA_V1)
- AFDB coordinate child cifUrl binding       (E420_AFDB_COORDINATE_CHILD_V1)
- ATLAS analysis ZIP member grammar          (E420_ATLAS_ANALYSIS_V1)
- ATLAS metadata envelope                    (E420_ATLAS_METADATA_V1)
- ATLAS parsable census member               (E420_ATLAS_CENSUS_UNIVERSE_V1)
- reference entry-files mmCIF reader         (E420_REFERENCE_COORDINATE_V1)

Every validator: raw-bytes hash/length binding, typed failures as frozen
dataclasses (no native exception crosses the boundary; validators return
:class:`GrammarOutcome` and never raise), and exact-key per-row provenance
envelopes per ``data/e418_retained_provenance_schema_20260830.json`` (v27).

This phase is OFFLINE: retained bytes and fixture-derived synthetic response
bodies only.  No network, source, API, or data fetch of any kind; no
candidate/coordinate/RMSF/pLDDT values are read or exported (pLDDT-family
fields are type-checked for presence only and never appear in any payload;
ATLAS RMSF lexemes are grammar-checked without exporting values).
"""
from __future__ import annotations

import datetime as _dt
import functools
import hashlib
import io
import json
import re
import zipfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

# --------------------------------------------------------------------------- #
# Retained contract pins (docs/e420_registration.md §3.2; read-only inputs).  #
# --------------------------------------------------------------------------- #

PROVENANCE_SCHEMA_REL = "data/e418_retained_provenance_schema_20260830.json"
PROVENANCE_SCHEMA_SHA256 = "a223f271162fa921cffc23d4b040a603403b76df288cbf439021984df6b14efe"
PROVENANCE_SCHEMA_MARKER = "e418-provenance-schema-generator-v27"
PROVENANCE_SCHEMA_NAME = "e418-retained-provenance-schema-v27"

PDBe_V2_OPENAPI_REL = "data/source_contracts/pdbe_v2_openapi_20260830.json"
PDBe_V2_OPENAPI_SHA256 = "e86ab8d3c2c9b6865f158526c93b742db344fa62ec7277dea1e4d17e6852899d"
AFDB_OPENAPI_REL = "data/source_contracts/afdb_api_openapi_20260830.json"
AFDB_OPENAPI_SHA256 = "714607265fd8edc581baf28df038ea804d96d871baa6034ff60d22d0cf893163"
ATLAS_OPENAPI_REL = "data/source_contracts/atlas_api_openapi_20260830.json"
ATLAS_OPENAPI_SHA256 = "824150cc10bd365c82176bc11884a168d02a64029da65b475620a44054d3ee44"
AFDB_OPTIONAL_TABLE_REL = "data/source_contracts/afdb_optional_types_generated_20260830.json"
AFDB_OPTIONAL_TABLE_SHA256 = "060de272d3a0de2d0d9346f67cb03efab49f00db15363e90f26811a791176eb3"
AFDB_OPTIONAL_TABLE_CONTENT_SHA256 = "32d205d21376616dd8fbad4d13963fc41e10c00edeadc6588b08f16994da4033"

ATLAS_CENSUS_ZIP_REL = "data/source_contracts/atlas_parsable_full_20260830.zip"
ATLAS_CENSUS_ZIP_SHA256 = "c9cba2b7190676814bb83daf11fcfe92ad0c3a3cdaca4eb4782ca6b98c67b75e"
ATLAS_CENSUS_MEMBER_PATH = "ATLAS_parsable_latest/2023_03_09_ATLAS_pdb.txt"
ATLAS_CENSUS_MEMBER_SHA256 = "667a5ebc28cd3f160a4fa313fe713d337dc1148aca877216dd9153f44a68cd0b"

EMITTED_BY = "e420-grammar-validate-v1"
FOLD_SALT = "g1-successor-fold-v1:"
SIFTS_GAP_DERIVATION_SOURCE = "SIFTS_GAP"
GAP_ROW_UNIPROT_NUM: None = None
HEX64_RE = re.compile(r"[0-9a-f]{64}")
PDB_ID_RE = re.compile(r"[0-9A-Za-z]{4}")
CHAIN_RE = re.compile(r"[A-Za-z0-9]+")
MODEL_CREATED_DATE_FORMAT = "%Y-%m-%dT%H:%M:%SZ"
REVISION_DATE_RE = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}")
MD5_RE = re.compile(r"[0-9a-f]{32}")
SEQ_LEXEME_RE = re.compile(r"-?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?")

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# --------------------------------------------------------------------------- #
# Typed failure vocabulary.                                                   #
# --------------------------------------------------------------------------- #


class SourceRole(str, Enum):
    """The eight registered e420 source roles (registration §3.5 literals)."""

    ATLAS_CENSUS_UNIVERSE = "E420_ATLAS_CENSUS_UNIVERSE_V1"
    PDBE_MAPPING_PRIMARY = "E420_PDBE_MAPPING_PRIMARY_V1"
    PDBE_MAPPING_SEGMENT = "E420_PDBE_MAPPING_SEGMENT_V1"
    AFDB_METADATA = "E420_AFDB_METADATA_V1"
    AFDB_COORDINATE_CHILD = "E420_AFDB_COORDINATE_CHILD_V1"
    ATLAS_ANALYSIS = "E420_ATLAS_ANALYSIS_V1"
    ATLAS_METADATA = "E420_ATLAS_METADATA_V1"
    REFERENCE_COORDINATE = "E420_REFERENCE_COORDINATE_V1"


class SegmentMode(str, Enum):
    """Segment-envelope use modes (registration §4.4 guards)."""

    RESOLUTION = "resolution"   # identity >= 0.90 and coverage >= 0.80
    FALLBACK = "fallback"       # identity == 100.0 exact


class GrammarFailureCode(str, Enum):
    """Every typed failure class the harness can emit."""

    BOUNDARY_INVALID = "boundary_invalid"
    RAW_BYTES_NOT_BYTES = "raw_bytes_not_bytes"
    RAW_BYTES_EMPTY = "raw_bytes_empty"
    RAW_BYTES_HASH_MISMATCH = "raw_bytes_hash_mismatch"
    RAW_LEN_MISMATCH = "raw_len_mismatch"
    PROVENANCE_SCHEMA_PIN_MISMATCH = "provenance_schema_pin_mismatch"
    CONTRACT_ARTIFACT_PIN_MISMATCH = "contract_artifact_pin_mismatch"
    ENVELOPE_KEYS_INVALID = "envelope_keys_invalid"
    ENVELOPE_SHA256_FORMAT_INVALID = "envelope_sha256_format_invalid"
    ENVELOPE_MEMBER_KEYS_INVALID = "envelope_member_keys_invalid"
    ENVELOPE_MEMBER_HASH_INVALID = "envelope_member_hash_invalid"
    ENVELOPE_CONTEXT_INVALID = "envelope_context_invalid"
    ENVELOPE_GAP_DERIVATION_INVALID = "envelope_gap_derivation_invalid"
    PDBE_PRIMARY_ROOT_KEY_INVALID = "pdbe_primary_root_key_invalid"
    PDBE_PRIMARY_ENTRY_KEYS_INVALID = "pdbe_primary_entry_keys_invalid"
    PDBE_PRIMARY_SEQUENCE_INVALID = "pdbe_primary_sequence_invalid"
    PDBE_PRIMARY_LENGTH_MISMATCH = "pdbe_primary_length_mismatch"
    PDBE_PRIMARY_DATA_EMPTY = "pdbe_primary_data_empty"
    PDBE_PRIMARY_RESIDUE_KEYS_INVALID = "pdbe_primary_residue_keys_invalid"
    PDBE_PRIMARY_RESIDUE_BOUNDS_INVALID = "pdbe_primary_residue_bounds_invalid"
    PDBE_PRIMARY_CODE_INVALID = "pdbe_primary_code_invalid"
    PDBE_PRIMARY_ADDITIONAL_KEYS_INVALID = "pdbe_primary_additional_keys_invalid"
    PDBE_PRIMARY_ENTITY_TYPE_INVALID = "pdbe_primary_entity_type_invalid"
    PDBE_PRIMARY_ENTITY_MISMATCH = "pdbe_primary_entity_mismatch"
    PDBE_PRIMARY_NOT_JSON = "pdbe_primary_not_json"
    PDBE_SEGMENT_ROOT_KEY_INVALID = "pdbe_segment_root_key_invalid"
    PDBE_SEGMENT_NOT_JSON = "pdbe_segment_not_json"
    PDBE_SEGMENT_OBJECT_KEYS_INVALID = "pdbe_segment_object_keys_invalid"
    PDBE_SEGMENT_MAPPING_KEYS_INVALID = "pdbe_segment_mapping_keys_invalid"
    PDBE_SEGMENT_MAPPING_TYPE_INVALID = "pdbe_segment_mapping_type_invalid"
    PDBE_SEGMENT_RESIDUE_KEYS_INVALID = "pdbe_segment_residue_keys_invalid"
    PDBE_SEGMENT_BOUNDS_INVALID = "pdbe_segment_bounds_invalid"
    PDBE_SEGMENT_ACCESSION_ABSENT = "pdbe_segment_accession_absent"
    PDBE_SEGMENT_ENTITY_AMBIGUOUS = "pdbe_segment_entity_ambiguous"
    PDBE_SEGMENT_ENTITY_RESOLUTION_MISMATCH = "pdbe_segment_entity_resolution_mismatch"
    PDBE_SEGMENT_THRESHOLD_FAILURE = "pdbe_segment_threshold_failure"
    AFDB_METADATA_NOT_JSON = "afdb_metadata_not_json"
    AFDB_METADATA_CARDINALITY_INVALID = "afdb_metadata_cardinality_invalid"
    AFDB_FIELD_REQUIRED_MISSING = "afdb_field_required_missing"
    AFDB_FIELD_UNKNOWN = "afdb_field_unknown"
    AFDB_FIELD_SUNSET = "afdb_field_sunset"
    AFDB_FIELD_TYPE_INVALID = "afdb_field_type_invalid"
    AFDB_MODEL_CREATED_DATE_INVALID = "afdb_model_created_date_invalid"
    AFDB_VERSION_INVALID = "afdb_version_invalid"
    AFDB_SEQUENCE_OFFSET_INVALID = "afdb_sequence_offset_invalid"
    AFDB_SEQUENCE_OFFSET_MISMATCH = "afdb_sequence_offset_mismatch"
    AFDB_SEQUENCE_CHECKSUM_INVALID = "afdb_sequence_checksum_invalid"
    AFDB_SEQUENCE_CHECKSUM_MISMATCH = "afdb_sequence_checksum_mismatch"
    AFDB_CIF_URL_MISSING = "afdb_cif_url_missing"
    AFDB_ACCESSION_MISMATCH = "afdb_accession_mismatch"
    AFDB_CHILD_URL_MISMATCH = "afdb_child_url_mismatch"
    AFDB_CHILD_URL_FORBIDDEN_SUBSTITUTE = "afdb_child_url_forbidden_substitute"
    AFDB_CHILD_BYTES_EMPTY = "afdb_child_bytes_empty"
    AFDB_CHILD_MODEL_CIF_HEADER_INVALID = "afdb_child_model_cif_header_invalid"
    AFDB_TEMPORAL_NOT_STRICT = "afdb_temporal_not_strict"
    ATLAS_CENSUS_ZIP_INVALID = "atlas_census_zip_invalid"
    ATLAS_CENSUS_MEMBER_PATH_INVALID = "atlas_census_member_path_invalid"
    ATLAS_CENSUS_ROW_GRAMMAR_INVALID = "atlas_census_row_grammar_invalid"
    ATLAS_ARCHIVE_ZIP_INVALID = "atlas_archive_zip_invalid"
    ATLAS_ARCHIVE_MEMBER_SET_INVALID = "atlas_archive_member_set_invalid"
    ATLAS_ARCHIVE_MEMBER_CRC_INVALID = "atlas_archive_member_crc_invalid"
    ATLAS_RMSF_HEADER_INVALID = "atlas_rmsf_header_invalid"
    ATLAS_CORRESP_HEADER_INVALID = "atlas_corresp_header_invalid"
    ATLAS_RMSF_ROW_INVALID = "atlas_rmsf_row_invalid"
    ATLAS_CORRESP_ROW_INVALID = "atlas_corresp_row_invalid"
    ATLAS_CORRESP_FLAG_INVALID = "atlas_corresp_flag_invalid"
    ATLAS_CORRESP_CIF_PDB_MISMATCH = "atlas_corresp_cif_pdb_mismatch"
    ATLAS_ROW_COUNT_MISMATCH = "atlas_row_count_mismatch"
    ATLAS_GAP_RATIO_EXCEEDED = "atlas_gap_ratio_exceeded"
    ATLAS_METADATA_NOT_JSON = "atlas_metadata_not_json"
    ATLAS_METADATA_CHAIN_KEY_INVALID = "atlas_metadata_chain_key_invalid"
    ATLAS_METADATA_KEYS_INVALID = "atlas_metadata_keys_invalid"
    ATLAS_METADATA_IDENTITY_MISMATCH = "atlas_metadata_identity_mismatch"
    REFERENCE_CIF_DECODE_FAILED = "reference_cif_decode_failed"
    REFERENCE_CIF_TOKEN_INVALID = "reference_cif_token_invalid"
    REFERENCE_CIF_DATA_BLOCK_INVALID = "reference_cif_data_block_invalid"
    REFERENCE_HISTORY_MISSING = "reference_history_missing"
    REFERENCE_HISTORY_AMBIGUOUS = "reference_history_ambiguous"
    REFERENCE_REVISION_DATE_INVALID = "reference_revision_date_invalid"
    RELEASE_LOCK_MISMATCH = "release_lock_mismatch"
    FOLD_PREIMAGE_AGGREGATE_MISMATCH = "fold_preimage_aggregate_mismatch"
    FOLD_PREIMAGE_ROW_INVALID = "fold_preimage_row_invalid"


@dataclass(frozen=True)
class GrammarFailure:
    """Typed failure (dataclass, never raised across the boundary)."""

    code: GrammarFailureCode
    role: SourceRole
    detail: str = ""
    expected_sha256: str | None = None
    observed_sha256: str | None = None


@dataclass(frozen=True)
class GrammarOutcome:
    """Result envelope of every boundary validator."""

    role: SourceRole
    ok: bool
    raw_sha256: str
    raw_len: int
    payload: Mapping[str, Any] | None = None
    failure: GrammarFailure | None = None


class E420GrammarPinError(ValueError):
    """Typed loader error (loaders only; validators map it to a typed failure)."""

    def __init__(self, message: str, code: GrammarFailureCode = GrammarFailureCode.CONTRACT_ARTIFACT_PIN_MISMATCH) -> None:
        self.code = code
        super().__init__(message)


class _GrammarError(Exception):
    """Internal typed error, converted to GrammarFailure at the boundary."""

    def __init__(self, code: GrammarFailureCode, detail: str = "", **kw: Any) -> None:
        self.code = code
        self.detail = detail
        self.extra = kw
        super().__init__(code.value, detail)


def _fail(code: GrammarFailureCode, detail: str = "", **kw: Any) -> Any:
    raise _GrammarError(code, detail, **kw)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_digest(value: Any) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    return _sha256(canonical.encode("utf-8"))


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _require_bytes(raw: Any, expected_raw_sha256: str | None, expected_raw_len: int | None) -> bytes:
    if isinstance(raw, (bytearray, memoryview)):
        raw = bytes(raw)
    if not isinstance(raw, bytes):
        _fail(GrammarFailureCode.RAW_BYTES_NOT_BYTES, type(raw).__name__)
    if not raw:
        _fail(GrammarFailureCode.RAW_BYTES_EMPTY)
    if expected_raw_len is not None:
        if not _is_int(expected_raw_len) or len(raw) != expected_raw_len:
            _fail(GrammarFailureCode.RAW_LEN_MISMATCH, f"observed={len(raw)}", raw_len=len(raw), raw_sha256=_sha256(raw))
    if expected_raw_sha256 is not None:
        observed = _sha256(raw)
        if not (isinstance(expected_raw_sha256, str) and HEX64_RE.fullmatch(expected_raw_sha256) and observed == expected_raw_sha256):
            _fail(GrammarFailureCode.RAW_BYTES_HASH_MISMATCH, "retained-bytes binding",
                  raw_len=len(raw), raw_sha256=observed, expected_sha256=str(expected_raw_sha256))
    return raw


def _grammar_boundary(role: SourceRole):
    """Boundary decorator: every failure becomes a typed GrammarOutcome."""

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> GrammarOutcome:
            empty_sha = _sha256(b"")
            try:
                return func(*args, **kwargs)
            except _GrammarError as exc:
                return GrammarOutcome(
                    role=role, ok=False,
                    raw_sha256=exc.extra.get("raw_sha256", empty_sha),
                    raw_len=exc.extra.get("raw_len", 0),
                    payload=None,
                    failure=GrammarFailure(code=exc.code, role=role, detail=exc.detail,
                                           expected_sha256=exc.extra.get("expected_sha256"),
                                           observed_sha256=exc.extra.get("raw_sha256")),
                )
            except E420GrammarPinError as exc:
                return GrammarOutcome(role=role, ok=False, raw_sha256=empty_sha, raw_len=0, payload=None,
                                      failure=GrammarFailure(code=exc.code, role=role, detail=str(exc)))
            except OSError as exc:
                return GrammarOutcome(role=role, ok=False, raw_sha256=empty_sha, raw_len=0, payload=None,
                                      failure=GrammarFailure(code=GrammarFailureCode.CONTRACT_ARTIFACT_PIN_MISMATCH,
                                                             role=role, detail=type(exc).__name__))
            except Exception as exc:  # noqa: BLE001 - boundary converts everything
                return GrammarOutcome(
                    role=role, ok=False, raw_sha256=empty_sha, raw_len=0, payload=None,
                    failure=GrammarFailure(code=GrammarFailureCode.BOUNDARY_INVALID, role=role,
                                           detail=f"{type(exc).__name__}"),
                )
        return wrapper
    return decorator


# --------------------------------------------------------------------------- #
# Retained contract loaders (pinned bytes; cached).                           #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ProvenanceSchemaBinding:
    """Verified v27 provenance schema binding + embedded PDBe field schemas."""

    schema_sha256: str
    required_keys: frozenset[str]
    conditional_keys: frozenset[str]
    member_role_keys: tuple[str, ...]
    primary_pdb_keys: frozenset[str]
    primary_data_keys: frozenset[str]
    primary_residue_required: frozenset[str]
    primary_residue_optional: frozenset[str]
    primary_additional_keys: frozenset[str]
    primary_residue_additional_keys: frozenset[str]
    segment_object_keys: frozenset[str]
    segment_mapping_keys: frozenset[str]
    segment_residue_required: frozenset[str]
    segment_residue_optional: frozenset[str]


@dataclass(frozen=True)
class AfdbGrammarContract:
    """Verified AFDB NewEntrySummary contract (retained OpenAPI + types table)."""

    openapi_sha256: str
    required_fields: frozenset[str]
    active_optional: frozenset[str]
    sunset_fields: frozenset[str]
    type_specs: Mapping[str, Any]


@dataclass(frozen=True)
class AtlasMetadataContract:
    """Verified ATLAS MetadataResponse contract (retained OpenAPI)."""

    openapi_sha256: str
    required_fields: frozenset[str]
    type_specs: Mapping[str, Any]


def _load_json(path: Any, expected_sha256: str) -> Any:
    raw = path.read_bytes()
    if _sha256(raw) != expected_sha256:
        raise E420GrammarPinError(f"{path.name}: sha256 pin mismatch")
    return json.loads(raw.decode("utf-8"))


@functools.lru_cache(maxsize=8)
def _load_provenance_schema_binding_cached(root_str: str) -> ProvenanceSchemaBinding:
    path = _project_path(root_str) / PROVENANCE_SCHEMA_REL
    if not path.is_file() or _sha256(path.read_bytes()) != PROVENANCE_SCHEMA_SHA256:
        raise E420GrammarPinError(f"{path.name}: sha256 pin mismatch",
                                  code=GrammarFailureCode.PROVENANCE_SCHEMA_PIN_MISMATCH)
    schema = json.loads(path.read_bytes().decode("utf-8"))
    if schema.get("schema") != PROVENANCE_SCHEMA_NAME or schema.get("parser_version_marker") != PROVENANCE_SCHEMA_MARKER:
        raise E420GrammarPinError(f"{path.name}: schema/marker pin mismatch",
                                  code=GrammarFailureCode.PROVENANCE_SCHEMA_PIN_MISMATCH)
    contract = schema["normalized_row_contract"]["per_row_provenance_required"]
    required_keys = frozenset(contract)
    conditional = frozenset({"gap_derivation", "modelcif_projection_sha256"})
    member_keys = tuple(schema["normalized_row_contract"]["per_row_provenance_required"]["source_member_hashes"])
    primary = schema["openapi_field_schemas"]["primary_route"]["schemas"]
    segment = schema["openapi_field_schemas"]["segment_route"]["schemas"]

    def keys_of(entry_map: Mapping[str, Any], name: str) -> tuple[frozenset[str], frozenset[str]]:
        fields = entry_map[name]["fields"]
        req = frozenset(k for k, v in fields.items() if v.get("required_field") is True)
        opt = frozenset(k for k, v in fields.items() if v.get("required_field") is False)
        return req, opt

    pdb_req, _ = keys_of(primary, "PDB")
    data_req, _ = keys_of(primary, "Data")
    res_req, res_opt = keys_of(primary, "Residue")
    add_keys = frozenset(primary["AdditionalData"]["fields"])
    resadd_keys = frozenset(primary["ResidueAdditionalData"]["fields"])
    seg_obj_keys = frozenset(segment["PdbUniprotSimpleMappings"]["fields"])
    seg_req, _ = keys_of(segment, "PdbUniprotSimpleMapping")
    segres_req, segres_opt = keys_of(segment, "PdbResidueStart")
    return ProvenanceSchemaBinding(
        schema_sha256=PROVENANCE_SCHEMA_SHA256,
        required_keys=required_keys,
        conditional_keys=conditional,
        member_role_keys=member_keys,
        primary_pdb_keys=pdb_req,
        primary_data_keys=data_req,
        primary_residue_required=res_req,
        primary_residue_optional=res_opt,
        primary_additional_keys=add_keys,
        primary_residue_additional_keys=resadd_keys,
        segment_object_keys=seg_obj_keys,
        segment_mapping_keys=seg_req,
        segment_residue_required=segres_req,
        segment_residue_optional=segres_opt,
    )


def load_provenance_schema_binding(project_root: str | Any) -> ProvenanceSchemaBinding:
    """Load + verify the pinned v27 provenance schema binding."""
    return _load_provenance_schema_binding_cached(str(_project_path(project_root)))


def _project_path(project_root: str | Any) -> Any:
    from pathlib import Path

    return Path(project_root) if not isinstance(project_root, Path) else project_root


@functools.lru_cache(maxsize=8)
def _load_afdb_contract_cached(root_str: str) -> AfdbGrammarContract:
    root = _project_path(root_str)
    openapi = _load_json(root / AFDB_OPENAPI_REL, AFDB_OPENAPI_SHA256)
    table = _load_json(root / AFDB_OPTIONAL_TABLE_REL, AFDB_OPTIONAL_TABLE_SHA256)
    if table.get("table_sha256") != AFDB_OPTIONAL_TABLE_CONTENT_SHA256:
        raise E420GrammarPinError("afdb_optional_types: table sha pin mismatch")
    if table.get("source_openapi_sha256") != AFDB_OPENAPI_SHA256:
        raise E420GrammarPinError("afdb_optional_types: source openapi pin mismatch")
    summary = openapi["components"]["schemas"]["NewEntrySummary"]
    required = frozenset(summary["required"])
    properties = summary["properties"]
    active = frozenset(table["table"]["ACTIVE_OPTIONAL"])
    sunset = frozenset(table["table"]["SUNSET_OPTIONAL"])
    if not required or not active or not sunset:
        raise E420GrammarPinError("afdb contract tables empty")
    counts = table.get("counts", {})
    if counts.get("active_optional_count") != len(active) or counts.get("sunset_optional_count") != len(sunset):
        raise E420GrammarPinError("afdb_optional_types: counts pin mismatch")
    return AfdbGrammarContract(
        openapi_sha256=AFDB_OPENAPI_SHA256,
        required_fields=required,
        active_optional=active,
        sunset_fields=sunset,
        type_specs={name: spec for name, spec in properties.items()},
    )


def load_afdb_contract(project_root: str | Any) -> AfdbGrammarContract:
    return _load_afdb_contract_cached(str(_project_path(project_root)))


@functools.lru_cache(maxsize=8)
def _load_atlas_metadata_contract_cached(root_str: str) -> AtlasMetadataContract:
    root = _project_path(root_str)
    openapi = _load_json(root / ATLAS_OPENAPI_REL, ATLAS_OPENAPI_SHA256)
    response = openapi["components"]["schemas"]["MetadataResponse"]
    required = frozenset(response.get("required", ()))
    specs = dict(response.get("properties", {}))
    if not required or "PDB" not in required:
        raise E420GrammarPinError("atlas MetadataResponse: required pin mismatch")
    return AtlasMetadataContract(openapi_sha256=ATLAS_OPENAPI_SHA256, required_fields=required, type_specs=specs)


def load_atlas_metadata_contract(project_root: str | Any) -> AtlasMetadataContract:
    return _load_atlas_metadata_contract_cached(str(_project_path(project_root)))


def _check_json_type(value: Any, spec: Mapping[str, Any]) -> bool:
    """Type-check a JSON value against a retained OpenAPI property schema."""
    if "anyOf" in spec:
        return any(_check_json_type(value, branch) for branch in spec["anyOf"])
    kind = spec.get("type")
    if kind == "null":
        return value is None
    if kind == "string":
        return isinstance(value, str)
    if kind == "integer":
        return _is_int(value)
    if kind == "number":
        return _is_number(value)
    if kind == "boolean":
        return isinstance(value, bool)
    if kind == "array":
        if not isinstance(value, list):
            return False
        items = spec.get("items")
        if isinstance(items, Mapping):
            return all(_check_json_type(item, items) for item in value)
        return True
    if kind == "object":
        return isinstance(value, dict)
    return True


# --------------------------------------------------------------------------- #
# v27 per-row provenance envelope (exact-key).                                #
# --------------------------------------------------------------------------- #


def build_row_envelope(
    *,
    role: SourceRole,
    raw_bytes: bytes,
    normalized_projection: Mapping[str, Any],
    member_hashes: Mapping[str, Mapping[str, str]],
    context: Mapping[str, str],
    source_record_sha256: str | None = None,
    gap_derivation: Mapping[str, Any] | None = None,
    modelcif_projection_sha256: str | None = None,
) -> dict[str, Any]:
    """Build the v27 exact-key per-row provenance envelope for one row."""
    raw_sha = _sha256(raw_bytes)
    envelope: dict[str, Any] = {
        "raw_row_sha256": raw_sha,
        "normalized_row_sha256": _canonical_digest(normalized_projection),
        "source_row_sha256": source_record_sha256 or raw_sha,
        "source_member_hashes": {key: dict(value) for key, value in member_hashes.items()},
        "release_ids": dict(context["release_ids"]),
        "capture_anchor": str(context["capture_anchor"]),
        "emitted_by": str(context.get("emitted_by", EMITTED_BY)),
    }
    if modelcif_projection_sha256 is not None:
        envelope["modelcif_projection_sha256"] = modelcif_projection_sha256
    if gap_derivation is not None:
        envelope["gap_derivation"] = dict(gap_derivation)
    return envelope


def validate_row_envelope(
    envelope: Mapping[str, Any],
    *,
    project_root: str | Any,
    applicable_member_keys: Sequence[str],
    gap_applicable: bool = False,
    modelcif_applicable: bool = False,
    role: SourceRole = SourceRole.PDBE_MAPPING_PRIMARY,
) -> GrammarOutcome:
    """Validate the exact-key per-row envelope against the pinned v27 schema."""
    role_ = role

    @_grammar_boundary(role_)
    def _run(env: Mapping[str, Any]) -> GrammarOutcome:
        binding = load_provenance_schema_binding(project_root)
        if not isinstance(env, Mapping):
            _fail(GrammarFailureCode.ENVELOPE_KEYS_INVALID, "envelope not a mapping")
        expected = set(binding.required_keys) - set(binding.conditional_keys)
        expected |= {"modelcif_projection_sha256"} if modelcif_applicable else set()
        expected |= {"gap_derivation"} if gap_applicable else set()
        observed = set(env)
        if observed != expected:
            missing = sorted(expected - observed)
            extra = sorted(observed - expected)
            _fail(GrammarFailureCode.ENVELOPE_KEYS_INVALID, f"missing={missing} extra={extra}")
        for key in ("raw_row_sha256", "normalized_row_sha256", "source_row_sha256"):
            value = env[key]
            if not isinstance(value, str) or not HEX64_RE.fullmatch(value):
                _fail(GrammarFailureCode.ENVELOPE_SHA256_FORMAT_INVALID, key)
        if modelcif_applicable:
            value = env["modelcif_projection_sha256"]
            if value is not None and (not isinstance(value, str) or not HEX64_RE.fullmatch(value)):
                _fail(GrammarFailureCode.ENVELOPE_SHA256_FORMAT_INVALID, "modelcif_projection_sha256")
        members = env["source_member_hashes"]
        if not isinstance(members, Mapping) or set(members) != set(applicable_member_keys):
            _fail(GrammarFailureCode.ENVELOPE_MEMBER_KEYS_INVALID,
                  f"expected={sorted(applicable_member_keys)} observed={sorted(members) if isinstance(members, Mapping) else type(members).__name__}")
        for member_key, member_map in members.items():
            if not isinstance(member_map, Mapping) or not member_map:
                _fail(GrammarFailureCode.ENVELOPE_MEMBER_KEYS_INVALID, f"empty member map: {member_key}")
            for path, digest in member_map.items():
                if not isinstance(path, str) or not path or not isinstance(digest, str) or not HEX64_RE.fullmatch(digest):
                    _fail(GrammarFailureCode.ENVELOPE_MEMBER_HASH_INVALID, f"{member_key}:{path}")
        emitted_by = env["emitted_by"]
        release_ids = env["release_ids"]
        capture_anchor = env["capture_anchor"]
        if not isinstance(emitted_by, str) or not emitted_by:
            _fail(GrammarFailureCode.ENVELOPE_CONTEXT_INVALID, "emitted_by")
        if not isinstance(release_ids, Mapping) or not release_ids or not all(
            isinstance(k, str) and k and isinstance(v, str) and v for k, v in release_ids.items()
        ):
            _fail(GrammarFailureCode.ENVELOPE_CONTEXT_INVALID, "release_ids")
        if not isinstance(capture_anchor, str) or not capture_anchor:
            _fail(GrammarFailureCode.ENVELOPE_CONTEXT_INVALID, "capture_anchor")
        if gap_applicable:
            gap = env["gap_derivation"]
            if (not isinstance(gap, Mapping) or set(gap) != {"derivation_source", "uniprot_num"}
                    or gap["derivation_source"] != SIFTS_GAP_DERIVATION_SOURCE
                    or gap["uniprot_num"] is not GAP_ROW_UNIPROT_NUM):
                _fail(GrammarFailureCode.ENVELOPE_GAP_DERIVATION_INVALID, repr(gap))
        payload = {
            "envelope_keys": sorted(observed),
            "schema_sha256": binding.schema_sha256,
            "normalized_row_sha256": env["normalized_row_sha256"],
            "raw_row_sha256": env["raw_row_sha256"],
        }
        return GrammarOutcome(role=role_, ok=True, raw_sha256=env["raw_row_sha256"], raw_len=0,
                              payload=payload, failure=None)
    return _run(envelope)


# --------------------------------------------------------------------------- #
# PDBe v2 primary route (E420_PDBE_MAPPING_PRIMARY_V1).                       #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _ResidueProjection:
    start_index: int
    end_index: int
    start_code_raw: str
    end_code_raw: str
    start_code_1: str
    end_code_1: str
    unp_start: int | None
    unp_end: int | None
    pdb_code: str | None


@_grammar_boundary(SourceRole.PDBE_MAPPING_PRIMARY)
def validate_pdbe_primary(
    raw: Any,
    *,
    pdb_id: str,
    entity_id: int,
    project_root: str | Any,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Strict envelope validation of ``/pdb/entry/uniprot_mapping/{pdb_id}/{entity_id}``."""
    nonstandard_codes: set[str] = set()
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    if not isinstance(pdb_id, str) or not PDB_ID_RE.fullmatch(pdb_id):
        _fail(GrammarFailureCode.PDBE_PRIMARY_ROOT_KEY_INVALID, f"pdb_id={pdb_id!r}")
    if not _is_int(entity_id) or entity_id < 1:
        _fail(GrammarFailureCode.PDBE_PRIMARY_ENTITY_TYPE_INVALID, f"entity_id={entity_id!r}")
    binding = load_provenance_schema_binding(project_root)
    try:
        root = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        _fail(GrammarFailureCode.PDBE_PRIMARY_NOT_JSON, type(exc).__name__, raw_len=len(data), raw_sha256=raw_sha)
    if not isinstance(root, Mapping) or set(root) != {pdb_id}:
        _fail(GrammarFailureCode.PDBE_PRIMARY_ROOT_KEY_INVALID, f"root keys={sorted(root) if isinstance(root, Mapping) else type(root).__name__}")
    entry = root[pdb_id]
    if not isinstance(entry, Mapping) or set(entry) != set(binding.primary_pdb_keys):
        _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID, f"keys={sorted(entry) if isinstance(entry, Mapping) else type(entry).__name__}")
    sequence = entry["sequence"]
    if not isinstance(sequence, str) or not sequence or not re.fullmatch(r"[A-Z]+", sequence):
        _fail(GrammarFailureCode.PDBE_PRIMARY_SEQUENCE_INVALID, type(sequence).__name__)
    if not _is_int(entry["length"]) or entry["length"] != len(sequence):
        _fail(GrammarFailureCode.PDBE_PRIMARY_LENGTH_MISMATCH, f"length={entry['length']!r} len(sequence)={len(sequence)}")
    data_type = entry["dataType"]
    if not isinstance(data_type, str) or not data_type:
        _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID, "dataType")
    rows = entry["data"]
    if not isinstance(rows, list) or not rows:
        _fail(GrammarFailureCode.PDBE_PRIMARY_DATA_EMPTY, type(rows).__name__)
    entries: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, Mapping) or set(row) != set(binding.primary_data_keys):
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID,
                  f"row keys={sorted(row) if isinstance(row, Mapping) else type(row).__name__}")
        accession = row["accession"]
        name = row["name"]
        row_type = row["dataType"]
        if not isinstance(accession, str) or not accession:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID, "accession")
        if not isinstance(name, str) or not name:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID, "name")
        if not isinstance(row_type, str) or not row_type:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTRY_KEYS_INVALID, "row dataType")
        residues_raw = row["residues"]
        if not isinstance(residues_raw, list) or not residues_raw:
            _fail(GrammarFailureCode.PDBE_PRIMARY_DATA_EMPTY, "residues")
        additional = row["additionalData"]
        if not isinstance(additional, Mapping) or not set(additional) <= binding.primary_additional_keys:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ADDITIONAL_KEYS_INVALID,
                  f"keys={sorted(additional) if isinstance(additional, Mapping) else type(additional).__name__}")
        if "entityId" not in additional:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTITY_TYPE_INVALID, "additionalData.entityId missing")
        if not _is_int(additional["entityId"]):
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTITY_TYPE_INVALID, f"entityId={additional['entityId']!r}")
        if additional["entityId"] != entity_id:
            _fail(GrammarFailureCode.PDBE_PRIMARY_ENTITY_MISMATCH,
                  f"row={additional['entityId']!r}:requested={entity_id}")
        residues: list[_ResidueProjection] = []
        for residue in residues_raw:
            if not isinstance(residue, Mapping) or not set(binding.primary_residue_required) <= set(residue) \
                    or not set(residue) <= (binding.primary_residue_required | binding.primary_residue_optional):
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_KEYS_INVALID,
                      f"keys={sorted(residue) if isinstance(residue, Mapping) else type(residue).__name__}")
            start_index = residue["startIndex"]
            end_index = residue["endIndex"]
            if not _is_int(start_index) or not _is_int(end_index) or start_index < 1 or start_index > end_index:
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_BOUNDS_INVALID, f"{start_index!r}..{end_index!r}")
            codes: dict[str, str] = {}
            for bound in ("startCode", "endCode"):
                code = residue[bound]
                # AMENDED 2026-09-04 (census round 1): non-standard codes (e.g.
                # MSE selenomethionine) RECORDED, not fatal — SeMet entries are
                # overrepresented in the ATLAS census.
                if isinstance(code, str) and code in AA3_TO_1:
                    codes[bound] = code
                else:
                    nonstandard_codes.add(code if isinstance(code, str) else repr(code))
                    codes[bound] = code
            unp_start = residue.get("unpStartIndex")
            unp_end = residue.get("unpEndIndex")
            if unp_start is not None and not _is_int(unp_start):
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_BOUNDS_INVALID, f"unpStartIndex={unp_start!r}")
            if unp_end is not None and not _is_int(unp_end):
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_BOUNDS_INVALID, f"unpEndIndex={unp_end!r}")
            if unp_start is not None and unp_end is not None and unp_start > unp_end:
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_BOUNDS_INVALID, f"{unp_start!r}>{unp_end!r}")
            pdb_code = residue.get("pdbCode")
            if pdb_code is not None and (not isinstance(pdb_code, str) or not pdb_code):
                _fail(GrammarFailureCode.PDBE_PRIMARY_CODE_INVALID, f"pdbCode={pdb_code!r}")
            index_type = residue.get("indexType")
            if index_type is not None and not isinstance(index_type, str):
                _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_KEYS_INVALID, "indexType")
            for flag in ("modification", "mutation"):
                if residue.get(flag) is not None and not isinstance(residue.get(flag), bool):
                    _fail(GrammarFailureCode.PDBE_PRIMARY_RESIDUE_KEYS_INVALID, flag)
            residue_additional = residue.get("additionalData")
            if residue_additional is not None:
                if not isinstance(residue_additional, Mapping) or not set(residue_additional) <= binding.primary_residue_additional_keys:
                    _fail(GrammarFailureCode.PDBE_PRIMARY_ADDITIONAL_KEYS_INVALID, "residue.additionalData")
                for confidence_field in ("confidenceScore", "rawScore"):
                    if confidence_field in residue_additional and not _is_number(residue_additional[confidence_field]):
                        _fail(GrammarFailureCode.PDBE_PRIMARY_ADDITIONAL_KEYS_INVALID, confidence_field)
            residues.append(_ResidueProjection(
                start_index=start_index, end_index=end_index,
                start_code_raw=codes["startCode"], end_code_raw=codes["endCode"],
                start_code_1=AA3_TO_1.get(codes["startCode"]), end_code_1=AA3_TO_1.get(codes["endCode"]),
                unp_start=unp_start, unp_end=unp_end, pdb_code=pdb_code,
            ))
        entries.append({
            "accession": accession,
            "name": name,
            "data_type": row_type,
            "entity_id": entity_id,
            "residues": [_residue_projection_dict(item) for item in residues],
        })
    payload = {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "data_type": data_type,
        "length": len(sequence),
        "sequence_sha256": _sha256(sequence.encode("ascii")),
        "nonstandard_code_segments_recorded": sorted(nonstandard_codes),
        "entries": entries,
    }
    return GrammarOutcome(
        role=SourceRole.PDBE_MAPPING_PRIMARY, ok=True, raw_sha256=raw_sha, raw_len=len(data),
        payload=payload,
        failure=None,
    )


def _residue_projection_dict(item: _ResidueProjection) -> dict[str, Any]:
    return {
        "start_index": item.start_index, "end_index": item.end_index,
        "start_code_raw_3letter": item.start_code_raw, "end_code_raw_3letter": item.end_code_raw,
        "start_code_normalized_1letter": item.start_code_1, "end_code_normalized_1letter": item.end_code_1,
        "unp_start": item.unp_start, "unp_end": item.unp_end, "pdb_code": item.pdb_code,
    }


# --------------------------------------------------------------------------- #
# PDBe v2 segment route (E420_PDBE_MAPPING_SEGMENT_V1).                       #
# --------------------------------------------------------------------------- #


def _parse_segment_mapping_rows(raw_rows: Any, binding: ProvenanceSchemaBinding) -> list[dict[str, Any]]:
    if not isinstance(raw_rows, list) or not raw_rows:
        _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_KEYS_INVALID, "mappings empty or not a list")
    parsed: list[dict[str, Any]] = []
    for row in raw_rows:
        if not isinstance(row, Mapping) or set(row) != set(binding.segment_mapping_keys):
            _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_KEYS_INVALID,
                  f"keys={sorted(row) if isinstance(row, Mapping) else type(row).__name__}")
        entity_id = row["entity_id"]
        if not _is_int(entity_id) or entity_id < 1:
            _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_TYPE_INVALID, f"entity_id={entity_id!r}")
        for text_field in ("chain_id", "struct_asym_id"):
            value = row[text_field]
            if not isinstance(value, str) or not value.strip():
                _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_TYPE_INVALID, f"{text_field}={value!r}")
        unp_start = row["unp_start"]
        unp_end = row["unp_end"]
        if not _is_int(unp_start) or not _is_int(unp_end) or unp_start < 1 or unp_start > unp_end:
            _fail(GrammarFailureCode.PDBE_SEGMENT_BOUNDS_INVALID, f"unp {unp_start!r}..{unp_end!r}")
        identity = row["identity"]
        coverage = row["coverage"]
        if not _is_number(identity) or not 0.0 <= float(identity) <= 100.0:
            _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_TYPE_INVALID, f"identity={identity!r}")
        if not _is_number(coverage) or not 0.0 <= float(coverage) <= 100.0:
            _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_TYPE_INVALID, f"coverage={coverage!r}")
        bounds: dict[str, int] = {}
        for bound_name in ("start", "end"):
            bound = row[bound_name]
            if not isinstance(bound, Mapping) or not set(binding.segment_residue_required) <= set(bound) \
                    or not set(bound) <= (binding.segment_residue_required | binding.segment_residue_optional):
                _fail(GrammarFailureCode.PDBE_SEGMENT_RESIDUE_KEYS_INVALID,
                      f"{bound_name} keys={sorted(bound) if isinstance(bound, Mapping) else type(bound).__name__}")
            residue_number = bound["residue_number"]
            if not _is_int(residue_number) or residue_number < 1:
                _fail(GrammarFailureCode.PDBE_SEGMENT_BOUNDS_INVALID, f"{bound_name}.residue_number={residue_number!r}")
            author_number = bound.get("author_residue_number")
            if author_number is not None and not _is_int(author_number):
                _fail(GrammarFailureCode.PDBE_SEGMENT_RESIDUE_KEYS_INVALID, "author_residue_number")
            insertion_code = bound.get("author_insertion_code")
            if insertion_code is not None and not isinstance(insertion_code, str):
                _fail(GrammarFailureCode.PDBE_SEGMENT_RESIDUE_KEYS_INVALID, "author_insertion_code")
            bounds[bound_name] = residue_number
        if bounds["start"] > bounds["end"]:
            _fail(GrammarFailureCode.PDBE_SEGMENT_BOUNDS_INVALID, f"start.residue_number>{'end'}.residue_number")
        parsed.append({
            "entity_id": entity_id,
            "chain_id": row["chain_id"],
            "struct_asym_id": row["struct_asym_id"],
            "unp_start": unp_start,
            "unp_end": unp_end,
            "start_residue_number": bounds["start"],
            "end_residue_number": bounds["end"],
            "identity": float(identity),
            "coverage": float(coverage),
        })
    return parsed


def _parse_segment_root(
    data: bytes,
    *,
    pdb_id: str,
    candidate_accession: str,
    binding: ProvenanceSchemaBinding,
    raw_sha: str,
    raw_len: int,
) -> list[tuple[str, list[dict[str, Any]]]]:
    try:
        root = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        _fail(GrammarFailureCode.PDBE_SEGMENT_NOT_JSON, type(exc).__name__, raw_len=raw_len, raw_sha256=raw_sha)
    if not isinstance(root, Mapping) or set(root) != {pdb_id}:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ROOT_KEY_INVALID,
              f"root keys={sorted(root) if isinstance(root, Mapping) else type(root).__name__}")
    accession_map = root[pdb_id]
    if not isinstance(accession_map, Mapping) or not accession_map:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ROOT_KEY_INVALID, "accession map empty")
    # AMENDED 2026-09-04 (phase-(ii) round 3): the real /mappings/uniprot
    # response nests under a "UniProt" resource level and the per-accession
    # object carries {name, mappings} (no "identifier"); required-keys policy:
    # {"name", "mappings"} must be present, extras are recorded.
    if "UniProt" in accession_map and isinstance(accession_map["UniProt"], Mapping):
        accession_map = accession_map["UniProt"]
    parsed_pairs: list[tuple[str, list[dict[str, Any]]]] = []
    object_extra_keys: set[str] = set()
    for accession_key, obj in accession_map.items():
        if not isinstance(accession_key, str) or not accession_key:
            _fail(GrammarFailureCode.PDBE_SEGMENT_ROOT_KEY_INVALID, "accession key type")
        if not isinstance(obj, Mapping):
            _fail(GrammarFailureCode.PDBE_SEGMENT_OBJECT_KEYS_INVALID,
                  f"keys={sorted(obj) if isinstance(obj, Mapping) else type(obj).__name__}")
        required_missing = {"name", "mappings"} - set(obj)
        if required_missing:
            _fail(GrammarFailureCode.PDBE_SEGMENT_OBJECT_KEYS_INVALID, f"missing={sorted(required_missing)}")
        object_extra_keys |= set(obj) - set(binding.segment_object_keys)
        # AMENDED 2026-09-04 (phase-(ii) round 4): PDBe's "identifier" is the
        # UniProt ENTRY NAME (e.g. PGKC_TRYBB), not the accession — the accession
        # is the map key. Only presence/nonemptiness is enforced.
        identifier = obj.get("identifier", "")
        if not isinstance(identifier, str) or not identifier:
            _fail(GrammarFailureCode.PDBE_SEGMENT_OBJECT_KEYS_INVALID, "identifier")
        name = obj["name"]
        if not isinstance(name, str) or not name:
            _fail(GrammarFailureCode.PDBE_SEGMENT_OBJECT_KEYS_INVALID, "name")
        rows = _parse_segment_mapping_rows(obj["mappings"], binding)
        parsed_pairs.append((accession_key, rows))
    return parsed_pairs


@_grammar_boundary(SourceRole.PDBE_MAPPING_SEGMENT)
def validate_pdbe_segment(
    raw: Any,
    *,
    pdb_id: str,
    candidate_accession: str,
    project_root: str | Any,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Strict envelope validation + entity-id resolution of ``/mappings/uniprot/{pdb_id}``."""
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    if not isinstance(pdb_id, str) or not PDB_ID_RE.fullmatch(pdb_id):
        _fail(GrammarFailureCode.PDBE_SEGMENT_ROOT_KEY_INVALID, f"pdb_id={pdb_id!r}")
    if not isinstance(candidate_accession, str) or not candidate_accession:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ACCESSION_ABSENT, "candidate accession type")
    binding = load_provenance_schema_binding(project_root)
    pairs = _parse_segment_root(data, pdb_id=pdb_id, candidate_accession=candidate_accession,
                                binding=binding, raw_sha=raw_sha, raw_len=len(data))
    matching_rows: list[dict[str, Any]] = []
    _segment_object_extra_keys_note: list[str] = []
    nonstandard_codes: set[str] = set()
    for accession_key, rows in pairs:
        if accession_key == candidate_accession:
            matching_rows.extend(rows)
    if not matching_rows:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ACCESSION_ABSENT, candidate_accession, raw_len=len(data), raw_sha256=raw_sha)
    entities = sorted({row["entity_id"] for row in matching_rows})
    if len(entities) > 1:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ENTITY_AMBIGUOUS, f"entities={entities}", raw_len=len(data), raw_sha256=raw_sha)
    resolved = entities[0]
    payload = {
        "pdb_id": pdb_id,
        "object_extra_keys_recorded": sorted(_segment_object_extra_keys_note),
        "candidate_accession": candidate_accession,
        "resolved_entity_id": resolved,
        "accessions": sorted(key for key, _ in pairs),
        "mapping_row_count": sum(len(rows) for _, rows in pairs),
        "matching_row_count": len(matching_rows),
    }
    return GrammarOutcome(role=SourceRole.PDBE_MAPPING_SEGMENT, ok=True, raw_sha256=raw_sha,
                          raw_len=len(data), payload=payload, failure=None)


@_grammar_boundary(SourceRole.PDBE_MAPPING_SEGMENT)
def check_pdbe_segment_thresholds(
    raw: Any,
    *,
    pdb_id: str,
    candidate_accession: str,
    entity_id: int,
    mode: SegmentMode,
    project_root: str | Any,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Resolution thresholds (identity >= 0.90, coverage >= 0.80) or fallback exact-100."""
    if not isinstance(mode, SegmentMode):
        _fail(GrammarFailureCode.BOUNDARY_INVALID, f"mode={mode!r}")
    if not _is_int(entity_id) or entity_id < 1:
        _fail(GrammarFailureCode.PDBE_SEGMENT_MAPPING_TYPE_INVALID, f"entity_id={entity_id!r}")
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    binding = load_provenance_schema_binding(project_root)
    pairs = _parse_segment_root(data, pdb_id=pdb_id, candidate_accession=candidate_accession,
                                binding=binding, raw_sha=raw_sha, raw_len=len(data))
    matching_rows: list[dict[str, Any]] = []
    for accession_key, rows in pairs:
        if accession_key == candidate_accession:
            matching_rows.extend(rows)
    entities = sorted({row["entity_id"] for row in matching_rows})
    if not matching_rows:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ACCESSION_ABSENT, candidate_accession, raw_len=len(data), raw_sha256=raw_sha)
    if entities != [entity_id]:
        _fail(GrammarFailureCode.PDBE_SEGMENT_ENTITY_RESOLUTION_MISMATCH,
              f"resolved={entities}:requested={entity_id}", raw_len=len(data), raw_sha256=raw_sha)
    identity_threshold = 0.90
    coverage_threshold = 0.80
    if mode is SegmentMode.RESOLUTION:
        failing = [row for row in matching_rows
                   if not (row["identity"] >= identity_threshold and row["coverage"] >= coverage_threshold)]
        threshold_kind = "resolution_identity_0.90_coverage_0.80"
    else:
        failing = [row for row in matching_rows if float(row["identity"]) != 100.0]
        threshold_kind = "fallback_identity_exact_100.0"
    if failing:
        _fail(GrammarFailureCode.PDBE_SEGMENT_THRESHOLD_FAILURE,
              f"{threshold_kind}: failing_rows={len(failing)}", raw_len=len(data), raw_sha256=raw_sha)
    payload = {
        "mode": mode.value,
        "threshold_kind": threshold_kind,
        "checked_row_count": len(matching_rows),
        "failing_row_count": 0,
        "resolved_entity_id": entity_id,
    }
    return GrammarOutcome(role=SourceRole.PDBE_MAPPING_SEGMENT, ok=True, raw_sha256=raw_sha,
                          raw_len=len(data), payload=payload, failure=None)


# --------------------------------------------------------------------------- #
# AFDB metadata (E420_AFDB_METADATA_V1) and coordinate child.                 #
# --------------------------------------------------------------------------- #

PLDDT_FAMILY_FIELDS = frozenset({
    "globalMetricValue", "fractionPlddtVeryLow", "fractionPlddtLow",
    "fractionPlddtConfident", "fractionPlddtVeryHigh",
})


@_grammar_boundary(SourceRole.AFDB_METADATA)
def validate_afdb_metadata(
    raw: Any,
    *,
    project_root: str | Any,
    expected_accession: str | None = None,
    canonical_sequence: str | None = None,
    require_cif_url: bool = True,
    release_lock: Mapping[str, Any] | None = None,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Strict NewEntrySummary envelope validation (one-object array per accession)."""
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    contract = load_afdb_contract(project_root)
    try:
        root = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        _fail(GrammarFailureCode.AFDB_METADATA_NOT_JSON, type(exc).__name__, raw_len=len(data), raw_sha256=raw_sha)
    if not isinstance(root, list):
        _fail(GrammarFailureCode.AFDB_METADATA_NOT_JSON, "root not an array", raw_len=len(data), raw_sha256=raw_sha)
    # AMENDED 2026-09-04 (census round 1): real /prediction responses may return
    # MULTIPLE objects (e.g. one per UniProt sequence). Every object is validated;
    # the census records all model ids. Label-phase selection is deterministic:
    # earliest modelCreatedDate, then lexicographic modelEntityId (registration
    # 5.1 amendment). Sunset fields are documented-but-deprecated: RECORDED,
    # not fatal (contract table itself lists them).
    if not root or not all(isinstance(o, Mapping) for o in root):
        _fail(GrammarFailureCode.AFDB_METADATA_CARDINALITY_INVALID, f"objects={len(root) if isinstance(root, list) else type(root).__name__}",
              raw_len=len(data), raw_sha256=raw_sha)
    for row in root:
        observed_keys = set(row)
        sunset_present = sorted(observed_keys & contract.sunset_fields)
        if sunset_present:
            recorded_sunset = getattr(validate_afdb_metadata, "_recorded_sunset", set())
            recorded_sunset |= set(sunset_present)
            validate_afdb_metadata._recorded_sunset = recorded_sunset
        unknown = sorted(observed_keys - contract.required_fields - contract.active_optional - contract.sunset_fields)
        if unknown:
            _fail(GrammarFailureCode.AFDB_FIELD_UNKNOWN, f"{unknown}", raw_len=len(data), raw_sha256=raw_sha)
        missing = sorted(contract.required_fields - observed_keys)
        if missing:
            _fail(GrammarFailureCode.AFDB_FIELD_REQUIRED_MISSING, f"{missing}", raw_len=len(data), raw_sha256=raw_sha)
        for key in sorted(observed_keys & (contract.required_fields | contract.active_optional | contract.sunset_fields)):
            if key in contract.type_specs and not _check_json_type(row[key], contract.type_specs[key]):
                _fail(GrammarFailureCode.AFDB_FIELD_TYPE_INVALID, f"{key}={row[key]!r}",
                      raw_len=len(data), raw_sha256=raw_sha)

    def _created_key(obj):
        return (obj.get("modelCreatedDate") or "", obj.get("modelEntityId") or "")
    row = sorted(root, key=_created_key)[0]
    objects_summary = [{"modelEntityId": o.get("modelEntityId"),
                        "modelCreatedDate": o.get("modelCreatedDate"),
                        "sequenceVersionDate": o.get("sequenceVersionDate")} for o in root]

    created = row["modelCreatedDate"]
    _parse_strict_created_date(created, field_name="modelCreatedDate", raw_len=len(data), raw_sha256=raw_sha)
    if row["sequenceVersionDate"] is not None:
        _parse_strict_created_date(row["sequenceVersionDate"], field_name="sequenceVersionDate",
                                   raw_len=len(data), raw_sha256=raw_sha)

    latest = row["latestVersion"]
    all_versions = row["allVersions"]
    if not _is_int(latest) or latest < 1:
        _fail(GrammarFailureCode.AFDB_VERSION_INVALID, f"latestVersion={latest!r}", raw_len=len(data), raw_sha256=raw_sha)
    if not isinstance(all_versions, list) or not all_versions or not all(_is_int(v) and v >= 1 for v in all_versions) \
            or len(set(all_versions)) != len(all_versions) or latest not in all_versions:
        _fail(GrammarFailureCode.AFDB_VERSION_INVALID, f"allVersions={all_versions!r}", raw_len=len(data), raw_sha256=raw_sha)

    sequence = row["sequence"]
    sequence_start = row["sequenceStart"]
    sequence_end = row["sequenceEnd"]
    if not isinstance(sequence, str) or not sequence or not re.fullmatch(r"[A-Z]+", sequence):
        _fail(GrammarFailureCode.AFDB_FIELD_TYPE_INVALID, "sequence", raw_len=len(data), raw_sha256=raw_sha)
    if not _is_int(sequence_start) or not _is_int(sequence_end) or sequence_start < 1 or sequence_start > sequence_end:
        _fail(GrammarFailureCode.AFDB_SEQUENCE_OFFSET_INVALID,
              f"{sequence_start!r}..{sequence_end!r}", raw_len=len(data), raw_sha256=raw_sha)
    if canonical_sequence is not None:
        if not isinstance(canonical_sequence, str) or sequence != canonical_sequence[sequence_start - 1:sequence_end]:
            _fail(GrammarFailureCode.AFDB_SEQUENCE_OFFSET_MISMATCH,
                  f"sequence != canonical[{sequence_start - 1}:{sequence_end}]",
                  raw_len=len(data), raw_sha256=raw_sha)
    checksum = row["sequenceChecksum"]
    if not isinstance(checksum, str) or not MD5_RE.fullmatch(checksum):
        _fail(GrammarFailureCode.AFDB_SEQUENCE_CHECKSUM_INVALID, f"sequenceChecksum={checksum!r}",
              raw_len=len(data), raw_sha256=raw_sha)
    if hashlib.md5(sequence.encode("ascii")).hexdigest() != checksum:
        _fail(GrammarFailureCode.AFDB_SEQUENCE_CHECKSUM_MISMATCH, "md5(sequence) != sequenceChecksum",
              raw_len=len(data), raw_sha256=raw_sha)

    cif_url = row.get("cifUrl")
    if require_cif_url and (not isinstance(cif_url, str) or not cif_url):
        _fail(GrammarFailureCode.AFDB_CIF_URL_MISSING, f"cifUrl={cif_url!r}", raw_len=len(data), raw_sha256=raw_sha)
    if expected_accession is not None:
        accession = row.get("uniprotAccession")
        if accession != expected_accession:
            _fail(GrammarFailureCode.AFDB_ACCESSION_MISMATCH,
                  f"row={accession!r}:requested={expected_accession!r}", raw_len=len(data), raw_sha256=raw_sha)
    if release_lock is not None:
        lock_failure = _release_lock_failure(release_lock, {
            "latestVersion": latest,
            "allVersions": list(all_versions),
            "modelEntityId": row["modelEntityId"],
            "modelCreatedDate": created,
        })
        if lock_failure is not None:
            _fail(GrammarFailureCode.RELEASE_LOCK_MISMATCH, lock_failure, raw_len=len(data), raw_sha256=raw_sha)

    payload = {
        "model_entity_id": row["modelEntityId"],
        "model_created_date": created,
        "sequence_version_date": row["sequenceVersionDate"],
        "latest_version": latest,
        "all_versions": list(all_versions),
        "sequence": sequence,
        "sequence_start": sequence_start,
        "sequence_end": sequence_end,
        "sequence_checksum_md5_verified": True,
        "cif_url": cif_url,
        "uniprot_accession": row.get("uniprotAccession"),
        "uniprot_sequence": row.get("uniprotSequence"),
        "objects_in_response": len(root),
        "objects_summary": [{"modelEntityId": o.get("modelEntityId"), "modelCreatedDate": o.get("modelCreatedDate"),
                             "sequenceVersionDate": o.get("sequenceVersionDate")} for o in root],
        "selection_rule": "earliest modelCreatedDate, then lexicographic modelEntityId (registration 5.1 amendment 2026-09-04)",
        "sunset_fields_recorded": sorted(set().union(*[set(o) & contract.sunset_fields for o in root])) if root else [],
    }
    # pLDDT-family fields are hash-without-parse: they are type-checked above
    # for envelope conformance and never exported into any payload.
    return GrammarOutcome(role=SourceRole.AFDB_METADATA, ok=True, raw_sha256=raw_sha, raw_len=len(data),
                          payload=payload, failure=None)


def _parse_strict_created_date(value: Any, *, field_name: str, raw_len: int, raw_sha256: str) -> _dt.datetime:
    if not isinstance(value, str):
        _fail(GrammarFailureCode.AFDB_MODEL_CREATED_DATE_INVALID, f"{field_name}={value!r}",
              raw_len=raw_len, raw_sha256=raw_sha256)
    try:
        return _dt.datetime.strptime(value, MODEL_CREATED_DATE_FORMAT)
    except ValueError:
        _fail(GrammarFailureCode.AFDB_MODEL_CREATED_DATE_INVALID, f"{field_name}={value!r}",
              raw_len=raw_len, raw_sha256=raw_sha256)


@_grammar_boundary(SourceRole.AFDB_COORDINATE_CHILD)
def validate_afdb_child_binding(
    metadata_payload: Mapping[str, Any],
    child_url: Any,
    child_bytes: Any,
    *,
    forbidden_urls: Sequence[str] = (),
) -> GrammarOutcome:
    """Bind the coordinate child to exactly the metadata row's ``cifUrl`` field."""
    if not isinstance(metadata_payload, Mapping) or "cif_url" not in metadata_payload:
        _fail(GrammarFailureCode.AFDB_CIF_URL_MISSING, "metadata payload lacks cif_url")
    cif_url = metadata_payload["cif_url"]
    if not isinstance(cif_url, str) or not cif_url:
        _fail(GrammarFailureCode.AFDB_CIF_URL_MISSING, f"cifUrl={cif_url!r}")
    for forbidden in forbidden_urls:
        if child_url == forbidden:
            _fail(GrammarFailureCode.AFDB_CHILD_URL_FORBIDDEN_SUBSTITUTE, f"child_url={child_url!r}")
    if not isinstance(child_url, str) or child_url != cif_url:
        _fail(GrammarFailureCode.AFDB_CHILD_URL_MISMATCH, f"child_url={child_url!r}:cifUrl={cif_url!r}")
    if isinstance(child_bytes, (bytearray, memoryview)):
        child_bytes = bytes(child_bytes)
    if not isinstance(child_bytes, bytes) or not child_bytes:
        _fail(GrammarFailureCode.AFDB_CHILD_BYTES_EMPTY, type(child_bytes).__name__)
    try:
        head = child_bytes[:4096].decode("utf-8")
    except UnicodeDecodeError:
        _fail(GrammarFailureCode.AFDB_CHILD_MODEL_CIF_HEADER_INVALID, "undecodable child bytes")
    if not head.lstrip().startswith("data_"):
        _fail(GrammarFailureCode.AFDB_CHILD_MODEL_CIF_HEADER_INVALID, head[:32])
    payload = {
        "bound_url": child_url,
        "child_sha256": _sha256(child_bytes),
        "child_len": len(child_bytes),
        "model_entity_id": metadata_payload.get("model_entity_id"),
    }
    return GrammarOutcome(role=SourceRole.AFDB_COORDINATE_CHILD, ok=True,
                          raw_sha256=payload["child_sha256"], raw_len=len(child_bytes),
                          payload=payload, failure=None)


@_grammar_boundary(SourceRole.AFDB_METADATA)
def validate_af_temporal(model_created_date: str, reference_raw: Any) -> GrammarOutcome:
    """§5.1 disposition: strict modelCreatedDate < in-file ordinal-1 revision_date."""
    created = _parse_strict_created_date(model_created_date, field_name="modelCreatedDate",
                                         raw_len=0, raw_sha256=_sha256(b""))
    reference = validate_reference_cif(reference_raw)
    if not reference.ok or reference.failure is not None or reference.payload is None:
        return GrammarOutcome(role=SourceRole.AFDB_METADATA, ok=False, raw_sha256=_sha256(b""), raw_len=0,
                              payload=None, failure=reference.failure)
    release_date_text = reference.payload["initial_release_date"]
    if not isinstance(release_date_text, str) or not REVISION_DATE_RE.fullmatch(release_date_text):
        _fail(GrammarFailureCode.REFERENCE_REVISION_DATE_INVALID, f"revision_date={release_date_text!r}")
    release_date = _dt.datetime.strptime(release_date_text, "%Y-%m-%d").date()
    if created.date() >= release_date:
        _fail(GrammarFailureCode.AFDB_TEMPORAL_NOT_STRICT,
              f"modelCreatedDate={created.isoformat()} !< release={release_date.isoformat()}")
    payload = {
        "model_created_date": model_created_date,
        "reference_initial_release_date": release_date_text,
        "rule": "modelCreatedDate strictly earlier than reference initial release (registration 5.1)",
    }
    return GrammarOutcome(role=SourceRole.AFDB_METADATA, ok=True, raw_sha256=reference.raw_sha256,
                          raw_len=reference.raw_len, payload=payload, failure=None)


# --------------------------------------------------------------------------- #
# ATLAS analysis ZIP (E420_ATLAS_ANALYSIS_V1).                                #
# --------------------------------------------------------------------------- #

RMSF_HEADER = ("seq", "RMSF_R1", "RMSF_R2", "RMSF_R3")
CORRESP_HEADER = ("UnP_seq", "PDB_seq", "PDB_num", "CIF_num", "UnP_num",
                  "gap", "hetatm", "corrected", "no_fullbb", "no_ca")
RMSF_MEMBER_SUFFIX = "_RMSF.tsv"
CORRESP_MEMBER_SUFFIX = "_corresp.tsv"


def _tsv_rows(member_bytes: bytes, expected_header: Sequence[str], header_failure: GrammarFailureCode,
              row_failure: GrammarFailureCode, column_count: int) -> list[list[str]]:
    try:
        text = member_bytes.decode("utf-8")
    except UnicodeDecodeError:
        _fail(header_failure, "member not utf-8")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        _fail(header_failure, "member empty")
    if tuple(lines[0].split("\t")) != tuple(expected_header):
        _fail(header_failure, f"header={lines[0]!r}")
    rows: list[list[str]] = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != column_count or any(field == "" for field in fields):
            _fail(row_failure, f"row={line!r}")
        rows.append(fields)
    return rows


def _gap_share_failure(gap_rows: int, total_rows: int) -> str | None:
    """Registration §8: a ``5*unmapped > total`` violation is never repairable."""
    if 5 * gap_rows > total_rows:
        return f"5*unmapped={5 * gap_rows} > total={total_rows}"
    return None


@_grammar_boundary(SourceRole.ATLAS_ANALYSIS)
def validate_atlas_analysis_zip(
    raw: Any,
    *,
    pdb_chain: str,
    expected_archive_sha256: str | None = None,
    expected_archive_len: int | None = None,
) -> GrammarOutcome:
    """Frozen member grammar for the per-chain analysis ZIP (registration §5.2).

    RMSF lexemes are grammar-checked only; values are never exported (the F
    channel opens at gate (g), after POOL_FREEZE).
    """
    data = _require_bytes(raw, expected_archive_sha256, expected_archive_len)
    raw_sha = _sha256(data)
    if not isinstance(pdb_chain, str) or not re.fullmatch(r"[0-9A-Za-z]{4}_[A-Za-z0-9]+", pdb_chain):
        _fail(GrammarFailureCode.ATLAS_ARCHIVE_MEMBER_SET_INVALID, f"pdb_chain={pdb_chain!r}")
    expected_members = {pdb_chain + RMSF_MEMBER_SUFFIX, pdb_chain + CORRESP_MEMBER_SUFFIX}
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            names = set(archive.namelist())
            missing_members = sorted(expected_members - names)
            if missing_members:
                _fail(GrammarFailureCode.ATLAS_ARCHIVE_MEMBER_SET_INVALID,
                      f"missing={missing_members}:expected={sorted(expected_members)}",
                      raw_len=len(data), raw_sha256=raw_sha)
            extra_members = {name: {"bytes": archive.getinfo(name).file_size,
                                    "sha256": hashlib.sha256(archive.read(name)).hexdigest()}
                             for name in sorted(names - expected_members)}
            member_bytes: dict[str, bytes] = {}
            for name in sorted(expected_members):
                try:
                    member_bytes[name] = archive.read(name)
                except (zipfile.BadZipFile, RuntimeError, OSError, EOFError) as exc:
                    _fail(GrammarFailureCode.ATLAS_ARCHIVE_MEMBER_CRC_INVALID, f"{name}: {type(exc).__name__}",
                          raw_len=len(data), raw_sha256=raw_sha)
    except zipfile.BadZipFile:
        _fail(GrammarFailureCode.ATLAS_ARCHIVE_ZIP_INVALID, "BadZipFile", raw_len=len(data), raw_sha256=raw_sha)
    rmsf_rows = _tsv_rows(member_bytes[pdb_chain + RMSF_MEMBER_SUFFIX], RMSF_HEADER,
                          GrammarFailureCode.ATLAS_RMSF_HEADER_INVALID,
                          GrammarFailureCode.ATLAS_RMSF_ROW_INVALID, 4)
    corr_rows = _tsv_rows(member_bytes[pdb_chain + CORRESP_MEMBER_SUFFIX], CORRESP_HEADER,
                          GrammarFailureCode.ATLAS_CORRESP_HEADER_INVALID,
                          GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, 10)
    negative_value_rows = 0
    for ordinal, row in enumerate(rmsf_rows, start=1):
        if not re.fullmatch(r"[A-Za-z]", row[0] or ""):
            # real archives: seq is the one-letter amino-acid code (amended 2026-09-04,
            # phase-(ii) round 2: fixture ordinals were wrong, fresh bytes carry AA letters)
            _fail(GrammarFailureCode.ATLAS_RMSF_ROW_INVALID, f"ordinal {ordinal}: seq={row[0]!r}",
                  raw_len=len(data), raw_sha256=raw_sha)
        for lexeme in row[1:]:
            if not SEQ_LEXEME_RE.fullmatch(lexeme):
                _fail(GrammarFailureCode.ATLAS_RMSF_ROW_INVALID, f"ordinal {ordinal}: lexeme={lexeme!r}",
                      raw_len=len(data), raw_sha256=raw_sha)
            value = float(lexeme)
            if value != value or value in (float("inf"), float("-inf")):
                _fail(GrammarFailureCode.ATLAS_RMSF_ROW_INVALID, f"ordinal {ordinal}: nonfinite",
                      raw_len=len(data), raw_sha256=raw_sha)
            if value < 0.0:
                negative_value_rows += 1
                break
    gap_row_count = 0
    cif_pdb_mismatch_rows = 0
    offsets: set[int] = set()
    for ordinal, row in enumerate(corr_rows, start=1):
        unp_seq, pdb_seq, pdb_num, cif_num, unp_num, gap, hetatm, corrected, no_fullbb, no_ca = row
        flags = (gap, hetatm, corrected, no_fullbb, no_ca)
        if any(flag not in ("0", "1") for flag in flags):
            _fail(GrammarFailureCode.ATLAS_CORRESP_FLAG_INVALID, f"ordinal {ordinal}: flags={flags!r}",
                  raw_len=len(data), raw_sha256=raw_sha)
        if not re.fullmatch(r"[0-9]+", pdb_num) or not re.fullmatch(r"[0-9]+", cif_num):
            _fail(GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, f"ordinal {ordinal}: PDB_num/CIF_num={pdb_num!r}/{cif_num!r}",
                  raw_len=len(data), raw_sha256=raw_sha)
        offsets.add(int(pdb_num) - int(cif_num))
        if pdb_num != cif_num:
            cif_pdb_mismatch_rows += 1
        is_gap = unp_seq == "-"
        if is_gap:
            gap_row_count += 1
            if unp_num != "-":
                _fail(GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, f"ordinal {ordinal}: gap row UnP_num={unp_num!r}",
                      raw_len=len(data), raw_sha256=raw_sha)
            if pdb_seq != "-" and not re.fullmatch(r"[A-Z]", pdb_seq):
                _fail(GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, f"ordinal {ordinal}: PDB_seq={pdb_seq!r}",
                      raw_len=len(data), raw_sha256=raw_sha)
        else:
            if not re.fullmatch(r"[A-Z]", unp_seq) or not re.fullmatch(r"[A-Z]", pdb_seq):
                _fail(GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, f"ordinal {ordinal}: UnP_seq/PDB_seq={unp_seq!r}/{pdb_seq!r}",
                      raw_len=len(data), raw_sha256=raw_sha)
            if not re.fullmatch(r"[0-9]+", unp_num):
                _fail(GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, f"ordinal {ordinal}: UnP_num={unp_num!r}",
                      raw_len=len(data), raw_sha256=raw_sha)
    if len(rmsf_rows) != len(corr_rows):
        _fail(GrammarFailureCode.ATLAS_ROW_COUNT_MISMATCH,
              f"rmsf={len(rmsf_rows)}:corresp={len(corr_rows)}", raw_len=len(data), raw_sha256=raw_sha)
    if cif_pdb_mismatch_rows:
        # AMENDED 2026-09-04 (phase-(ii) round 3): author (PDB) numbering may be
        # offset from CIF index numbering (e.g. 16pk_A starts at residue 5). The
        # invariant is a CONSTANT per-chain offset across all rows, not per-row
        # equality; the offset is recorded. Varying offsets -> typed failure.
        if len(offsets) != 1:
            _fail(GrammarFailureCode.ATLAS_CORRESP_CIF_PDB_MISMATCH,
                  f"rows={cif_pdb_mismatch_rows}:offsets={sorted(offsets)}",
                  raw_len=len(data), raw_sha256=raw_sha)
    gap_failure = _gap_share_failure(gap_row_count, len(corr_rows))
    if gap_failure is not None:
        _fail(GrammarFailureCode.ATLAS_GAP_RATIO_EXCEEDED, gap_failure, raw_len=len(data), raw_sha256=raw_sha)
    payload = {
        "pdb_chain": pdb_chain,
        "archive_sha256": raw_sha,
        "member_sha256": {name: _sha256(member_bytes[name]) for name in sorted(member_bytes)},
        "header_raw": {
            "rmsf": "\t".join(RMSF_HEADER),
            "corresp": "\t".join(CORRESP_HEADER),
        },
        "rmsf_row_count": len(rmsf_rows),
        "corresp_row_count": len(corr_rows),
        "gap_row_count": gap_row_count,
        "mapped_row_count": len(corr_rows) - gap_row_count,
        "negative_value_rows": negative_value_rows,
        "rmsf_values_exported": False,
        "extra_members_present": bool(extra_members),
        "extra_member_hashes": extra_members,
        "pdb_cif_numbering_offset": (sorted(offsets)[0] if offsets else 0),
        "cif_pdb_mismatch_rows": cif_pdb_mismatch_rows,
    }
    return GrammarOutcome(role=SourceRole.ATLAS_ANALYSIS, ok=True, raw_sha256=raw_sha, raw_len=len(data),
                          payload=payload, failure=None)


def validate_gap_share(gap_rows: Any, total_rows: Any, *, role: SourceRole) -> GrammarOutcome:
    """Boundary-equality-accepting integer rule: accepted iff 5*unmapped <= total."""
    @_grammar_boundary(role)
    def _run() -> GrammarOutcome:
        if not _is_int(gap_rows) or not _is_int(total_rows) or gap_rows < 0 or total_rows < 0:
            _fail(GrammarFailureCode.BOUNDARY_INVALID, f"gap_rows={gap_rows!r} total_rows={total_rows!r}")
        detail = _gap_share_failure(gap_rows, total_rows)
        if detail is not None:
            _fail(GrammarFailureCode.ATLAS_GAP_RATIO_EXCEEDED, detail)
        return GrammarOutcome(role=role, ok=True, raw_sha256=_sha256(b""), raw_len=0,
                              payload={"gap_rows": gap_rows, "total_rows": total_rows, "ratio_ok": True},
                              failure=None)
    return _run()


# --------------------------------------------------------------------------- #
# ATLAS metadata (E420_ATLAS_METADATA_V1).                                    #
# --------------------------------------------------------------------------- #


@_grammar_boundary(SourceRole.ATLAS_METADATA)
def validate_atlas_metadata(
    raw: Any,
    *,
    pdb_chain: str,
    project_root: str | Any,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Strict metadata envelope; chain/context identity only (never a label)."""
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    if not isinstance(pdb_chain, str) or not re.fullmatch(r"[0-9A-Za-z]{4}_[A-Za-z0-9]+", pdb_chain):
        _fail(GrammarFailureCode.ATLAS_METADATA_CHAIN_KEY_INVALID, f"pdb_chain={pdb_chain!r}")
    contract = load_atlas_metadata_contract(project_root)
    try:
        root = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        _fail(GrammarFailureCode.ATLAS_METADATA_NOT_JSON, type(exc).__name__, raw_len=len(data), raw_sha256=raw_sha)
    if not isinstance(root, Mapping) or set(root) != {pdb_chain}:
        _fail(GrammarFailureCode.ATLAS_METADATA_CHAIN_KEY_INVALID,
              f"keys={sorted(root) if isinstance(root, Mapping) else type(root).__name__}",
              raw_len=len(data), raw_sha256=raw_sha)
    entry = root[pdb_chain]
    if not isinstance(entry, Mapping):
        _fail(GrammarFailureCode.ATLAS_METADATA_KEYS_INVALID, type(entry).__name__, raw_len=len(data), raw_sha256=raw_sha)
    unknown = sorted(set(entry) - set(contract.type_specs))
    missing = sorted(contract.required_fields - set(entry))
    if missing:
        _fail(GrammarFailureCode.ATLAS_METADATA_KEYS_INVALID, f"missing={missing}", raw_len=len(data), raw_sha256=raw_sha)
    type_mismatch_fields = []
    for key, value in entry.items():
        if key in contract.type_specs and key != "PDB" and not _check_json_type(value, contract.type_specs[key]):
            # amended 2026-09-04 phase-(ii) round 2: non-identity descriptor fields whose
            # real-world JSON type differs from the retained OpenAPI type are RECORDED,
            # not fatal (observed: alpha% numeric vs documented string). Identity stays strict.
            type_mismatch_fields.append(key)
    identity = entry["PDB"]
    if identity not in (pdb_chain, pdb_chain.replace("_", "")):
        # AMENDED 2026-09-04 (phase-(ii) round 3): real metadata carries the
        # chain-qualified form ('1a62_A'); both forms are accepted.
        _fail(GrammarFailureCode.ATLAS_METADATA_IDENTITY_MISMATCH,
              f"PDB={identity!r}:pdb_chain={pdb_chain!r}", raw_len=len(data), raw_sha256=raw_sha)
    payload = {
        "pdb_chain": pdb_chain,
        "pdb_identity": identity,
        "present_fields": sorted(entry),
        "unknown_fields_recorded": unknown,
        "type_mismatch_fields": sorted(set(type_mismatch_fields)),
    }
    return GrammarOutcome(role=SourceRole.ATLAS_METADATA, ok=True, raw_sha256=raw_sha, raw_len=len(data),
                          payload=payload, failure=None)


# --------------------------------------------------------------------------- #
# ATLAS parsable census member (E420_ATLAS_CENSUS_UNIVERSE_V1).               #
# --------------------------------------------------------------------------- #


@_grammar_boundary(SourceRole.ATLAS_CENSUS_UNIVERSE)
def validate_atlas_census_zip(
    raw: Any,
    *,
    expected_zip_sha256: str | None = ATLAS_CENSUS_ZIP_SHA256,
    expected_member_sha256: str | None = ATLAS_CENSUS_MEMBER_SHA256,
    member_path: str = ATLAS_CENSUS_MEMBER_PATH,
    expected_zip_len: int | None = None,
) -> GrammarOutcome:
    """Retained-census member grammar + pin binding (registration §3.4)."""
    data = _require_bytes(raw, expected_zip_sha256, expected_zip_len)
    raw_sha = _sha256(data)
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            names = set(archive.namelist())
            if member_path not in names:
                _fail(GrammarFailureCode.ATLAS_CENSUS_MEMBER_PATH_INVALID,
                      f"member={member_path!r}:archive={sorted(names)}", raw_len=len(data), raw_sha256=raw_sha)
            try:
                member = archive.read(member_path)
            except (zipfile.BadZipFile, RuntimeError, OSError, EOFError) as exc:
                _fail(GrammarFailureCode.ATLAS_CENSUS_ZIP_INVALID, f"member read: {type(exc).__name__}",
                      raw_len=len(data), raw_sha256=raw_sha)
    except zipfile.BadZipFile:
        _fail(GrammarFailureCode.ATLAS_CENSUS_ZIP_INVALID, "BadZipFile", raw_len=len(data), raw_sha256=raw_sha)
    member_sha = _sha256(member)
    if expected_member_sha256 is not None and member_sha != expected_member_sha256:
        _fail(GrammarFailureCode.RAW_BYTES_HASH_MISMATCH, "census member binding",
              raw_len=len(member), raw_sha256=member_sha, expected_sha256=str(expected_member_sha256))
    try:
        text = member.decode("ascii")
    except UnicodeDecodeError:
        _fail(GrammarFailureCode.ATLAS_CENSUS_ROW_GRAMMAR_INVALID, "member not ascii", raw_len=len(member), raw_sha256=member_sha)
    rows = text.split("\n")
    if rows and rows[-1] == "":
        rows.pop()
    seen: set[str] = set()
    entries: set[str] = set()
    for ordinal, row in enumerate(rows, start=1):
        if not re.fullmatch(r"[0-9A-Za-z]{4}_[A-Za-z0-9]+", row):
            _fail(GrammarFailureCode.ATLAS_CENSUS_ROW_GRAMMAR_INVALID, f"ordinal {ordinal}: {row!r}",
                  raw_len=len(member), raw_sha256=member_sha)
        if row in seen:
            _fail(GrammarFailureCode.ATLAS_CENSUS_ROW_GRAMMAR_INVALID, f"duplicate row: {row!r}",
                  raw_len=len(member), raw_sha256=member_sha)
        seen.add(row)
        entries.add(row.split("_", 1)[0])
    payload = {
        "member_path": member_path,
        "member_sha256": member_sha,
        "zip_sha256": raw_sha,
        "row_count": len(rows),
        "unique_row_count": len(seen),
        "unique_entry_count": len(entries),
    }
    return GrammarOutcome(role=SourceRole.ATLAS_CENSUS_UNIVERSE, ok=True, raw_sha256=raw_sha, raw_len=len(data),
                          payload=payload, failure=None)


# --------------------------------------------------------------------------- #
# Reference entry-files mmCIF reader (E420_REFERENCE_COORDINATE_V1).          #
# --------------------------------------------------------------------------- #

AUDIT_HISTORY_CATEGORY = "_pdbx_audit_revision_history"


def _tokenize_cif(text: str) -> list[str]:
    tokens: list[str] = []
    index = 0
    length = len(text)
    while index < length:
        char = text[index]
        if char in " \t\r\n":
            index += 1
            continue
        if char == "#":
            while index < length and text[index] != "\n":
                index += 1
            continue
        if char == ";" and (index == 0 or text[index - 1] == "\n"):
            end = text.find("\n;", index)
            if end == -1:
                _fail(GrammarFailureCode.REFERENCE_CIF_TOKEN_INVALID, "unterminated text field")
            tokens.append(text[index + 1:end])
            index = end + 2
            continue
        if char in "'\"":
            quote = char
            end = index + 1
            while True:
                end = text.find(quote, end)
                if end == -1:
                    _fail(GrammarFailureCode.REFERENCE_CIF_TOKEN_INVALID, "unterminated quoted value")
                after = text[end + 1:end + 2]
                if after in ("", " ", "\t", "\r", "\n"):
                    break
                end += 1
            tokens.append(text[index + 1:end])
            index = end + 1
            continue
        start = index
        while index < length and text[index] not in " \t\r\n":
            index += 1
        tokens.append(text[start:index])
    return tokens


@_grammar_boundary(SourceRole.REFERENCE_COORDINATE)
def validate_reference_cif(
    raw: Any,
    *,
    expected_raw_sha256: str | None = None,
    expected_raw_len: int | None = None,
) -> GrammarOutcome:
    """Entry-files mmCIF reader: data block, entry id, audit-history ordinal-1 date.

    Grammar-only: coordinate values are never read or exported.
    """
    data = _require_bytes(raw, expected_raw_sha256, expected_raw_len)
    raw_sha = _sha256(data)
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        _fail(GrammarFailureCode.REFERENCE_CIF_DECODE_FAILED, "not utf-8", raw_len=len(data), raw_sha256=raw_sha)
    tokens = _tokenize_cif(text)
    data_blocks = [token for token in tokens if token.startswith("data_")]
    if not data_blocks:
        _fail(GrammarFailureCode.REFERENCE_CIF_DATA_BLOCK_INVALID, "no data_ block", raw_len=len(data), raw_sha256=raw_sha)
    if len(data_blocks) > 1:
        _fail(GrammarFailureCode.REFERENCE_CIF_DATA_BLOCK_INVALID, f"blocks={data_blocks}",
              raw_len=len(data), raw_sha256=raw_sha)
    # Category-targeted extraction (amended 2026-09-04, phase-(ii) round 2):
    # the reader extracts exactly the registered fields (data block, _entry.id,
    # audit-history ordinal-1 revision_date). Full-file mmCIF walking is out of
    # scope: legacy files carry unquoted "_"-prefixed VALUE tokens inside loops
    # (e.g. _atom_site. in _pdbx_audit_revision_item.item), which make a strict
    # whole-file walk ambiguous. Coordinates are never read.
    data_block_tokens = [token for token in tokens if token.startswith("data_")]
    if not data_block_tokens:
        _fail(GrammarFailureCode.REFERENCE_CIF_DATA_BLOCK_INVALID, "no data_ block", raw_len=len(data), raw_sha256=raw_sha)
    entry_id = None
    for index, token in enumerate(tokens):
        if token == "_entry.id" and index + 1 < len(tokens):
            entry_id = tokens[index + 1]
            break
    if entry_id is None:
        _fail(GrammarFailureCode.REFERENCE_CIF_TOKEN_INVALID, "entry id absent", raw_len=len(data), raw_sha256=raw_sha)
    audit_prefix = AUDIT_HISTORY_CATEGORY + "."
    scalar_dates = {}
    for index, token in enumerate(tokens):
        if token in (f"{AUDIT_HISTORY_CATEGORY}.ordinal_id", f"{AUDIT_HISTORY_CATEGORY}.revision_date") \
                and index + 1 < len(tokens):
            scalar_dates[token] = tokens[index + 1]
    history_names: list[str] = []
    history_values: list[str] = []
    for index, token in enumerate(tokens):
        if token.startswith(audit_prefix):
            history_names.append(token)
            continue
        if history_names and not token.startswith("_") and token not in ("loop_", "data_", "stop_"):
            history_values.append(token)
            continue
        if history_values:
            break
    if not history_names:
        if scalar_dates.get(f"{AUDIT_HISTORY_CATEGORY}.ordinal_id") is not None and \
                scalar_dates.get(f"{AUDIT_HISTORY_CATEGORY}.revision_date") is not None:
            history_names = [f"{AUDIT_HISTORY_CATEGORY}.ordinal_id", f"{AUDIT_HISTORY_CATEGORY}.revision_date"]
            history_values = [scalar_dates[f"{AUDIT_HISTORY_CATEGORY}.ordinal_id"],
                              scalar_dates[f"{AUDIT_HISTORY_CATEGORY}.revision_date"]]
        else:
            _fail(GrammarFailureCode.REFERENCE_HISTORY_MISSING, "category absent", raw_len=len(data), raw_sha256=raw_sha)
    if len(history_values) % len(history_names) != 0 or not history_values:
        _fail(GrammarFailureCode.REFERENCE_CIF_TOKEN_INVALID,
              f"audit-history loop arity {len(history_values)}%{len(history_names)}",
              raw_len=len(data), raw_sha256=raw_sha)
    history_rows: list[dict[str, str]] = []
    columns = [name.split(".", 1)[1] for name in history_names]
    for start in range(0, len(history_values), len(history_names)):
        history_rows.append(dict(zip(columns, history_values[start:start + len(history_names)])))
    ordinal_col = next((name.split(".", 1)[1] for name in history_names
                        if name.split(".", 1)[1] in ("ordinal", "ordinal_id")), None)
    if ordinal_col is None or f"{AUDIT_HISTORY_CATEGORY}.revision_date" not in history_names:
        _fail(GrammarFailureCode.REFERENCE_HISTORY_MISSING, "ordinal/revision_date column absent",
              raw_len=len(data), raw_sha256=raw_sha)
    ordinal_one = [row for row in history_rows if re.fullmatch(r"0*1", row.get(ordinal_col, ""))]
    if not ordinal_one:
        _fail(GrammarFailureCode.REFERENCE_HISTORY_MISSING, "no ordinal_id=1 row", raw_len=len(data), raw_sha256=raw_sha)
    if len(ordinal_one) > 1:
        _fail(GrammarFailureCode.REFERENCE_HISTORY_AMBIGUOUS, f"ordinal-1 rows={len(ordinal_one)}",
              raw_len=len(data), raw_sha256=raw_sha)
    revision_date = ordinal_one[0].get("revision_date", "")
    if not REVISION_DATE_RE.fullmatch(revision_date):
        _fail(GrammarFailureCode.REFERENCE_REVISION_DATE_INVALID, f"revision_date={revision_date!r}",
              raw_len=len(data), raw_sha256=raw_sha)
    try:
        _dt.datetime.strptime(revision_date, "%Y-%m-%d")
    except ValueError:
        _fail(GrammarFailureCode.REFERENCE_REVISION_DATE_INVALID, f"revision_date={revision_date!r}",
              raw_len=len(data), raw_sha256=raw_sha)
    observed_categories = sorted({name.split(".", 1)[0] for name in history_names})
    payload = {
        "data_block": data_blocks[0][len("data_"):] or None,
        "entry_id": entry_id,
        "initial_release_date": revision_date,
        "history_row_count": len(history_rows),
        "observed_categories": observed_categories,
    }
    return GrammarOutcome(role=SourceRole.REFERENCE_COORDINATE, ok=True, raw_sha256=raw_sha, raw_len=len(data),
                          payload=payload, failure=None)


# --------------------------------------------------------------------------- #
# Release/version lock (registration §7 mutation check 4).                    #
# --------------------------------------------------------------------------- #


def _release_lock_failure(pinned: Mapping[str, Any], observed: Mapping[str, Any]) -> str | None:
    pinned_keys = set(pinned)
    observed_keys = set(observed)
    if pinned_keys != observed_keys:
        return f"keys pinned={sorted(pinned_keys)} observed={sorted(observed_keys)}"
    for key in sorted(pinned_keys & observed_keys):
        expected_value = pinned[key]
        observed_value = observed[key]
        if isinstance(expected_value, list):
            expected_value = tuple(expected_value)
        if isinstance(observed_value, list):
            observed_value = tuple(observed_value)
        if isinstance(expected_value, bool) != isinstance(observed_value, bool) or expected_value != observed_value:
            return f"{key}: pinned={pinned[key]!r} observed={observed[key]!r}"
        if isinstance(expected_value, float) or isinstance(observed_value, float):
            if not (_is_number(expected_value) and _is_number(observed_value)):
                return f"{key}: type pinned={type(pinned[key]).__name__} observed={type(observed[key]).__name__}"
        elif type(expected_value) is not type(observed_value) and not (
            isinstance(expected_value, tuple) and isinstance(observed_value, tuple)
        ):
            return f"{key}: type pinned={type(pinned[key]).__name__} observed={type(observed[key]).__name__}"
    return None


@_grammar_boundary(SourceRole.AFDB_METADATA)
def validate_release_lock(pinned: Any, observed: Any, *, role: SourceRole) -> GrammarOutcome:
    """Exact release/version scalar lock; any scalar change is a typed rejection."""
    @_grammar_boundary(role)
    def _run() -> GrammarOutcome:
        if not isinstance(pinned, Mapping) or not pinned:
            _fail(GrammarFailureCode.BOUNDARY_INVALID, "pinned not a nonempty mapping")
        if not isinstance(observed, Mapping) or not observed:
            _fail(GrammarFailureCode.BOUNDARY_INVALID, "observed not a nonempty mapping")
        detail = _release_lock_failure(pinned, observed)
        if detail is not None:
            _fail(GrammarFailureCode.RELEASE_LOCK_MISMATCH, detail)
        return GrammarOutcome(role=role, ok=True, raw_sha256=_sha256(b""), raw_len=0,
                              payload={"fields_locked": sorted(pinned)}, failure=None)
    return _run()


# --------------------------------------------------------------------------- #
# Fold/group preimage (registration §6/§7 mutation check 5).                  #
# --------------------------------------------------------------------------- #


def fold_hex(group_id: str) -> str:
    """Registered fold order: SHA256("g1-successor-fold-v1:" || group_id)."""
    if not isinstance(group_id, str) or not group_id:
        raise E420GrammarPinError("group_id invalid")
    return _sha256(f"{FOLD_SALT}{group_id}".encode("utf-8"))


def build_fold_preimage_row(group_id: str, accession: str) -> dict[str, str]:
    return {"group_id": group_id, "accession": accession, "fold_hex": fold_hex(group_id)}


def aggregate_fold_preimage_hash(rows: Sequence[Mapping[str, str]]) -> str:
    return _canonical_digest([dict(row) for row in rows])


def validate_fold_preimage(rows: Any, expected_aggregate_sha256: Any, *, role: SourceRole) -> GrammarOutcome:
    """Row-consistency + aggregate-hash validation of the fold/group preimage."""
    @_grammar_boundary(role)
    def _run() -> GrammarOutcome:
        if not isinstance(rows, (list, tuple)) or not rows:
            _fail(GrammarFailureCode.FOLD_PREIMAGE_ROW_INVALID, "rows empty or not a sequence")
        if not isinstance(expected_aggregate_sha256, str) or not HEX64_RE.fullmatch(expected_aggregate_sha256):
            _fail(GrammarFailureCode.FOLD_PREIMAGE_AGGREGATE_MISMATCH, "expected digest format")
        canonical_rows: list[dict[str, str]] = []
        for row in rows:
            if not isinstance(row, Mapping) or set(row) != {"group_id", "accession", "fold_hex"}:
                _fail(GrammarFailureCode.FOLD_PREIMAGE_ROW_INVALID,
                      f"keys={sorted(row) if isinstance(row, Mapping) else type(row).__name__}")
            if not isinstance(row["group_id"], str) or not row["group_id"] \
                    or not isinstance(row["accession"], str) or not row["accession"]:
                _fail(GrammarFailureCode.FOLD_PREIMAGE_ROW_INVALID, "field types")
            if not isinstance(row["fold_hex"], str) or not HEX64_RE.fullmatch(row["fold_hex"]):
                _fail(GrammarFailureCode.FOLD_PREIMAGE_ROW_INVALID, "fold_hex format")
            if row["fold_hex"] != fold_hex(row["group_id"]):
                _fail(GrammarFailureCode.FOLD_PREIMAGE_ROW_INVALID,
                      f"group={row['group_id']!r}: fold_hex mismatch", raw_sha256=row["fold_hex"],
                      expected_sha256=fold_hex(row["group_id"]))
            canonical_rows.append(dict(row))
        aggregate = aggregate_fold_preimage_hash(canonical_rows)
        if aggregate != expected_aggregate_sha256:
            _fail(GrammarFailureCode.FOLD_PREIMAGE_AGGREGATE_MISMATCH, "aggregate digest mismatch",
                  raw_sha256=aggregate, expected_sha256=expected_aggregate_sha256)
        return GrammarOutcome(role=role, ok=True, raw_sha256=aggregate, raw_len=len(canonical_rows),
                              payload={"row_count": len(canonical_rows), "aggregate_sha256": aggregate},
                              failure=None)
    return _run()


__all__ = [
    "SourceRole", "SegmentMode", "GrammarFailureCode", "GrammarFailure", "GrammarOutcome",
    "E420GrammarPinError", "ProvenanceSchemaBinding", "AfdbGrammarContract", "AtlasMetadataContract",
    "load_provenance_schema_binding", "load_afdb_contract", "load_atlas_metadata_contract",
    "build_row_envelope", "validate_row_envelope",
    "validate_pdbe_primary", "validate_pdbe_segment", "check_pdbe_segment_thresholds",
    "validate_afdb_metadata", "validate_afdb_child_binding", "validate_af_temporal",
    "validate_atlas_analysis_zip", "validate_atlas_metadata", "validate_atlas_census_zip",
    "validate_reference_cif", "validate_release_lock", "validate_gap_share",
    "fold_hex", "build_fold_preimage_row", "aggregate_fold_preimage_hash", "validate_fold_preimage",
    "EMITTED_BY", "PLDDT_FAMILY_FIELDS", "RMSF_HEADER", "CORRESP_HEADER",
    "ATLAS_CENSUS_ZIP_SHA256", "ATLAS_CENSUS_MEMBER_SHA256", "ATLAS_CENSUS_MEMBER_PATH",
    "PDBe_V2_OPENAPI_SHA256", "AFDB_OPENAPI_SHA256", "ATLAS_OPENAPI_SHA256",
    "AFDB_OPTIONAL_TABLE_SHA256", "PROVENANCE_SCHEMA_SHA256", "SIFTS_GAP_DERIVATION_SOURCE",
]
