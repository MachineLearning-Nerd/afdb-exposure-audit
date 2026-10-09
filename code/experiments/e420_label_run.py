#!/usr/bin/env python3
"""e420 gate (g): label-construction RUNNER (bounded live fetches; POOL_FREEZE sealed).

Registration: docs/e420_registration.md sections 3.5, 4.4, 5, 7, 10(g), 11,
12. Schedule/chronology: docs/g1_successor_analysis_plan.md sections 4-6.
Frozen thresholds: docs/g1_successor_label_rules.json. Lanes:
docs/g1_successor_source_allowlist.json.

The runner PERFORMS the registered gate-(g) fetches for the admitted pool
accessions and builds labeled rows.

Chronology (registration 10(g); analysis plan section 6), per scope:

  (0) AFDB metadata is NEVER refetched.  The census-retained bytes
      (results/e420/source_availability_manifest.json raw_refs +
      data/e420/census_raw/<sha256>) are re-validated through the phase-(i)
      validator (gv.validate_afdb_metadata); the SELECTED object (earliest
      modelCreatedDate, then lexicographic modelEntityId; registration
      15.2) supplies the frozen cifUrl and the per-entry pins.
  (a) ONE AF coordinate-child GET per admitted accession, exactly the
      selected metadata row's cifUrl field (pdbUrl/bcifUrl are forbidden
      substitutes), bound by gv.validate_afdb_child_binding.
  (b) ONE reference-coordinate GET per admitted census entry
      (https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif), validated by
      gv.validate_reference_cif, then bound temporally per registration 5.1
      by gv.validate_af_temporal: modelCreatedDate STRICTLY earlier than
      the in-file ordinal-1 revision_date; mismatch = typed temporal
      failure row (no labels for that accession/entry; ATLAS is not opened
      for it).
  (c) ONE ATLAS analysis GET per admitted census {entry}_{chain}
      (https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_chain});
      the frozen member grammar must validate (gv.validate_atlas_analysis_zip)
      or the scope is a typed failure.

Labels (registration 5; thresholds frozen in the label rules):

  E channel -- chain-scoped CA rows (reference: the census row's chain via
  auth_asym_id; AF child: the first polymer chain in file order via
  label_asym_id, falling back to auth_asym_id) through
  labels.kabsch_fit + labels.lddt_per_residue + labels.e_residue_state,
  matched on UniProt positions through the PDBe primary mapping bytes
  retained at census (reused, never refetched) with the registration-4.4
  guards: one-to-one indexType=='PDB' segments only, observed CA count >= 30
  on both sides and on the matched set, amino-acid agreement >= 0.97 over
  the canonical AFDB sequence span, rank-degenerate fits -> per-residue U.

  F channel -- ATLAS RMSF.tsv three-replicate mean (labels.f_residue_state;
  math.fsum/3) joined to UniProt positions through the corresp table + the
  primary mapping (SIFTS path, registration 5.2); gap rows and SIFTS-absent
  rows stay per-row U with full no-drop accounting
  (f_reason=missing_atlas_join / atlas_gap); the 5*unmapped > total
  profile rejection is a typed scope failure (registration 5.4).

  Combined -- labels.combine_residue_labels (registration 5.3): each
  residue exactly one of OC/A/B/O/U; any U endpoint yields U with a
  retained reason.

Every ledger row carries the registration-7 provenance envelope
(gv.build_row_envelope; per-member raw sha refs under the v27 member keys
pdbe/sifts/afdb/reference/atlas) validated with gv.validate_row_envelope.
The per-row raw_row_sha256 is the raw bytes the row identity comes from:
the census-retained primary mapping bytes for combined/E rows, the ATLAS
analysis archive bytes for F-only U rows.

Schedule/budget/receipts exactly like experiments/e420_census.py:
ThreadPoolExecutor <= 2 workers; a global-lock inter-start gate enforces
>= 250 ms +/- 100 ms jitter between request starts; on 429/5xx the runner
backs off 2 s doubling, max 3 attempts, then a typed failure row; raw
response bytes under --raw-dir (default data/e420/labels_raw/) named by
content sha256; every response attempt receipted (per-accession tagged for
resume attribution).

Resume: atomic progress (ledger rows + typed failures + per-accession
receipts) rewritten every 10 completed accessions; --resume skips
completed accessions without refetching.

Outputs: results/e420/label_ledger.json,
results/e420/four_cell_census.json (the labels.four_cell_census payload,
digest intact; written only when every admitted accession completed),
results/e420/label_run_receipt.json (request counts, typed failures,
flags: labels_written=true, embeddings_read=false, models_run=0,
e412_started=false).

STRICT write scope: nothing is written outside --raw-dir and the three
output paths.  No embeddings, no models, no e412 -- fetch + labels only.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import sys
import threading
import time
import zipfile
from enum import Enum
from typing import Any, Callable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census_run  # noqa: E402  (request/receipt/resume pattern; REUSE)
import e420_grammar_validate as gv  # noqa: E402  (phase-(i) validators; USE)
import e420_labels as labels  # noqa: E402  (frozen label constructors; USE)
import e420_pool as pool  # noqa: E402  (freeze verification + digests)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --------------------------------------------------------------------------
# Registered paths (registration 12; raw cache under data/ per fence 7).
# --------------------------------------------------------------------------
POOL_FREEZE_REL = os.path.join("results", "e420", "pool_freeze.json")
MANIFEST_REL = os.path.join("results", "e420", "source_availability_manifest.json")
CENSUS_RAW_REL = os.path.join("data", "e420", "census_raw")
LEDGER_REL = os.path.join("results", "e420", "label_ledger.json")
CENSUS_OUT_REL = os.path.join("results", "e420", "four_cell_census.json")
RECEIPT_REL = os.path.join("results", "e420", "label_run_receipt.json")
LABELS_RAW_REL = os.path.join("data", "e420", "labels_raw")

SCHEMA = "e420-label-ledger-v1"
RECEIPT_SCHEMA = "e420-label-run-receipt-v1"
PROGRESS_SCHEMA = "e420-label-run-progress-v1"
PROGRESS_NAME = "e420_label_run_progress.json"
EMITTED_BY = "e420-label-run-v1"
POOL_FREEZE_SCHEMA = "e420-pool-freeze-v1"
MANIFEST_SCHEMA = census_run.SCHEMA

PROGRESS_EVERY_ACCESSIONS = 10
CONCURRENCY_MAX = census_run.CONCURRENCY_MAX  # 2 (analysis plan section 5)

# Registration 4.4 identity/mapping guards (always apply, over retained bytes).
OBSERVED_CA_MIN = 30
AA_AGREEMENT_MIN = 0.97

# Registered fetch endpoints (allowlist source roles; nothing else is fetched).
REFERENCE_URL = "https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"
ATLAS_ANALYSIS_URL = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_chain}"

ROLE_CHILD = gv.SourceRole.AFDB_COORDINATE_CHILD.value
ROLE_REFERENCE = gv.SourceRole.REFERENCE_COORDINATE.value
ROLE_ATLAS = gv.SourceRole.ATLAS_ANALYSIS.value
ROLE_PRIMARY = gv.SourceRole.PDBE_MAPPING_PRIMARY.value
ROLE_SEGMENT = gv.SourceRole.PDBE_MAPPING_SEGMENT.value
ROLE_AFDB_META = gv.SourceRole.AFDB_METADATA.value

# Envelope member keys (provenance schema v27 normalized_row_contract).
MEMBER_PDBE = "pdbe"
MEMBER_SIFTS = "sifts"
MEMBER_AFDB = "afdb"
MEMBER_REFERENCE = "reference"
MEMBER_ATLAS = "atlas"

# Typed unknown carried from the allowlist (fence 5: never invented).
SIFTS_RELEASE_TYPED_UNKNOWN = "TYPED_UNKNOWN:PDBe_SIFTS_MAPPING_RELEASE_NOT_DOCUMENTED"

CHILD_UNAVAILABLE = "E420_AF_CHILD_UNAVAILABLE"
REFERENCE_UNAVAILABLE = "E420_REFERENCE_COORDINATE_ABSENT"
ATLAS_UNAVAILABLE = "E420_ATLAS_ANALYSIS_ABSENT"


class LabelRunFailureCode(str, Enum):
    """Runner-level typed failure codes (validator/label codes pass through
    verbatim)."""

    POOL_FREEZE_INVALID = "E420_POOL_FREEZE_INVALID"
    CENSUS_MANIFEST_INVALID = "E420_CENSUS_MANIFEST_INVALID"
    CENSUS_MANIFEST_SHA_MISMATCH = "E420_CENSUS_MANIFEST_SHA_MISMATCH"
    ADMITTED_ROW_MISSING = "E420_ADMITTED_CENSUS_ROW_MISSING"
    AFDB_BYTES_UNAVAILABLE = "E420_AFDB_CENSUS_BYTES_UNAVAILABLE"
    MAPPING_BYTES_UNAVAILABLE = "E420_MAPPING_CENSUS_BYTES_UNAVAILABLE"
    AF_TEMPORAL_FAILURE = "E420_AF_TEMPORAL_FAILURE"
    REFERENCE_CHAIN_NO_CA = "E420_REFERENCE_CHAIN_NO_CA"
    AF_CHILD_NO_CA = "E420_AF_CHILD_NO_CA"
    MAPPING_GUARD_FAILURE = "E420_MAPPING_GUARD_FAILURE"
    CIF_PARSE_FAILURE = "E420_CIF_PARSE_FAILURE"
    CIF_CHAIN_TAG_MISSING = "E420_CIF_CHAIN_TAG_MISSING"
    ATLAS_MEMBER_PARSE_FAILURE = "E420_ATLAS_MEMBER_PARSE_FAILURE"
    ATLAS_JOIN_REJECTED = "E420_ATLAS_JOIN_REJECTED"
    ENVELOPE_INVALID = "E420_ENVELOPE_INVALID"
    INTERNAL_ERROR = "E420_INTERNAL_ERROR"


class LabelRunError(RuntimeError):
    """Typed pre-network setup error (pool freeze, census manifest, args)."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


class CifParseError(Exception):
    """Typed chain-scoped CA parse failure (converted to a failure row)."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


class ScopeFailure(Exception):
    """Typed per-scope failure (guard / join rejection), retained, no repair."""

    def __init__(self, code: str, detail: str, kind: str = "guard") -> None:
        self.code = code
        self.detail = detail
        self.kind = kind
        super().__init__(f"{code}: {detail}")


def _utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _clip(text: Any, limit: int = 300) -> str:
    return "" if text is None else str(text)[:limit]


# --------------------------------------------------------------------------
# Chain-scoped CA parsing (frozen selection mirrors labels.parse_mmcif_ca:
# label_atom_id == CA, label_comp_id in the standard AA set, model 1; the
# chain id is retained instead of dropped).
# --------------------------------------------------------------------------
_CA_BASE_TAGS = ("label_atom_id", "label_comp_id", "auth_seq_id",
                 "cartn_x", "cartn_y", "cartn_z", "pdbx_pdb_model_num")


def parse_ca_rows(raw: Any, *, chain_field: str) -> list[dict[str, Any]]:
    """Chain-scoped standard-AA model-1 CA rows from entry-files mmCIF /
    ModelCIF: [{"auth_seq_id", "comp3", "chain", "x", "y", "z"}, ...] in
    file order.  No imputation, renumbering, or silent drops; malformed
    input is a typed CifParseError."""
    if not isinstance(raw, (bytes, bytearray)) or not bytes(raw):
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                            f"raw bytes type={type(raw).__name__}")
    try:
        text = bytes(raw).decode("utf-8")
    except UnicodeDecodeError:
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value, "not utf-8")
    try:
        tokens = gv._tokenize_cif(text)
        loops = labels._extract_cif_loops(tokens)
    except gv._GrammarError as exc:
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                            _clip(exc.detail)) from None
    atom_loops = [(tags, values) for tags, values in loops
                  if any(tag.lower().startswith("_atom_site.") for tag in tags)]
    if not atom_loops:
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                            "no _atom_site loop")
    if len(atom_loops) > 1:
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                            f"multiple _atom_site loops ({len(atom_loops)})")
    tags, values = atom_loops[0]
    lower_tags = [t.lower().split(".", 1)[1] if "." in t.lower() else t.lower()
                  for t in tags]
    required = set(_CA_BASE_TAGS) | {chain_field}
    missing = sorted(required - set(lower_tags))
    if missing:
        code = (LabelRunFailureCode.CIF_CHAIN_TAG_MISSING.value
                if missing == [chain_field]
                else LabelRunFailureCode.CIF_PARSE_FAILURE.value)
        raise CifParseError(code, f"missing {missing}")
    if not lower_tags or len(values) % len(lower_tags) != 0:
        raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                            f"loop arity {len(values)} % {len(lower_tags) or 0} != 0")
    col = {name: pos for pos, name in enumerate(lower_tags)}
    out: list[dict[str, Any]] = []
    total_rows = len(values) // len(lower_tags)
    for row_index in range(total_rows):
        row = values[row_index * len(lower_tags):(row_index + 1) * len(lower_tags)]
        if row[col["pdbx_pdb_model_num"]] != "1":
            continue
        if row[col["label_atom_id"]] != "CA":
            continue
        comp3 = row[col["label_comp_id"]].upper()
        if comp3 not in labels.STANDARD_AA3:
            continue
        try:
            auth_seq = int(row[col["auth_seq_id"]])
        except ValueError:
            raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                                f"row {row_index}: auth_seq_id={row[col['auth_seq_id']]!r}") from None
        xyz: list[float] = []
        for tag in ("cartn_x", "cartn_y", "cartn_z"):
            try:
                value = float(row[col[tag]])
            except ValueError:
                raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                                    f"row {row_index}: {tag}={row[col[tag]]!r}") from None
            if value != value or value in (float("inf"), float("-inf")):
                raise CifParseError(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                                    f"row {row_index}: {tag} nonfinite") from None
            xyz.append(value)
        out.append({"auth_seq_id": auth_seq, "comp3": comp3,
                    "chain": row[col[chain_field]],
                    "x": xyz[0], "y": xyz[1], "z": xyz[2]})
    return out


def scope_ca_rows(rows: list[dict[str, Any]], chain: str) -> list[dict[str, Any]]:
    """Rows of exactly one chain (order preserved)."""
    return [row for row in rows if row["chain"] == chain]


def first_polymer_chain(rows: list[dict[str, Any]]) -> tuple[str, list[dict[str, Any]]]:
    """First distinct chain id in file order among CA rows (AF child: the
    first polymer chain; CA rows exist only for polymers)."""
    for row in rows:
        return row["chain"], scope_ca_rows(rows, row["chain"])
    return "", []


# --------------------------------------------------------------------------
# Author-number -> UniProt position map (registration 4.4 declared one-to-one
# path) from the VALIDATED census-retained PDBe primary bytes.  Only
# identity fields are read (accession, entityId, bestChainId, residue
# ranges, indexType) per the frozen parse whitelist
# E420_PDBE_MAPPING_PRIMARY_V1; one-to-one indexType=='PDB' segments only.
# --------------------------------------------------------------------------
class AuthorMap:
    STATUS_MAPPED = "mapped"
    STATUS_UNMAPPED = "unmapped"
    STATUS_AMBIGUOUS = "ambiguous"

    def __init__(self, segments: list[tuple[int, int, int]], residue_count: int) -> None:
        self._segments = sorted(segments)  # (start, end, unp_start), one-to-one
        self.residue_count = int(residue_count)

    def map(self, author_num: int) -> tuple[str, int | None]:
        hits = {unp_start + (author_num - start)
                for start, end, unp_start in self._segments
                if start <= author_num <= end}
        if len(hits) == 1:
            return self.STATUS_MAPPED, hits.pop()
        if not hits:
            return self.STATUS_UNMAPPED, None
        return self.STATUS_AMBIGUOUS, None

    @property
    def segment_count(self) -> int:
        return len(self._segments)


def build_author_map(primary_root: Any, *, pdb_id: str, accession: str,
                     entity_id: int, chain: str) -> AuthorMap:
    """Identity-field scan of the validated primary envelope (raw JSON root)
    for (accession, entity_id, chain)."""
    segments: list[tuple[int, int, int]] = []
    residue_count = 0
    entry = primary_root.get(pdb_id) if isinstance(primary_root, dict) else None
    data_rows = entry.get("data") if isinstance(entry, dict) else None
    if isinstance(data_rows, list):
        for data_row in data_rows:
            if not isinstance(data_row, dict) or data_row.get("accession") != accession:
                continue
            additional = data_row.get("additionalData")
            if not isinstance(additional, dict) or additional.get("entityId") != entity_id:
                continue
            best_chain = additional.get("bestChainId")
            if best_chain is not None and best_chain != chain:
                continue
            residues = data_row.get("residues")
            if not isinstance(residues, list):
                continue
            for segment in residues:
                if not isinstance(segment, dict):
                    continue
                if segment.get("indexType") != "PDB":
                    continue
                start, end = segment.get("startIndex"), segment.get("endIndex")
                unp_start, unp_end = (segment.get("unpStartIndex"),
                                      segment.get("unpEndIndex"))
                if any(isinstance(v, bool) or not isinstance(v, int)
                       for v in (start, end, unp_start, unp_end)):
                    continue
                if start <= end and unp_start <= unp_end \
                        and (end - start) == (unp_end - unp_start):
                    segments.append((start, end, unp_start))
                    residue_count += end - start + 1
    return AuthorMap(segments, residue_count)


# --------------------------------------------------------------------------
# Pool inputs (pool freeze seal + census manifest binding).
# --------------------------------------------------------------------------
def load_pool_inputs(
    project_root: str,
    *,
    pool_freeze_path: str | None = None,
    manifest_path: str | None = None,
) -> dict[str, Any]:
    """Verify the sealed pool freeze (digest + fold preimage aggregate) and
    bind the census manifest by its sealed sha256.  Returns the admitted
    accessions, the frozen group map, and the per-accession census scopes."""
    pool_freeze_path = pool_freeze_path or os.path.join(project_root, POOL_FREEZE_REL)
    manifest_path = manifest_path or os.path.join(project_root, MANIFEST_REL)
    try:
        with open(pool_freeze_path, "rb") as handle:
            freeze_bytes = handle.read()
    except OSError as exc:
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            f"unreadable: {exc}") from None
    try:
        artifact = json.loads(freeze_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            type(exc).__name__) from None
    if not isinstance(artifact, dict) or artifact.get("schema") != POOL_FREEZE_SCHEMA:
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            "schema invalid")
    if not pool.verify_pool_freeze(artifact):
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            "digest / fold-preimage aggregate mismatch")
    selection = artifact.get("selection")
    admitted = selection.get("admitted") if isinstance(selection, dict) else None
    if not isinstance(admitted, list) or not admitted or \
            not all(isinstance(acc, str) and acc for acc in admitted):
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            "selection.admitted missing/empty")
    group_map = artifact.get("group_map")
    if not isinstance(group_map, dict) or not group_map:
        raise LabelRunError(LabelRunFailureCode.POOL_FREEZE_INVALID.value,
                            "group_map missing/empty")
    try:
        with open(manifest_path, "rb") as handle:
            manifest_bytes = handle.read()
    except OSError as exc:
        raise LabelRunError(LabelRunFailureCode.CENSUS_MANIFEST_INVALID.value,
                            f"unreadable: {exc}") from None
    manifest_sha256 = _sha256_hex(manifest_bytes)
    if manifest_sha256 != artifact.get("manifest_sha256"):
        raise LabelRunError(
            LabelRunFailureCode.CENSUS_MANIFEST_SHA_MISMATCH.value,
            f"manifest={manifest_sha256}:freeze={artifact.get('manifest_sha256')}")
    try:
        manifest = json.loads(manifest_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LabelRunError(LabelRunFailureCode.CENSUS_MANIFEST_INVALID.value,
                            type(exc).__name__) from None
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA \
            or not isinstance(manifest.get("rows"), list):
        raise LabelRunError(LabelRunFailureCode.CENSUS_MANIFEST_INVALID.value,
                            "schema/rows invalid")
    admitted_set = set(admitted)
    scopes: dict[str, list[dict[str, Any]]] = {acc: [] for acc in admitted}
    for position, row in enumerate(manifest["rows"]):
        if not isinstance(row, dict) or row.get("inclusion_status") != "included":
            continue
        entry, chain = row.get("entry"), row.get("chain")
        entity_id = row.get("entity_id")
        accessions = row.get("accessions")
        if not isinstance(entry, str) or not isinstance(chain, str) \
                or isinstance(entity_id, bool) or not isinstance(entity_id, int) \
                or not isinstance(accessions, list):
            raise LabelRunError(LabelRunFailureCode.CENSUS_MANIFEST_INVALID.value,
                                f"row {position} identity fields invalid")
        raw_refs = row.get("raw_refs")
        if not isinstance(raw_refs, dict):
            raw_refs = {}
        for acc in accessions:
            if isinstance(acc, str) and acc in admitted_set:
                scopes[acc].append({
                    "entry": entry,
                    "chain": chain,
                    "pdb_chain": f"{entry}_{chain}",
                    "entity_id": int(entity_id),
                    "accession": acc,
                    "row_key": row.get("row_key"),
                    "raw_refs": dict(raw_refs),
                })
    selection = artifact.get("selection") or {}
    cap_bound = bool(selection.get("cap_bound", False)) if isinstance(selection, dict) else False
    return {
        "admitted": sorted(admitted),
        "group_map": dict(group_map),
        "scopes": scopes,
        "manifest_sha256": manifest_sha256,
        "pool_freeze_sha256": _sha256_hex(freeze_bytes),
        "pool_freeze_digest": artifact.get("pool_freeze_digest"),
        "pool_cap_bound": cap_bound,
    }


# --------------------------------------------------------------------------
# Tagged scheduler (per-accession receipt attribution for resume; schedule
# discipline is the census module's own, unchanged).
# --------------------------------------------------------------------------
class TaggedRequestScheduler(census_run.RequestScheduler):
    """RequestScheduler + thread-local accession tags on receipts."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._tag_local = threading.local()

    def set_tag(self, accession: str | None) -> None:
        self._tag_local.accession = accession

    def get(self, role: gv.SourceRole, url: str) -> census_run.RequestResult:
        result = super().get(role, url)
        tag = getattr(self._tag_local, "accession", None)
        for receipt in result.receipts:
            receipt["accession"] = tag
        return result


def _sorted_receipts(receipts: list[dict]) -> list[dict]:
    def key(r: dict):
        sha = r.get("sha256") or ""
        status = r["status"] if isinstance(r.get("status"), int) else -1
        return (r.get("accession") or "", r["role"], r["url"], int(r["attempt"]),
                sha, status, int(r["bytes_len"]))
    return sorted(receipts, key=key)


# --------------------------------------------------------------------------
# Per-fetch steps (validator-first; census-retained bytes reused, never
# refetched).
# --------------------------------------------------------------------------
def _transport_failure(res: census_run.RequestResult) -> dict[str, Any]:
    return {"ok": False, "kind": "transport", "status": res.status,
            "raw": None, "sha256": None, "url": res.url, "payload": None,
            "code": res.failure_code, "detail": res.failure_detail,
            "raw_ref": res.raw_ref}


def _http_failure(res: census_run.RequestResult, *, absent_code: str) -> dict[str, Any]:
    if res.status == 404:
        code, kind = absent_code, "absent"
    else:
        code, kind = census_run.CensusFailureCode.HTTP_UNEXPECTED_STATUS.value, "unavailable"
    return {"ok": False, "kind": kind, "status": res.status, "raw": None,
            "sha256": None, "url": res.url, "payload": None, "code": code,
            "detail": _clip(f"status={res.status}"), "raw_ref": res.raw_ref}


def load_afdb_selected(census_raw_dir: str, raw_ref: dict | None,
                       accession: str, *, project_root: str) -> dict[str, Any]:
    """Re-validate the census-retained AFDB metadata bytes; the selected
    object (earliest modelCreatedDate, then modelEntityId; amendment 15.2)
    is the canonical identity/temporal/coordinate source."""
    sha = (raw_ref or {}).get("sha256")
    if not sha:
        return {"ok": False, "kind": "guard",
                "code": LabelRunFailureCode.AFDB_BYTES_UNAVAILABLE.value,
                "detail": f"no {ROLE_AFDB_META} raw ref for {accession}"}
    path = os.path.join(census_raw_dir, sha)
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError:
        return {"ok": False, "kind": "guard",
                "code": LabelRunFailureCode.AFDB_BYTES_UNAVAILABLE.value,
                "detail": f"census raw bytes missing: sha256={sha}"}
    outcome = gv.validate_afdb_metadata(raw, project_root=project_root,
                                        expected_accession=accession)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "code": outcome.failure.code.value,
                "detail": _clip(outcome.failure.detail)}
    return {"ok": True, "kind": "ok", "payload": outcome.payload,
            "sha256": outcome.raw_sha256, "raw_ref": dict(raw_ref or {})}


def fetch_child(scheduler: TaggedRequestScheduler, *, url: str,
                metadata_payload: dict) -> dict[str, Any]:
    """Chronology (a): exactly the selected metadata row's cifUrl, bound by
    the phase-(i) child-binding validator (no forbidden substitute)."""
    res = scheduler.get(gv.SourceRole.AFDB_COORDINATE_CHILD, url)
    if not res.ok:
        return _transport_failure(res)
    if res.status != 200:
        return _http_failure(res, absent_code=CHILD_UNAVAILABLE)
    outcome = gv.validate_afdb_child_binding(metadata_payload, url, res.content)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "status": res.status, "raw": None,
                "sha256": None, "url": url, "payload": None,
                "code": outcome.failure.code.value,
                "detail": _clip(outcome.failure.detail), "raw_ref": res.raw_ref}
    return {"ok": True, "kind": "ok", "status": res.status, "raw": res.content,
            "sha256": outcome.payload["child_sha256"], "url": url,
            "payload": outcome.payload, "raw_ref": res.raw_ref}


def fetch_reference(scheduler: TaggedRequestScheduler, entry: str) -> dict[str, Any]:
    """Chronology (b): the entry-files reference mmCIF, grammar-validated
    (coordinates are read later, only through the chain-scoped CA parser)."""
    url = REFERENCE_URL.format(entry=entry)
    res = scheduler.get(gv.SourceRole.REFERENCE_COORDINATE, url)
    if not res.ok:
        return _transport_failure(res)
    if res.status != 200:
        return _http_failure(res, absent_code=REFERENCE_UNAVAILABLE)
    outcome = gv.validate_reference_cif(res.content)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "status": res.status, "raw": None,
                "sha256": None, "url": url, "payload": None,
                "code": outcome.failure.code.value,
                "detail": _clip(outcome.failure.detail), "raw_ref": res.raw_ref}
    return {"ok": True, "kind": "ok", "status": res.status, "raw": res.content,
            "sha256": outcome.raw_sha256, "url": url, "payload": outcome.payload,
            "raw_ref": res.raw_ref}


def fetch_atlas(scheduler: TaggedRequestScheduler, pdb_chain: str) -> dict[str, Any]:
    """Chronology (c): the per-chain ATLAS analysis archive, member grammar
    validated (values are read later, only through build_f_evidence)."""
    url = ATLAS_ANALYSIS_URL.format(pdb_chain=pdb_chain)
    res = scheduler.get(gv.SourceRole.ATLAS_ANALYSIS, url)
    if not res.ok:
        return _transport_failure(res)
    if res.status != 200:
        return _http_failure(res, absent_code=ATLAS_UNAVAILABLE)
    outcome = gv.validate_atlas_analysis_zip(res.content, pdb_chain=pdb_chain)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "status": res.status, "raw": None,
                "sha256": None, "url": url, "payload": None,
                "code": outcome.failure.code.value,
                "detail": _clip(outcome.failure.detail), "raw_ref": res.raw_ref}
    return {"ok": True, "kind": "ok", "status": res.status, "raw": res.content,
            "sha256": outcome.raw_sha256, "url": url, "payload": outcome.payload,
            "raw_ref": res.raw_ref}


def load_primary(census_raw_dir: str, raw_ref: dict | None, *, entry: str,
                 entity_id: int, project_root: str) -> dict[str, Any]:
    """Re-validate the census-retained primary mapping bytes (reused, never
    refetched); keep the raw bytes (row-identity source) and the parsed
    identity root for the AuthorMap."""
    sha = (raw_ref or {}).get("sha256")
    if not sha:
        return {"ok": False, "kind": "guard",
                "code": LabelRunFailureCode.MAPPING_BYTES_UNAVAILABLE.value,
                "detail": f"no {ROLE_PRIMARY} raw ref for {entry}/{entity_id}"}
    path = os.path.join(census_raw_dir, sha)
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError:
        return {"ok": False, "kind": "guard",
                "code": LabelRunFailureCode.MAPPING_BYTES_UNAVAILABLE.value,
                "detail": f"census raw bytes missing: sha256={sha}"}
    outcome = gv.validate_pdbe_primary(raw, pdb_id=entry, entity_id=entity_id,
                                       project_root=project_root)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "code": outcome.failure.code.value,
                "detail": _clip(outcome.failure.detail)}
    try:
        root = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        return {"ok": False, "kind": "grammar",
                "code": LabelRunFailureCode.INTERNAL_ERROR.value,
                "detail": f"validated bytes reparse failed: {type(exc).__name__}"}
    return {"ok": True, "kind": "ok", "root": root, "raw": raw, "sha256": sha,
            "raw_ref": dict(raw_ref or {})}


class EntryCache:
    """Thread-safe fetch-once caches for entry-scoped / chain-scoped
    artifacts shared across accessions (the census fetch-once pattern)."""

    def __init__(self) -> None:
        self._reference = census_run.AfdbCache()
        self._atlas = census_run.AfdbCache()
        self._primary = census_run.AfdbCache()

    def reference(self, entry: str, fetch: Callable[[], dict]) -> dict:
        return self._reference.get(entry, fetch)

    def atlas(self, pdb_chain: str, fetch: Callable[[], dict]) -> dict:
        return self._atlas.get(pdb_chain, fetch)

    def primary(self, key: tuple[str, int], fetch: Callable[[], dict]) -> dict:
        return self._primary.get(key, fetch)


# --------------------------------------------------------------------------
# E / F evidence construction (registration 5.1 / 5.2; 4.4 guards).
# --------------------------------------------------------------------------
def build_e_evidence(
    *,
    child_raw: bytes,
    reference_raw: bytes,
    primary_root: Any,
    afdb_payload: dict,
    accession: str,
    entity_id: int,
    chain: str,
    entry: str,
) -> dict[str, Any]:
    """E-channel evidence for one census scope.  Raises ScopeFailure on the
    typed guard failures; deterministic; no repair."""
    # Reference side: the census row's chain (author numbering = join key).
    try:
        ref_all = parse_ca_rows(reference_raw, chain_field="auth_asym_id")
    except CifParseError as exc:
        raise ScopeFailure(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                           f"reference {entry}: {exc.detail}") from None
    ref_rows = scope_ca_rows(ref_all, chain)
    if not ref_rows:
        raise ScopeFailure(LabelRunFailureCode.REFERENCE_CHAIN_NO_CA.value,
                           f"entry={entry} chain={chain}: 0 model-1 standard-AA CA rows")
    # Predictor side: the AF child's first polymer chain (file order).
    try:
        try:
            child_all = parse_ca_rows(child_raw, chain_field="label_asym_id")
        except CifParseError as exc:
            if exc.code != LabelRunFailureCode.CIF_CHAIN_TAG_MISSING.value:
                raise
            child_all = parse_ca_rows(child_raw, chain_field="auth_asym_id")
    except CifParseError as exc:
        raise ScopeFailure(LabelRunFailureCode.CIF_PARSE_FAILURE.value,
                           f"af child {accession}: {exc.detail}") from None
    child_chain, child_rows = first_polymer_chain(child_all)
    if not child_rows:
        raise ScopeFailure(LabelRunFailureCode.AF_CHILD_NO_CA.value,
                           f"accession={accession}: 0 model-1 standard-AA CA rows")
    # Declared one-to-one mapping path (registration 4.4): author -> UniProt.
    author_map = build_author_map(primary_root, pdb_id=entry, accession=accession,
                                  entity_id=entity_id, chain=chain)
    sequence = afdb_payload["sequence"]
    sequence_start = int(afdb_payload["sequence_start"])
    sequence_end = int(afdb_payload["sequence_end"])
    af_by_unp: dict[int, dict[str, Any]] = {}
    af_collisions = 0
    for row in child_rows:
        unp = sequence_start + int(row["auth_seq_id"]) - 1
        if unp < sequence_start or unp > sequence_end:
            continue
        if unp in af_by_unp:
            af_collisions += 1
            continue
        af_by_unp[unp] = row
    ref_by_unp: dict[int, dict[str, Any]] = {}
    ref_conflicts: set[int] = set()
    status_counts = {AuthorMap.STATUS_MAPPED: 0, AuthorMap.STATUS_UNMAPPED: 0,
                     AuthorMap.STATUS_AMBIGUOUS: 0}
    aa_compared = 0
    aa_matches = 0
    for row in ref_rows:
        status, unp = author_map.map(int(row["auth_seq_id"]))
        status_counts[status] += 1
        if status != AuthorMap.STATUS_MAPPED:
            continue
        if unp in ref_by_unp:
            ref_conflicts.add(unp)  # two author numbers -> one UniProt position
            continue
        ref_by_unp[unp] = row
        if sequence_start <= unp <= sequence_end:
            code1 = gv.AA3_TO_1.get(row["comp3"])
            if code1 is not None:
                aa_compared += 1
                if code1 == sequence[unp - sequence_start]:
                    aa_matches += 1
    usable_ref = {u: r for u, r in ref_by_unp.items() if u not in ref_conflicts}
    matched_unps = sorted(set(usable_ref) & set(af_by_unp))
    # Registration 4.4 guards (always apply, over retained bytes only).
    if len(ref_rows) < OBSERVED_CA_MIN:
        raise ScopeFailure(LabelRunFailureCode.MAPPING_GUARD_FAILURE.value,
                           f"reference chain CA count {len(ref_rows)} < {OBSERVED_CA_MIN}")
    if len(child_rows) < OBSERVED_CA_MIN:
        raise ScopeFailure(LabelRunFailureCode.MAPPING_GUARD_FAILURE.value,
                           f"af child CA count {len(child_rows)} < {OBSERVED_CA_MIN}")
    if len(matched_unps) < OBSERVED_CA_MIN:
        raise ScopeFailure(LabelRunFailureCode.MAPPING_GUARD_FAILURE.value,
                           f"matched CA count {len(matched_unps)} < {OBSERVED_CA_MIN}")
    if aa_compared and (aa_matches / aa_compared) < AA_AGREEMENT_MIN:
        raise ScopeFailure(
            LabelRunFailureCode.MAPPING_GUARD_FAILURE.value,
            f"amino-acid agreement {aa_matches}/{aa_compared} < {AA_AGREEMENT_MIN}")
    mobile = labels.ca_coords_array(
        [(0, af_by_unp[u]["x"], af_by_unp[u]["y"], af_by_unp[u]["z"])
         for u in matched_unps])
    reference = labels.ca_coords_array(
        [(0, usable_ref[u]["x"], usable_ref[u]["y"], usable_ref[u]["z"])
         for u in matched_unps])
    fit = labels.kabsch_fit(mobile, reference)
    lddt = labels.lddt_per_residue(fit["superposed"], reference)
    e_map: dict[int, tuple[str, str]] = {}
    for index, unp in enumerate(matched_unps):
        if not fit["rank_ok"]:
            e_map[unp] = ("U", "u_rank_degenerate")
        elif lddt["statuses"][index] != "ok":
            e_map[unp] = ("U", lddt["statuses"][index])
        else:
            e_map[unp] = labels.e_residue_state(float(lddt["lddt"][index]),
                                                float(fit["distances"][index]))
    return {
        "e_map": e_map,
        "author_map": author_map,
        "child_chain": child_chain,
        "stats": {
            "reference_chain_ca": len(ref_rows),
            "af_child_chain_ca": len(child_rows),
            "matched_ca": len(matched_unps),
            "aa_compared": aa_compared,
            "aa_matches": aa_matches,
            "rank_ok": bool(fit["rank_ok"]),
            "author_segments": author_map.segment_count,
            "author_mapped": status_counts[AuthorMap.STATUS_MAPPED],
            "author_unmapped": status_counts[AuthorMap.STATUS_UNMAPPED],
            "author_ambiguous": status_counts[AuthorMap.STATUS_AMBIGUOUS],
            "ref_position_conflicts": len(ref_conflicts),
            "af_position_collisions": af_collisions,
            "residual_max_a": fit["residual_max"],
        },
    }


def _atlas_member_rows(atlas_raw: bytes, pdb_chain: str) -> tuple[list, list]:
    """Re-parse the two validated members (the validator already enforced
    the frozen grammar; this re-read yields the F-channel row values)."""
    try:
        with zipfile.ZipFile(zipfile.io.BytesIO(atlas_raw)) as archive:
            rmsf_bytes = archive.read(pdb_chain + gv.RMSF_MEMBER_SUFFIX)
            corr_bytes = archive.read(pdb_chain + gv.CORRESP_MEMBER_SUFFIX)
        rmsf_rows = gv._tsv_rows(rmsf_bytes, gv.RMSF_HEADER,
                                 gv.GrammarFailureCode.ATLAS_RMSF_HEADER_INVALID,
                                 gv.GrammarFailureCode.ATLAS_RMSF_ROW_INVALID, 4)
        corr_rows = gv._tsv_rows(corr_bytes, gv.CORRESP_HEADER,
                                 gv.GrammarFailureCode.ATLAS_CORRESP_HEADER_INVALID,
                                 gv.GrammarFailureCode.ATLAS_CORRESP_ROW_INVALID, 10)
    except (zipfile.BadZipFile, KeyError, gv._GrammarError) as exc:
        raise ScopeFailure(LabelRunFailureCode.ATLAS_MEMBER_PARSE_FAILURE.value,
                           _clip(f"{type(exc).__name__}: {exc}"), kind="grammar") from None
    return rmsf_rows, corr_rows


def build_f_evidence(*, atlas_raw: bytes, pdb_chain: str,
                     author_map: AuthorMap) -> dict[str, Any]:
    """F-channel evidence for one census scope (registration 5.2/5.4).
    Raises ScopeFailure on the typed profile rejection."""
    rmsf_rows, corr_rows = _atlas_member_rows(atlas_raw, pdb_chain)
    f_map: dict[int, tuple[str, str]] = {}
    u_rows: list[dict[str, Any]] = []
    unmapped = 0
    gap_rows = 0
    join_ambiguous = 0
    for ordinal, corr_row in enumerate(corr_rows, start=1):
        rmsf_row = rmsf_rows[ordinal - 1]
        flags = dict(zip(labels.CORRESP_FLAG_NAMES, corr_row[5:10]))
        unp_seq = corr_row[0]
        try:
            triple = [float(rmsf_row[1]), float(rmsf_row[2]), float(rmsf_row[3])]
            pnum = int(corr_row[2])
        except ValueError as exc:
            raise ScopeFailure(LabelRunFailureCode.ATLAS_MEMBER_PARSE_FAILURE.value,
                               f"ordinal {ordinal}: {type(exc).__name__}",
                               kind="grammar") from None
        if unp_seq == "-":
            gap_rows += 1
            unmapped += 1
            u_rows.append({"pdb_num": pnum, "f_state": "U",
                           "f_reason": "atlas_gap", "gap": True,
                           "ordinal": ordinal})
            continue
        status, unp = author_map.map(pnum)
        if status == AuthorMap.STATUS_MAPPED:
            f_map[unp] = labels.f_residue_state(triple, flags)
        else:
            unmapped += 1
            if status == AuthorMap.STATUS_AMBIGUOUS:
                join_ambiguous += 1
            u_rows.append({"pdb_num": pnum, "f_state": "U",
                           "f_reason": "missing_atlas_join", "gap": False,
                           "ordinal": ordinal})
    share = labels.gap_share_status(unmapped, len(corr_rows))
    if not share["accepted"]:
        raise ScopeFailure(LabelRunFailureCode.ATLAS_JOIN_REJECTED.value,
                           f"pdb_chain={pdb_chain}: {share['detail']}",
                           kind="join_rejected")
    return {
        "f_map": f_map,
        "u_rows": u_rows,
        "stats": {
            "atlas_rows": len(corr_rows),
            "f_positions": len(f_map),
            "gap_rows": gap_rows,
            "join_unmapped_rows": unmapped - gap_rows,
            "join_ambiguous_rows": join_ambiguous,
        },
    }


# --------------------------------------------------------------------------
# The runner.
# --------------------------------------------------------------------------
class LabelRunRunner:
    """Gate-(g) label construction over the admitted pool accessions."""

    def __init__(
        self,
        *,
        ledger_path: str,
        census_out_path: str,
        receipt_path: str,
        raw_dir: str,
        project_root: str = ROOT,
        inputs: dict | None = None,
        pool_freeze_path: str | None = None,
        manifest_path: str | None = None,
        census_raw_dir: str | None = None,
        transport: census_run.Transport | None = None,
        resume: bool = False,
        concurrency: int = CONCURRENCY_MAX,
        min_delay_s: float = census_run.INTER_REQUEST_DELAY_MS_MIN / 1000.0,
        jitter_s: float = census_run.INTER_REQUEST_JITTER_MS / 1000.0,
        backoff_initial_s: float = census_run.BACKOFF_INITIAL_S,
        backoff_factor: float = census_run.BACKOFF_FACTOR,
        max_attempts: int = census_run.MAX_ATTEMPTS,
        timeout_s: float = census_run.REQUEST_TIMEOUT_S,
        seed: int = 20260904,
        sleeper: Callable[[float], None] = time.sleep,
        flush_hook: Callable[[str, str], None] | None = None,
    ) -> None:
        if concurrency > CONCURRENCY_MAX:
            raise LabelRunError("E420_CONCURRENCY_ABOVE_REGISTERED_BOUND",
                                f"concurrency={concurrency} > {CONCURRENCY_MAX}")
        if max_attempts < 1:
            raise LabelRunError("E420_INVALID_ARGUMENT", "max_attempts < 1")
        self._project_root = project_root
        self._census_raw_dir = census_raw_dir or os.path.join(project_root, CENSUS_RAW_REL)
        self._ledger_path = ledger_path
        self._census_out_path = census_out_path
        self._receipt_path = receipt_path
        self._raw_dir = raw_dir
        self._resume = bool(resume)
        self._concurrency = int(concurrency)
        self._seed = int(seed)
        self._flush_hook = flush_hook
        self._inputs = inputs if inputs is not None else load_pool_inputs(
            project_root, pool_freeze_path=pool_freeze_path,
            manifest_path=manifest_path)
        os.makedirs(self._raw_dir, exist_ok=True)
        self._progress_path = os.path.join(self._raw_dir, PROGRESS_NAME)
        self._scheduler = TaggedRequestScheduler(
            transport if transport is not None else census_run.RequestsTransport(),
            self._raw_dir,
            min_delay_s=min_delay_s, jitter_s=jitter_s,
            backoff_initial_s=backoff_initial_s, backoff_factor=backoff_factor,
            max_attempts=max_attempts, timeout_s=timeout_s, seed=seed,
            sleeper=sleeper)
        self._caches = EntryCache()
        self._results: dict[str, dict] = {}
        self._results_lock = threading.Lock()
        self._flush_lock = threading.Lock()
        self._completed_since_flush = 0
        self._flush_count = 0
        self._started_utc = _utc_now()
        self._resumed_accessions = 0

    # -- inputs ---------------------------------------------------------------
    @property
    def inputs(self) -> dict:
        return self._inputs

    def _context(self, afdb_payload: dict, entry: str,
                 release_date: str | None) -> dict:
        release_ids = {
            "afdb": (f"{afdb_payload.get('model_entity_id')}:"
                     f"v{afdb_payload.get('latest_version')}:"
                     f"{afdb_payload.get('model_created_date')}:"
                     f"{afdb_payload.get('sequence_version_date')}"),
            "sifts": SIFTS_RELEASE_TYPED_UNKNOWN,
        }
        if release_date:
            release_ids["reference"] = f"pdb_entry:{entry}:initial_release:{release_date}"
        return {
            "emitted_by": EMITTED_BY,
            "release_ids": release_ids,
            "capture_anchor": ("e420-census-manifest-sha256:"
                               + self._inputs["manifest_sha256"]),
        }

    # -- row envelope (registration 7) ----------------------------------------
    def _envelope(self, *, members: dict[str, dict[str, str]], raw_bytes: bytes,
                  normalized: dict, gap_derivation: dict | None,
                  modelcif_sha: str | None, context: dict) -> dict:
        envelope = gv.build_row_envelope(
            role=gv.SourceRole.PDBE_MAPPING_PRIMARY,
            raw_bytes=raw_bytes,
            normalized_projection=normalized,
            member_hashes=members,
            context=context,
            gap_derivation=gap_derivation,
            modelcif_projection_sha256=modelcif_sha,
        )
        outcome = gv.validate_row_envelope(
            envelope, project_root=self._project_root,
            applicable_member_keys=sorted(members),
            gap_applicable=gap_derivation is not None,
            modelcif_applicable=modelcif_sha is not None)
        if not outcome.ok:
            raise ScopeFailure(LabelRunFailureCode.ENVELOPE_INVALID.value,
                               _clip(outcome.failure.detail), kind="internal")
        return envelope

    def _failure_row(self, scope: dict, code: str, detail: str,
                     kind: str) -> dict:
        return {"accession": scope["accession"], "entry": scope["entry"],
                "chain": scope["chain"], "pdb_chain": scope["pdb_chain"],
                "code": code, "detail": _clip(detail), "kind": kind}

    # -- scope pipeline ---------------------------------------------------------
    def _process_scope(self, scope: dict, *, afdb_payload: dict,
                       child: dict) -> tuple[list, list]:
        """E+F labels for one census scope.  Returns (ledger_rows,
        typed_failures); typed failures are returned, never raised past the
        accession boundary."""
        accession = scope["accession"]
        entry = scope["entry"]
        chain = scope["chain"]
        pdb_chain = scope["pdb_chain"]
        entity_id = scope["entity_id"]
        rows: list[dict] = []
        failures: list[dict] = []
        raw_refs = scope.get("raw_refs") or {}

        # Chronology (b): reference coordinates (shared per entry).
        reference = self._caches.reference(
            entry, lambda: fetch_reference(self._scheduler, entry))
        if not reference["ok"]:
            return [], [self._failure_row(scope, reference["code"],
                                          reference.get("detail") or "",
                                          reference.get("kind", "transport"))]
        temporal = gv.validate_af_temporal(afdb_payload["model_created_date"],
                                           reference["raw"])
        if not temporal.ok:
            return [], [self._failure_row(
                scope, LabelRunFailureCode.AF_TEMPORAL_FAILURE.value,
                f"{temporal.failure.code.value}: {temporal.failure.detail}",
                kind="temporal")]
        # Chronology (c): ATLAS analysis archive (shared per pdb_chain).
        atlas = self._caches.atlas(
            pdb_chain, lambda: fetch_atlas(self._scheduler, pdb_chain))
        if not atlas["ok"]:
            return [], [self._failure_row(scope, atlas["code"],
                                          atlas.get("detail") or "",
                                          atlas.get("kind", "transport"))]
        # Census-retained primary mapping (reused, never refetched).
        primary = self._caches.primary(
            (entry, entity_id),
            lambda: load_primary(self._census_raw_dir, raw_refs.get(ROLE_PRIMARY),
                                 entry=entry, entity_id=entity_id,
                                 project_root=self._project_root))
        if not primary["ok"]:
            return [], [self._failure_row(scope, primary["code"],
                                          primary.get("detail") or "",
                                          primary.get("kind", "guard"))]
        try:
            e_out = build_e_evidence(
                child_raw=child["raw"], reference_raw=reference["raw"],
                primary_root=primary["root"], afdb_payload=afdb_payload,
                accession=accession, entity_id=entity_id, chain=chain,
                entry=entry)
            f_out = build_f_evidence(atlas_raw=atlas["raw"], pdb_chain=pdb_chain,
                                     author_map=e_out["author_map"])
        except ScopeFailure as exc:
            return [], [self._failure_row(scope, exc.code, exc.detail, exc.kind)]

        # Row envelope members (registration 7; per-member raw sha refs).
        members: dict[str, dict[str, str]] = {
            MEMBER_PDBE: {f"{entry}/entity_{entity_id}_uniprot_mapping.json":
                          primary["sha256"]},
            MEMBER_ATLAS: {f"{pdb_chain}_analysis.zip": atlas["sha256"]},
            MEMBER_AFDB: {f"{accession}/prediction.json":
                          scope.get("afdb_meta_sha256", ""),
                          f"{accession}/modelcif_child.cif": child["sha256"]},
            MEMBER_REFERENCE: {f"{entry}.cif": reference["sha256"]},
        }
        if not members[MEMBER_AFDB][f"{accession}/prediction.json"]:
            members[MEMBER_AFDB].pop(f"{accession}/prediction.json")
        segment_ref = raw_refs.get(ROLE_SEGMENT)
        if isinstance(segment_ref, dict) and segment_ref.get("sha256"):
            members[MEMBER_SIFTS] = {f"{entry}_mappings_uniprot.json":
                                     segment_ref["sha256"]}
        context = self._context(afdb_payload, entry,
                                (reference.get("payload") or {})
                                .get("initial_release_date"))

        # F-only U rows (gap / SIFTS-join miss): uniprot_num null, cell U,
        # full no-drop accounting; raw row bytes = the ATLAS archive.
        for u_row in f_out["u_rows"]:
            gap_derivation = None
            if u_row["gap"]:
                gap_derivation = {"derivation_source": gv.SIFTS_GAP_DERIVATION_SOURCE,
                                  "uniprot_num": None}
                u_members = {MEMBER_ATLAS: members[MEMBER_ATLAS]}
            else:
                u_members = {MEMBER_ATLAS: members[MEMBER_ATLAS],
                             MEMBER_PDBE: members[MEMBER_PDBE]}
            try:
                envelope = self._envelope(
                    members=u_members, raw_bytes=atlas["raw"],
                    normalized={"accession": accession, "entry": entry,
                                "chain": chain, "uniprot_num": None,
                                "cell": "U", "e_state": "U",
                                "f_state": u_row["f_state"],
                                "evidence_e": False, "evidence_f": True},
                    gap_derivation=gap_derivation, modelcif_sha=None,
                    context=context)
            except ScopeFailure as exc:
                failures.append(self._failure_row(scope, exc.code, exc.detail,
                                                  exc.kind))
                continue
            rows.append({
                "row_key": _sha256_hex(
                    f"{accession}|{entry}|{chain}|gap|{u_row['ordinal']}"
                    .encode("utf-8")),
                "accession": accession, "entry": entry, "chain": chain,
                "pdb_chain": pdb_chain, "uniprot_num": None,
                "e_state": "U", "f_state": u_row["f_state"], "cell": "U",
                "e_reason": "", "f_reason": u_row["f_reason"],
                "reason": f"u_f_endpoint:{u_row['f_reason']}",
                "atlas_pdb_num": u_row["pdb_num"],
                "envelope": envelope,
            })
        # Combined rows over the union of E and F UniProt positions.
        e_map, f_map = e_out["e_map"], f_out["f_map"]
        positions = sorted(set(e_map) | set(f_map))
        e_items = [e_map.get(u, ("U", "missing_e_match")) for u in positions]
        f_items = [f_map.get(u, ("U", "missing_f_profile_position"))
                   for u in positions]
        combined = labels.combine_residue_labels(e_items, f_items)
        for index, unp in enumerate(positions):
            item = combined[index]
            e_present = unp in e_map
            f_present = unp in f_map
            try:
                envelope = self._envelope(
                    members=members, raw_bytes=primary["raw"],
                    normalized={"accession": accession, "entry": entry,
                                "chain": chain, "uniprot_num": int(unp),
                                "cell": item["cell"], "e_state": item["e_state"],
                                "f_state": item["f_state"],
                                "evidence_e": e_present, "evidence_f": f_present},
                    gap_derivation=None,
                    modelcif_sha=child["sha256"] if e_present else None,
                    context=context)
            except ScopeFailure as exc:
                failures.append(self._failure_row(scope, exc.code, exc.detail,
                                                  exc.kind))
                continue
            rows.append({
                "row_key": _sha256_hex(
                    f"{accession}|{entry}|{chain}|{unp}".encode("utf-8")),
                "accession": accession, "entry": entry, "chain": chain,
                "pdb_chain": pdb_chain, "uniprot_num": int(unp),
                "e_state": item["e_state"], "f_state": item["f_state"],
                "cell": item["cell"],
                "e_reason": item["e_reason"], "f_reason": item["f_reason"],
                "reason": item["reason"],
                "envelope": envelope,
            })
        stats = dict(e_out["stats"])
        stats.update(f_out["stats"])
        stats["rows"] = len(rows)
        return rows, failures

    # -- accession pipeline ------------------------------------------------------
    def _process_accession(self, accession: str) -> dict:
        self._scheduler.set_tag(accession)
        try:
            scopes = self._inputs["scopes"].get(accession) or []
            rows: list[dict] = []
            failures: list[dict] = []
            stats: dict[str, Any] = {}
            if not scopes:
                failures.append({"accession": accession, "entry": None,
                                 "chain": None, "pdb_chain": None,
                                 "code": LabelRunFailureCode.ADMITTED_ROW_MISSING.value,
                                 "detail": "no included census row carries the accession",
                                 "kind": "guard"})
                return self._pack(accession, rows, failures, stats)
            afdb_ref = None
            for scope in sorted(scopes, key=lambda s: (s["entry"], s["chain"])):
                ref = (scope.get("raw_refs") or {}).get(ROLE_AFDB_META)
                if isinstance(ref, dict) and ref.get("sha256"):
                    afdb_ref = ref
                    break
            afdb = load_afdb_selected(self._census_raw_dir, afdb_ref, accession,
                                      project_root=self._project_root)
            if not afdb["ok"]:
                failures.append({"accession": accession, "entry": None,
                                 "chain": None, "pdb_chain": None,
                                 "code": afdb["code"],
                                 "detail": afdb.get("detail") or "",
                                 "kind": afdb.get("kind", "guard")})
                return self._pack(accession, rows, failures, stats)
            payload = afdb["payload"]
            stats["afdb_sha256"] = afdb["sha256"]
            # Chronology (a): the selected object's cifUrl, bound per 5.1.
            child = fetch_child(self._scheduler, url=payload["cif_url"],
                                metadata_payload=afdb["payload"])
            if not child["ok"]:
                failures.append({"accession": accession, "entry": None,
                                 "chain": None, "pdb_chain": None,
                                 "code": child["code"],
                                 "detail": child.get("detail") or "",
                                 "kind": child.get("kind", "transport")})
                return self._pack(accession, rows, failures, stats)
            stats["child_sha256"] = child["sha256"]
            stats["scopes"] = len(scopes)
            stats["scopes_with_rows"] = 0
            for scope in sorted(scopes, key=lambda s: (s["entry"], s["chain"])):
                scope = dict(scope)
                scope["afdb_meta_sha256"] = afdb["sha256"]
                scope_rows, scope_failures = self._process_scope(
                    scope, afdb_payload=payload, child=child)
                rows.extend(scope_rows)
                failures.extend(scope_failures)
                if scope_rows:
                    stats["scopes_with_rows"] += 1
            return self._pack(accession, rows, failures, stats)
        except Exception as exc:  # never lose an accession to an internal error
            return {"accession": accession, "rows": [],
                    "typed_failures": [{"accession": accession, "entry": None,
                                        "chain": None, "pdb_chain": None,
                                        "code": LabelRunFailureCode.INTERNAL_ERROR.value,
                                        "detail": _clip(f"{type(exc).__name__}: {exc}"),
                                        "kind": "internal"}],
                    "receipts": [], "stats": {}}
        finally:
            self._scheduler.set_tag(None)

    def _pack(self, accession: str, rows: list, failures: list,
              stats: dict) -> dict:
        receipts = [dict(r) for r in self._scheduler.receipts
                    if r.get("accession") == accession]
        return {"accession": accession, "rows": rows, "typed_failures": failures,
                "receipts": receipts, "stats": stats}

    # -- progress / artifacts ---------------------------------------------------
    def _all_rows(self) -> list[dict]:
        out: list[dict] = []
        for accession in sorted(self._results):
            out.extend(self._results[accession]["rows"])
        return out

    def _all_failures(self) -> list[dict]:
        out: list[dict] = []
        for accession in sorted(self._results):
            out.extend(self._results[accession]["typed_failures"])
        return out

    def _ledger_obj(self, finished_utc: str | None, *, complete: bool) -> dict:
        rows = sorted(
            self._all_rows(),
            key=lambda r: (r["accession"], r["entry"], r["chain"],
                           r["uniprot_num"] if isinstance(r["uniprot_num"], int) else -1,
                           r["row_key"]))
        census_rows = [r for r in rows if isinstance(r["uniprot_num"], int)]
        cells: dict[str, int] = {c: 0 for c in labels.ALL_CELLS}
        for row in rows:
            cells[row["cell"]] = cells.get(row["cell"], 0) + 1
        payload = {
            "schema": SCHEMA,
            "emitted_by": EMITTED_BY,
            "registration": "docs/e420_registration.md",
            "label_rules": "docs/g1_successor_label_rules.json",
            "analysis_plan": "docs/g1_successor_analysis_plan.md",
            "started_utc": self._started_utc,
            "finished_utc": finished_utc,
            "inputs": {
                "pool_freeze_sha256": self._inputs["pool_freeze_sha256"],
                "pool_freeze_digest": self._inputs["pool_freeze_digest"],
                "census_manifest_sha256": self._inputs["manifest_sha256"],
            },
            "rows": rows,
            "summary": {
                "accessions_admitted": len(self._inputs["admitted"]),
                "accessions_completed": len(self._results),
                "ledger_rows": len(rows),
                "census_rows": len(census_rows),
                "u_rows_without_uniprot": len(rows) - len(census_rows),
                "cells": {c: cells.get(c, 0) for c in labels.ALL_CELLS},
                "complete": bool(complete),
            },
        }
        payload["ledger_digest"] = pool.canonical_digest(payload)
        return payload

    def _receipt_obj(self, *, complete: bool) -> dict:
        requests = _sorted_receipts(self._scheduler.receipts)
        by_role: dict[str, int] = {}
        for receipt in requests:
            by_role[receipt["role"]] = by_role.get(receipt["role"], 0) + 1
        failures = sorted(
            self._all_failures(),
            key=lambda f: (f.get("accession") or "", f.get("entry") or "",
                           f.get("chain") or "", f["code"], f.get("detail") or ""))
        rows = self._all_rows()
        census_rows = [r for r in rows if isinstance(r["uniprot_num"], int)]
        payload = {
            "schema": RECEIPT_SCHEMA,
            "emitted_by": EMITTED_BY,
            "started_utc": self._started_utc,
            "finished_utc": _utc_now(),
            "registration": {
                "registration": "docs/e420_registration.md",
                "label_rules": "docs/g1_successor_label_rules.json",
                "source_allowlist": "docs/g1_successor_source_allowlist.json",
                "analysis_plan": "docs/g1_successor_analysis_plan.md",
                "pool_freeze_sha256": self._inputs["pool_freeze_sha256"],
                "pool_freeze_digest": self._inputs["pool_freeze_digest"],
                "census_manifest_sha256": self._inputs["manifest_sha256"],
            },
            "schedule": {
                "concurrency_max": self._concurrency,
                "inter_request_delay_ms_min": self._scheduler._min_delay_s * 1000.0,
                "inter_request_delay_jitter_ms": self._scheduler._jitter_s * 1000.0,
                "backoff_initial_s": self._scheduler._backoff_initial_s,
                "backoff_factor": self._scheduler._backoff_factor,
                "max_attempts": self._scheduler._max_attempts,
                "request_timeout_s": self._scheduler._timeout_s,
            },
            "requests": requests,
            "request_counts": {"total": len(requests), "by_role": by_role,
                               "afdb_metadata_refetched": by_role.get(ROLE_AFDB_META, 0)},
            "typed_failures": failures,
            "flags": {"labels_written": True, "embeddings_read": False,
                      "models_run": 0, "e412_started": False},
            "summary": {
                "accessions_admitted": len(self._inputs["admitted"]),
                "accessions_completed": len(self._results),
                "accessions_pending": len(self._inputs["admitted"]) - len(self._results),
                "resumed_accessions": self._resumed_accessions,
                "ledger_rows": len(rows),
                "census_rows": len(census_rows),
                "u_rows_without_uniprot": len(rows) - len(census_rows),
                "typed_failure_rows": len(failures),
                "complete": bool(complete),
                "census_written": bool(complete),
                "flush_count": self._flush_count,
                "ledger_path": os.path.abspath(self._ledger_path),
                "census_path": os.path.abspath(self._census_out_path),
                "raw_dir": os.path.abspath(self._raw_dir),
            },
        }
        payload["receipt_digest"] = pool.canonical_digest(payload)
        return payload

    def _flush(self, *, final: bool, complete: bool) -> None:
        progress = {
            "schema": PROGRESS_SCHEMA,
            "results": {acc: self._results[acc] for acc in sorted(self._results)},
            "requests": _sorted_receipts(self._scheduler.receipts),
        }
        census_run._atomic_write_json(self._progress_path, progress)
        census_run._atomic_write_json(self._ledger_path,
                                      self._ledger_obj(_utc_now() if final else None,
                                                       complete=complete))
        if final:
            census_run._atomic_write_json(self._receipt_path,
                                          self._receipt_obj(complete=complete))
            if complete:
                census_rows = [{"base_accession": row["accession"],
                                "uniprot_num": row["uniprot_num"],
                                "cell": row["cell"]}
                               for row in self._all_rows()
                               if isinstance(row["uniprot_num"], int)]
                # cap_bound derives from the bound pool_freeze (registration 15.3):
                # admitted == cap => cap_bound true; never hardcoded.
                cap_bound = bool(self._inputs.get("pool_cap_bound", False))
                census = labels.four_cell_census(
                    census_rows, dict(self._inputs["group_map"]), cap_bound=cap_bound)
                census_run._atomic_write_json(self._census_out_path, census)
        self._flush_count += 1
        if self._flush_hook is not None:
            self._flush_hook(self._progress_path, self._ledger_path)

    # -- resume -----------------------------------------------------------------
    def _load_progress(self) -> set[str]:
        completed: set[str] = set()
        if not (self._resume and os.path.exists(self._progress_path)):
            return completed
        try:
            with open(self._progress_path, "r", encoding="utf-8") as handle:
                data = json.load(handle)
        except (OSError, json.JSONDecodeError):
            return completed
        if data.get("schema") != PROGRESS_SCHEMA:
            return completed
        for accession, result in (data.get("results") or {}).items():
            self._results[accession] = result
            completed.add(accession)
        self._scheduler.receipts.extend(data.get("requests") or [])
        return completed

    # -- entry point --------------------------------------------------------------
    def run(self, *, max_accessions: int | None = None) -> dict:
        admitted = list(self._inputs["admitted"])
        completed = self._load_progress()
        self._resumed_accessions = sum(1 for acc in admitted if acc in completed)
        pending = [acc for acc in admitted if acc not in completed]
        if max_accessions is not None:
            pending = pending[:max(0, int(max_accessions))]
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=self._concurrency) as executor:
            futures = [executor.submit(self._process_accession, accession)
                       for accession in pending]
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                with self._results_lock:
                    self._results[result["accession"]] = result
                    self._completed_since_flush += 1
                    due = self._completed_since_flush >= PROGRESS_EVERY_ACCESSIONS
                    if due:
                        self._completed_since_flush = 0
                if due:
                    with self._flush_lock:
                        self._flush(final=False, complete=False)
        complete = all(acc in self._results for acc in admitted)
        with self._flush_lock:
            self._flush(final=True, complete=complete)
        summary = dict(self._receipt_obj(complete=complete)["summary"])
        return summary


# --------------------------------------------------------------------------
# CLI.
# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="e420 gate (g) label construction (POOL_FREEZE sealed).")
    parser.add_argument("--pool-freeze", default=os.path.join(ROOT, POOL_FREEZE_REL))
    parser.add_argument("--manifest", default=os.path.join(ROOT, MANIFEST_REL))
    parser.add_argument("--census-raw-dir", default=os.path.join(ROOT, CENSUS_RAW_REL))
    parser.add_argument("--ledger", default=os.path.join(ROOT, LEDGER_REL))
    parser.add_argument("--census-out", default=os.path.join(ROOT, CENSUS_OUT_REL))
    parser.add_argument("--receipt", default=os.path.join(ROOT, RECEIPT_REL))
    parser.add_argument("--raw-dir", default=os.path.join(ROOT, LABELS_RAW_REL))
    parser.add_argument("--resume", action="store_true",
                        help="skip accessions completed in the progress file")
    parser.add_argument("--max-accessions", type=int, default=None,
                        help="bound this invocation to the next N pending accessions")
    parser.add_argument("--project-root", default=ROOT)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--workers", type=int, default=CONCURRENCY_MAX)
    parser.add_argument("--min-delay-ms", type=float,
                        default=census_run.INTER_REQUEST_DELAY_MS_MIN)
    parser.add_argument("--jitter-ms", type=float,
                        default=census_run.INTER_REQUEST_JITTER_MS)
    parser.add_argument("--backoff-initial-s", type=float,
                        default=census_run.BACKOFF_INITIAL_S)
    parser.add_argument("--max-attempts", type=int, default=census_run.MAX_ATTEMPTS)
    parser.add_argument("--registered-schedule-override", action="store_true",
                        help="required when deviating from the registered schedule (audit flag)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    min_delay_s = args.min_delay_ms / 1000.0
    jitter_s = args.jitter_ms / 1000.0
    registered = (min_delay_s == census_run.INTER_REQUEST_DELAY_MS_MIN / 1000.0
                  and jitter_s == census_run.INTER_REQUEST_JITTER_MS / 1000.0
                  and args.backoff_initial_s == census_run.BACKOFF_INITIAL_S
                  and args.max_attempts == census_run.MAX_ATTEMPTS
                  and args.workers <= CONCURRENCY_MAX)
    if not registered and not args.registered_schedule_override:
        print("refusing: schedule deviates from the registered bound without "
              "--registered-schedule-override", file=sys.stderr)
        return 2
    try:
        runner = LabelRunRunner(
            ledger_path=args.ledger, census_out_path=args.census_out,
            receipt_path=args.receipt, raw_dir=args.raw_dir,
            project_root=args.project_root, pool_freeze_path=args.pool_freeze,
            manifest_path=args.manifest, census_raw_dir=args.census_raw_dir,
            resume=args.resume, concurrency=args.workers,
            min_delay_s=min_delay_s, jitter_s=jitter_s,
            backoff_initial_s=args.backoff_initial_s,
            max_attempts=args.max_attempts, seed=args.seed)
        summary = runner.run(max_accessions=args.max_accessions)
    except LabelRunError as exc:
        print(json.dumps({"typed_error": exc.code, "detail": exc.detail},
                         sort_keys=True, separators=(",", ":")), file=sys.stderr)
        return 3
    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
