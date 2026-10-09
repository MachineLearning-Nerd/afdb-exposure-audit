#!/usr/bin/env python3
"""e420 gate-(e) phase (iii): predictor-blind availability-census runner.

Registration: docs/e420_registration.md sections 4.1/4.2/10(e)-(f)/12.
Schedule: docs/g1_successor_analysis_plan.md section 5. Endpoints/roles:
docs/g1_successor_source_allowlist.json.

Pre-registered behavior:
  * candidate universe = the retained census member only; the zip and member
    SHA-256 pins are verified through the phase-(i) validator
    (e420_grammar_validate.validate_atlas_census_zip) BEFORE any network
    activity;
  * one segment GET per unique entry (1,735 expected); per row the row's chain
    is resolved to (entity_id, accession) via the frozen O3 route; one primary
    GET per (entry, entity); one AFDB prediction GET per unique accession;
  * ThreadPoolExecutor max_workers=2; a global-lock inter-start gate enforces
    >= 250 ms with +/- 100 ms jitter between request starts; on 429/5xx the
    runner backs off 2 s doubling, max 3 attempts, then emits a typed failure
    row; a 429 Retry-After is honored (capped at 600 s so no retry can
    violate the bounded-timeout discipline);
  * every response is receipted (url, method, status, headers, byte length,
    sha256, UTC) and its raw bytes are saved under --raw-dir with hash
    filenames;
  * per-row typed availability record per registration 4.1 (opaque row key,
    entry/chain, source roles touched, status, raw sha refs, accession(s),
    entity_id, inclusion/exclusion status, typed failure code + detail). NO
    predictor/coordinate/RMSF/B-factor/pLDDT values are read; AFDB metadata is
    parsed ONLY by the phase-(i) validator (identity/temporal fields;
    fraction_plddt* fields are hash-without-parse);
  * resume: a progress file (rows completed + raw refs) lives inside
    --raw-dir and is rewritten atomically every 25 rows; --resume skips
    completed rows;
  * output: the manifest JSON (sorted keys, compact) is rewritten
    incrementally and carries a summary block (rows attempted/succeeded/
    typed-failed, unique entries/accessions, request count).

STRICT write scope: nothing is written outside --raw-dir and --manifest.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import hashlib
import io
import json
import os
import random
import sys
import tempfile
import threading
import time
import zipfile
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import requests  # noqa: E402  (default transport only; tests inject fakes)

import e420_grammar_validate as gv  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --------------------------------------------------------------------------
# Registered pins (registration 3.4 / analysis plan section 2).
# --------------------------------------------------------------------------
CENSUS_ZIP_REL = os.path.join("data", "source_contracts", "atlas_parsable_full_20260830.zip")
CENSUS_MEMBER_PATH = gv.ATLAS_CENSUS_MEMBER_PATH
CENSUS_ZIP_SHA256 = gv.ATLAS_CENSUS_ZIP_SHA256
CENSUS_MEMBER_SHA256 = gv.ATLAS_CENSUS_MEMBER_SHA256

SCHEMA = "e420-census-availability-manifest-v1"
EMITTED_BY = "e420-census-v1"
PROGRESS_SCHEMA = "e420-census-progress-v1"
PROGRESS_NAME = "e420_census_progress.json"

# --------------------------------------------------------------------------
# Registered request schedule (analysis plan section 5 / allowlist
# census_request_schedule).  Defaults are the registered conservative bound.
# --------------------------------------------------------------------------
CONCURRENCY_MAX = 2
INTER_REQUEST_DELAY_MS_MIN = 250.0
INTER_REQUEST_JITTER_MS = 100.0
BACKOFF_INITIAL_S = 2.0
BACKOFF_FACTOR = 2.0
MAX_ATTEMPTS = 3
REQUEST_TIMEOUT_S = 60.0
RETRY_AFTER_CAP_S = 600.0
PROGRESS_EVERY_ROWS = 25
USER_AGENT = "moluq-e420-census/1.0"

# Registered endpoints (allowlist source_roles; no other route is permitted).
PDBE_SEGMENT_URL = "https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}"
PDBE_PRIMARY_URL = "https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/uniprot_mapping/{pdb_id}/{entity_id}"
AFDB_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

ROLE_SEGMENT = gv.SourceRole.PDBE_MAPPING_SEGMENT.value
ROLE_PRIMARY = gv.SourceRole.PDBE_MAPPING_PRIMARY.value
ROLE_AFDB = gv.SourceRole.AFDB_METADATA.value

# AFDB identity/temporal fields retained per section 4.1 + parse whitelist.
AFDB_IDENTITY_FIELDS = (
    "model_entity_id",
    "model_created_date",
    "sequence_version_date",
    "latest_version",
    "all_versions",
    "cif_url",
    "uniprot_accession",
)

ROW_FIELDS = frozenset({
    "row_key",          # opaque candidate key (sha256 of the census row)
    "entry",
    "chain",
    "source_roles_touched",
    "status",           # "ok" | "typed_failure"
    "raw_refs",         # role -> {url,status,sha256,bytes_len}
    "accessions",       # resolved canonical accession(s)
    "entity_id",
    "availability",     # {"pdbe_mapping": bool, "afdb_prediction": bool}
    "inclusion_status",  # "included" | "excluded" (availability-level only)
    "exclusion_reason",
    "failure_code",
    "failure_detail",
    "afdb_identity",    # identity/temporal fields from the phase-(i) validator
})

RECEIPT_FIELDS = frozenset({
    "url", "method", "status", "headers", "bytes_len", "sha256", "utc",
    "role", "attempt", "raw_name", "retry_after_s", "next_backoff_s", "error",
})


class CensusFailureCode(str, Enum):
    """Census-level typed failure codes (grammar codes come from the
    phase-(i) validator and are passed through verbatim)."""

    TRANSPORT_ERROR = "E420_TRANSPORT_ERROR"
    HTTP_429_EXHAUSTED = "E420_HTTP_429_MAX_ATTEMPTS"
    HTTP_5XX_EXHAUSTED = "E420_HTTP_5XX_MAX_ATTEMPTS"
    HTTP_UNEXPECTED_STATUS = "E420_HTTP_UNEXPECTED_STATUS"
    SEGMENT_UNAVAILABLE = "E420_SEGMENT_UNAVAILABLE"
    PRIMARY_UNAVAILABLE = "E420_PRIMARY_UNAVAILABLE"
    AFDB_PREDICTION_ABSENT = "E420_AFDB_PREDICTION_ABSENT"
    CHAIN_UNRESOLVED = "E420_CHAIN_UNRESOLVED"
    CHAIN_AMBIGUOUS = "E420_CHAIN_AMBIGUOUS"
    PRIMARY_ACCESSION_MISMATCH = "E420_PRIMARY_ACCESSION_MISMATCH"
    SEGMENT_THRESHOLD_FAILURE = "E420_RESOLUTION_THRESHOLD_FAILURE"
    INTERNAL_ERROR = "E420_INTERNAL_ERROR"


class CensusError(RuntimeError):
    """Typed pre-network setup error (census universe, pins, arguments)."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


def _utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _clip(text: Any, limit: int = 300) -> str:
    return "" if text is None else str(text)[:limit]


# --------------------------------------------------------------------------
# Transport abstraction (injectable; tests use offline fakes).
# --------------------------------------------------------------------------
class HttpResponse:
    __slots__ = ("status", "headers", "content")

    def __init__(self, status: int, headers: dict, content: bytes) -> None:
        self.status = status
        self.headers = headers
        self.content = content


class Transport:
    def get(self, url: str, timeout_s: float) -> HttpResponse:  # pragma: no cover
        raise NotImplementedError


class RequestsTransport(Transport):
    """Default GET-only transport (public endpoints, no credentials)."""

    def __init__(self, user_agent: str = USER_AGENT) -> None:
        self._user_agent = user_agent

    def get(self, url: str, timeout_s: float) -> HttpResponse:
        r = requests.get(url, headers={"User-Agent": self._user_agent}, timeout=timeout_s)
        return HttpResponse(r.status_code, dict(r.headers), r.content)


@dataclass
class RequestResult:
    """Final outcome of one scheduled request (all attempts receipted)."""

    role: str
    url: str
    ok: bool = False
    status: int | None = None
    sha256: str | None = None
    bytes_len: int | None = None
    content: bytes | None = None
    raw_name: str | None = None
    receipts: list = field(default_factory=list)
    failure_code: str | None = None
    failure_detail: str = ""

    @property
    def raw_ref(self) -> dict | None:
        if self.sha256 is None:
            return None
        return {"url": self.url, "status": self.status, "sha256": self.sha256,
                "bytes_len": self.bytes_len}


def _retry_after_seconds(headers: Any) -> float | None:
    if not isinstance(headers, dict):
        return None
    value = None
    for key in ("Retry-After", "retry-after"):
        if key in headers:
            value = headers[key]
            break
    if value is None:
        return None
    text = str(value).strip()
    try:
        return max(0.0, float(text))
    except ValueError:
        pass
    for fmt in ("%a, %b %d %H:%M:%S %Y", "%A, %d-%b-%y %H:%M:%S %Z"):
        try:
            dt = datetime.datetime.strptime(text, fmt)
            delta = (dt.replace(tzinfo=datetime.timezone.utc)
                     - datetime.datetime.now(datetime.timezone.utc)).total_seconds()
            return max(0.0, delta)
        except ValueError:
            continue
    return None


class RequestScheduler:
    """Global-lock inter-start gate + registered 429/5xx retry discipline.

    Concurrency, delay floor and backoff come from the registered schedule;
    every response (including each failed attempt) is receipted and its raw
    bytes are persisted under the raw dir with hash filenames.
    """

    def __init__(
        self,
        transport: Transport,
        raw_dir: str,
        *,
        min_delay_s: float = INTER_REQUEST_DELAY_MS_MIN / 1000.0,
        jitter_s: float = INTER_REQUEST_JITTER_MS / 1000.0,
        backoff_initial_s: float = BACKOFF_INITIAL_S,
        backoff_factor: float = BACKOFF_FACTOR,
        max_attempts: int = MAX_ATTEMPTS,
        timeout_s: float = REQUEST_TIMEOUT_S,
        seed: int = 20260904,
        sleeper: Callable[[float], None] = time.sleep,
    ) -> None:
        self._transport = transport
        self._raw_dir = raw_dir
        self._min_delay_s = float(min_delay_s)
        self._jitter_s = float(jitter_s)
        self._backoff_initial_s = float(backoff_initial_s)
        self._backoff_factor = float(backoff_factor)
        self._max_attempts = int(max_attempts)
        self._timeout_s = float(timeout_s)
        self._sleeper = sleeper
        self._rng = random.Random(seed)
        self._gate_lock = threading.Lock()
        self._next_start: float | None = None
        self._receipts_lock = threading.Lock()
        self.receipts: list[dict] = []

    # -- schedule -----------------------------------------------------------
    def _sample_delay(self) -> float:
        if self._min_delay_s <= 0.0 and self._jitter_s <= 0.0:
            return 0.0
        return max(0.0, self._min_delay_s + self._rng.uniform(-self._jitter_s, self._jitter_s))

    def _enter_gate(self) -> None:
        """Serialize request STARTS: >= sampled delay between consecutive
        starts across all worker threads (global lock)."""
        with self._gate_lock:
            now = time.monotonic()
            wait = 0.0 if self._next_start is None else self._next_start - now
            if wait > 0.0:
                self._sleeper(wait)
            self._next_start = time.monotonic() + self._sample_delay()

    # -- raw persistence ----------------------------------------------------
    def _save_raw(self, content: bytes) -> str:
        name = _sha256_hex(content)
        path = os.path.join(self._raw_dir, name)
        if not os.path.exists(path):
            _atomic_write_bytes(path, content)
        return name

    def _add_receipt(self, receipt: dict) -> None:
        with self._receipts_lock:
            self.receipts.append(receipt)

    # -- request ------------------------------------------------------------
    def get(self, role: gv.SourceRole, url: str) -> RequestResult:
        result = RequestResult(role=role.value, url=url)
        for attempt in range(1, self._max_attempts + 1):
            self._enter_gate()
            try:
                resp = self._transport.get(url, self._timeout_s)
            except Exception as exc:  # transport-level error: typed, no retry
                receipt = {"url": url, "method": "GET", "status": None, "headers": {},
                           "bytes_len": 0, "sha256": None, "utc": _utc_now(),
                           "role": role.value, "attempt": attempt, "raw_name": None,
                           "retry_after_s": None, "next_backoff_s": None,
                           "error": type(exc).__name__}
                self._add_receipt(receipt)
                result.receipts.append(receipt)
                result.failure_code = CensusFailureCode.TRANSPORT_ERROR.value
                result.failure_detail = _clip(f"{type(exc).__name__}: {exc}")
                return result
            content = resp.content
            sha = _sha256_hex(content)
            raw_name = self._save_raw(content)
            receipt = {"url": url, "method": "GET", "status": int(resp.status),
                       "headers": dict(resp.headers), "bytes_len": len(content),
                       "sha256": sha, "utc": _utc_now(), "role": role.value,
                       "attempt": attempt, "raw_name": raw_name,
                       "retry_after_s": None, "next_backoff_s": None, "error": None}
            self._add_receipt(receipt)
            result.receipts.append(receipt)
            result.status = int(resp.status)
            result.sha256 = sha
            result.bytes_len = len(content)
            result.content = content
            result.raw_name = raw_name
            if resp.status == 429 or 500 <= resp.status < 600:
                retry_after = _retry_after_seconds(resp.headers)
                if retry_after is not None and retry_after > RETRY_AFTER_CAP_S:
                    retry_after = RETRY_AFTER_CAP_S
                receipt["retry_after_s"] = retry_after
                if attempt >= self._max_attempts:
                    result.failure_code = (CensusFailureCode.HTTP_429_EXHAUSTED.value
                                           if resp.status == 429
                                           else CensusFailureCode.HTTP_5XX_EXHAUSTED.value)
                    result.failure_detail = _clip(f"status={resp.status} after {attempt} attempts")
                    return result
                delay = self._backoff_initial_s * (self._backoff_factor ** (attempt - 1))
                if retry_after is not None:
                    delay = max(delay, retry_after)
                receipt["next_backoff_s"] = delay
                self._sleeper(delay)  # backoff outside the start gate
                continue
            if not (200 <= resp.status < 300):
                # e423 repair (DECISIONS 2026-09-10): a non-2xx status
                # outside the 429/5xx retry class (400/404/...) was marked
                # ok and its error body parsed as data — the silent-drop
                # defect the e422 dispatch shakedown caught.  Typed here.
                result.failure_code = CensusFailureCode.HTTP_UNEXPECTED_STATUS.value
                result.failure_detail = _clip(
                    f"status={resp.status} (non-retryable); body head: "
                    + content[:200].decode("utf-8", errors="replace"))
                return result
            result.ok = True
            return result
        return result  # pragma: no cover - loop always returns


# --------------------------------------------------------------------------
# Atomic IO helpers.
# --------------------------------------------------------------------------
def _atomic_write_bytes(path: str, data: bytes) -> None:
    dirp = os.path.dirname(os.path.abspath(path))
    os.makedirs(dirp, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dirp, prefix=".tmp_e420_", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _atomic_write_json(path: str, obj: Any) -> None:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8") + b"\n"
    _atomic_write_bytes(path, payload)


# --------------------------------------------------------------------------
# Census universe loading (phase-(i) validator enforces the pins).
# --------------------------------------------------------------------------
def load_census_rows(
    project_root: str,
    *,
    zip_path: str | None = None,
    zip_sha256: str | None = None,
    member_sha256: str | None = None,
    member_path: str | None = None,
) -> tuple[list[str], dict]:
    """Load + verify the retained census member; return sorted rows and pins."""
    zip_path = zip_path if zip_path is not None else os.path.join(project_root, CENSUS_ZIP_REL)
    zip_sha256 = zip_sha256 if zip_sha256 is not None else CENSUS_ZIP_SHA256
    member_sha256 = member_sha256 if member_sha256 is not None else CENSUS_MEMBER_SHA256
    member_path = member_path if member_path is not None else CENSUS_MEMBER_PATH
    with open(zip_path, "rb") as handle:
        raw = handle.read()
    outcome = gv.validate_atlas_census_zip(
        raw, expected_zip_sha256=zip_sha256, expected_member_sha256=member_sha256,
        member_path=member_path)
    if not outcome.ok:
        failure = outcome.failure
        raise CensusError("E420_CENSUS_UNIVERSE_INVALID",
                          f"{failure.code.value}: {_clip(failure.detail)}")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        member = archive.read(member_path)
    rows = sorted({line for line in member.decode("ascii").split("\n") if line})
    pins = {
        "census_zip_sha256": zip_sha256,
        "census_member_sha256": member_sha256,
        "census_member_path": member_path,
        "census_row_count": int(outcome.payload["row_count"]),
        "census_unique_entry_count": int(outcome.payload["unique_entry_count"]),
    }
    return rows, pins


# --------------------------------------------------------------------------
# Chain -> (entity_id, accession) resolution via the frozen O3 route.
# --------------------------------------------------------------------------
def scan_segment_identity(raw: bytes, pdb_id: str, chain: str) -> list[tuple[str, int]]:
    """Allowlist-bounded identity scan of the segment envelope: accession map
    keys + per-mapping chain_id/entity_id ONLY (identity/coverage numerics are
    read later, only by the phase-(i) threshold check)."""
    data = json.loads(raw.decode("utf-8"))
    entry = data.get(pdb_id) if isinstance(data, dict) else None
    uniprot = entry.get("UniProt") if isinstance(entry, dict) else None
    matches: list[tuple[str, int]] = []
    if isinstance(uniprot, dict):
        for accession, obj in uniprot.items():
            if not isinstance(obj, dict):
                continue
            mappings = obj.get("mappings")
            if not isinstance(mappings, list):
                continue
            for mapping in mappings:
                if not isinstance(mapping, dict):
                    continue
                if mapping.get("chain_id") == chain and isinstance(mapping.get("entity_id"), int):
                    matches.append((str(accession), int(mapping["entity_id"])))
    return matches


def resolve_chain(raw: bytes, pdb_id: str, chain: str, project_root: str) -> dict:
    """Resolve one row's chain to (accession, entity_id) + threshold guard.

    All envelope grammar passes through the phase-(i) validators
    (validate_pdbe_segment / check_pdbe_segment_thresholds).
    """
    try:
        matches = scan_segment_identity(raw, pdb_id, chain)
    except Exception:
        # Route the raw bytes through the validator to obtain the typed code.
        probe = gv.validate_pdbe_segment(raw, pdb_id=pdb_id, candidate_accession="SCANPROBE",
                                         project_root=project_root)
        code = probe.failure.code.value if probe.failure else gv.GrammarFailureCode.BOUNDARY_INVALID.value
        detail = probe.failure.detail if probe.failure else "segment identity scan failed"
        return {"ok": False, "kind": "grammar", "accession": None, "entity_id": None,
                "failure_code": code, "failure_detail": _clip(detail)}
    accessions = sorted({a for a, _ in matches})
    entities = sorted({e for _, e in matches})
    if not matches:
        return {"ok": False, "kind": "unresolved", "accession": None, "entity_id": None,
                "failure_code": CensusFailureCode.CHAIN_UNRESOLVED.value,
                "failure_detail": _clip(f"chain={chain!r} absent from segment mappings")}
    if len(accessions) > 1 or len(entities) > 1:
        return {"ok": False, "kind": "ambiguous", "accession": None, "entity_id": None,
                "failure_code": CensusFailureCode.CHAIN_AMBIGUOUS.value,
                "failure_detail": _clip(f"accessions={accessions}:entities={entities}")}
    accession, entity_id = accessions[0], entities[0]
    outcome = gv.validate_pdbe_segment(raw, pdb_id=pdb_id, candidate_accession=accession,
                                       project_root=project_root)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "accession": accession, "entity_id": entity_id,
                "failure_code": outcome.failure.code.value,
                "failure_detail": _clip(outcome.failure.detail)}
    thresholds = gv.check_pdbe_segment_thresholds(
        raw, pdb_id=pdb_id, candidate_accession=accession, entity_id=entity_id,
        mode=gv.SegmentMode.RESOLUTION, project_root=project_root)
    if not thresholds.ok:
        return {"ok": False, "kind": "threshold", "accession": accession, "entity_id": entity_id,
                "failure_code": CensusFailureCode.SEGMENT_THRESHOLD_FAILURE.value,
                "failure_detail": _clip(f"{thresholds.failure.code.value}: {thresholds.failure.detail}")}
    return {"ok": True, "kind": "ok", "accession": accession, "entity_id": entity_id,
            "failure_code": None, "failure_detail": ""}


# --------------------------------------------------------------------------
# Per-role fetch+validate steps (validator-first; no value reads beyond the
# registered census parse whitelist).
# --------------------------------------------------------------------------
def fetch_primary(scheduler: RequestScheduler, entry: str, entity_id: int,
                  project_root: str) -> dict:
    url = PDBE_PRIMARY_URL.format(pdb_id=entry, entity_id=entity_id)
    res = scheduler.get(gv.SourceRole.PDBE_MAPPING_PRIMARY, url)
    # e423 repair: classify 404 directly on the status — under the repaired
    # scheduler a 404 is ok=False (it used to sneak through as ok), so the
    # 404-specific "unavailable" semantics must be tested before the hard
    # failure gate.  Transport errors (status None) still fall to "hard".
    if res.status == 404:
        return {"ok": False, "kind": "unavailable", "status": res.status, "accessions": [],
                "raw_ref": res.raw_ref,
                "failure_code": CensusFailureCode.PRIMARY_UNAVAILABLE.value,
                "failure_detail": _clip("status=404")}
    if not res.ok:
        return {"ok": False, "kind": "hard", "status": res.status, "accessions": [],
                "raw_ref": res.raw_ref, "failure_code": res.failure_code,
                "failure_detail": res.failure_detail}
    if res.status != 200:
        return {"ok": False, "kind": "unavailable", "status": res.status, "accessions": [],
                "raw_ref": res.raw_ref,
                "failure_code": CensusFailureCode.HTTP_UNEXPECTED_STATUS.value,
                "failure_detail": _clip(f"status={res.status}")}
    outcome = gv.validate_pdbe_primary(res.content, pdb_id=entry, entity_id=entity_id,
                                       project_root=project_root)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "status": res.status, "accessions": [],
                "raw_ref": res.raw_ref, "failure_code": outcome.failure.code.value,
                "failure_detail": _clip(outcome.failure.detail)}
    accessions = sorted({str(item["accession"]) for item in outcome.payload.get("entries", [])
                         if isinstance(item, dict) and item.get("accession")})
    return {"ok": True, "kind": "ok", "status": res.status, "accessions": accessions,
            "raw_ref": res.raw_ref, "failure_code": None, "failure_detail": ""}


def fetch_afdb(scheduler: RequestScheduler, accession: str, project_root: str) -> dict:
    url = AFDB_PREDICTION_URL.format(accession=accession)
    res = scheduler.get(gv.SourceRole.AFDB_METADATA, url)
    # e423 repair: 404 classified directly on the status (see fetch_primary)
    if res.status == 404:
        return {"ok": False, "kind": "absent", "status": res.status, "identity": None,
                "raw_ref": res.raw_ref,
                "failure_code": CensusFailureCode.AFDB_PREDICTION_ABSENT.value,
                "failure_detail": _clip("status=404")}
    if not res.ok:
        return {"ok": False, "kind": "hard", "status": res.status, "identity": None,
                "raw_ref": res.raw_ref, "failure_code": res.failure_code,
                "failure_detail": res.failure_detail}
    if res.status != 200:
        return {"ok": False, "kind": "unavailable", "status": res.status, "identity": None,
                "raw_ref": res.raw_ref,
                "failure_code": CensusFailureCode.HTTP_UNEXPECTED_STATUS.value,
                "failure_detail": _clip(f"status={res.status}")}
    outcome = gv.validate_afdb_metadata(res.content, project_root=project_root,
                                        expected_accession=accession)
    if not outcome.ok:
        return {"ok": False, "kind": "grammar", "status": res.status, "identity": None,
                "raw_ref": res.raw_ref, "failure_code": outcome.failure.code.value,
                "failure_detail": _clip(outcome.failure.detail)}
    payload = outcome.payload
    identity = {name: payload[name] for name in AFDB_IDENTITY_FIELDS if name in payload}
    return {"ok": True, "kind": "ok", "status": res.status, "identity": identity,
            "raw_ref": res.raw_ref, "failure_code": None, "failure_detail": ""}


class AfdbCache:
    """Global one-GET-per-accession cache (thread-safe fetch-once)."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._events: dict[str, threading.Event] = {}
        self._results: dict[str, dict | None] = {}

    def preload(self, accession: str, result: dict) -> None:
        event = threading.Event()
        event.set()
        with self._lock:
            self._events[accession] = event
            self._results[accession] = result

    def get(self, accession: str, fetch: Callable[[], dict]) -> dict:
        with self._lock:
            if accession in self._results:
                return self._results[accession]
            if accession in self._events:
                event, owner = self._events[accession], False
            else:
                event = threading.Event()
                self._events[accession] = event
                owner = True
        if not owner:
            event.wait()
            with self._lock:
                return self._results.get(accession)
        try:
            result = fetch()
        except BaseException as exc:
            result = {"ok": False, "kind": "hard", "status": None, "identity": None,
                      "raw_ref": None, "failure_code": CensusFailureCode.INTERNAL_ERROR.value,
                      "failure_detail": _clip(f"{type(exc).__name__}: {exc}")}
        finally:
            with self._lock:
                self._results[accession] = result
            event.set()
        return result


# --------------------------------------------------------------------------
# The runner.
# --------------------------------------------------------------------------
class CensusRunner:
    def __init__(
        self,
        *,
        manifest_path: str,
        raw_dir: str,
        project_root: str = ROOT,
        transport: Transport | None = None,
        resume: bool = False,
        concurrency: int = CONCURRENCY_MAX,
        min_delay_s: float = INTER_REQUEST_DELAY_MS_MIN / 1000.0,
        jitter_s: float = INTER_REQUEST_JITTER_MS / 1000.0,
        backoff_initial_s: float = BACKOFF_INITIAL_S,
        backoff_factor: float = BACKOFF_FACTOR,
        max_attempts: int = MAX_ATTEMPTS,
        timeout_s: float = REQUEST_TIMEOUT_S,
        seed: int = 20260904,
        sleeper: Callable[[float], None] = time.sleep,
        flush_hook: Callable[[str, str], None] | None = None,
    ) -> None:
        if concurrency > CONCURRENCY_MAX:
            raise CensusError("E420_CONCURRENCY_ABOVE_REGISTERED_BOUND",
                              f"concurrency={concurrency} > {CONCURRENCY_MAX}")
        if max_attempts < 1:
            raise CensusError("E420_INVALID_ARGUMENT", "max_attempts < 1")
        self._manifest_path = manifest_path
        self._raw_dir = raw_dir
        self._project_root = project_root
        self._resume = bool(resume)
        self._concurrency = int(concurrency)
        self._seed = int(seed)
        self._flush_hook = flush_hook
        os.makedirs(self._raw_dir, exist_ok=True)
        self._progress_path = os.path.join(self._raw_dir, PROGRESS_NAME)
        self._scheduler = RequestScheduler(
            transport if transport is not None else RequestsTransport(),
            self._raw_dir,
            min_delay_s=min_delay_s, jitter_s=jitter_s,
            backoff_initial_s=backoff_initial_s, backoff_factor=backoff_factor,
            max_attempts=max_attempts, timeout_s=timeout_s, seed=seed, sleeper=sleeper)
        self._records: dict[str, dict] = {}
        self._records_lock = threading.Lock()
        self._flush_lock = threading.Lock()
        self._completed_since_flush = 0
        self._flush_count = 0
        self._started_utc = _utc_now()
        self._base_pins: dict = {}
        self._schedule_pins: dict | None = None
        self._primary_known: dict[tuple[str, int], dict] = {}
        self._afdb = AfdbCache()

    # -- progress / manifest ------------------------------------------------
    def _sorted_receipts(self) -> list[dict]:
        def key(r: dict):
            sha = r.get("sha256") or ""
            status = r["status"] if isinstance(r.get("status"), int) else -1
            return (r["role"], r["url"], int(r["attempt"]), sha, status, int(r["bytes_len"]))
        return sorted(self._scheduler.receipts, key=key)

    def _manifest_obj(self, finished_utc: str | None) -> dict:
        rows = sorted(self._records.values(), key=lambda r: r["row_key"])
        requests = self._sorted_receipts()
        accessions = {acc for row in rows for acc in row["accessions"]}
        summary = {
            "rows_attempted": len(rows),
            "rows_succeeded": sum(1 for row in rows if row["status"] == "ok"),
            "rows_typed_failed": sum(1 for row in rows if row["status"] == "typed_failure"),
            "unique_entries": len({row["entry"] for row in rows}),
            "unique_accessions": len(accessions),
            "request_count": len(requests),
        }
        registration = dict(self._base_pins)
        if self._schedule_pins is not None:
            registration["schedule"] = dict(self._schedule_pins)
        return {
            "schema": SCHEMA,
            "emitted_by": EMITTED_BY,
            "started_utc": self._started_utc,
            "finished_utc": finished_utc,
            "registration": registration,
            "rows": rows,
            "requests": requests,
            "summary": summary,
        }

    def _flush(self, *, final: bool) -> None:
        progress = {
            "schema": PROGRESS_SCHEMA,
            "rows": sorted(self._records.values(), key=lambda r: r["row_key"]),
            "requests": self._sorted_receipts(),
        }
        _atomic_write_json(self._progress_path, progress)
        _atomic_write_json(self._manifest_path,
                           self._manifest_obj(_utc_now() if final else None))
        self._flush_count += 1
        if self._flush_hook is not None:
            self._flush_hook(self._progress_path, self._manifest_path)

    def _maybe_flush(self) -> None:
        with self._flush_lock:
            if self._completed_since_flush >= PROGRESS_EVERY_ROWS:
                self._completed_since_flush = 0
                with self._records_lock:
                    self._flush(final=False)

    # -- resume -------------------------------------------------------------
    def _load_progress(self) -> set[str]:
        completed: set[str] = set()
        if not (self._resume and os.path.exists(self._progress_path)):
            return completed
        with open(self._progress_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        for record in data.get("rows", []):
            self._records[record["row_key"]] = record
            completed.add(record["row_key"])
            self._seed_caches(record)
        self._scheduler.receipts.extend(data.get("requests", []))
        return completed

    def _seed_caches(self, record: dict) -> None:
        availability = record.get("availability") or {}
        raw_refs = record.get("raw_refs") or {}
        if (availability.get("afdb_prediction") and record.get("afdb_identity")
                and ROLE_AFDB in raw_refs and record["accessions"]):
            self._afdb.preload(record["accessions"][0], {
                "ok": True, "kind": "ok",
                "status": raw_refs[ROLE_AFDB].get("status", 200),
                "identity": dict(record["afdb_identity"]),
                "raw_ref": dict(raw_refs[ROLE_AFDB]),
                "failure_code": None, "failure_detail": ""})

    # -- per-entry pipeline -------------------------------------------------
    def _process_entry(self, entry: str, chains: list[str]) -> None:
        records: dict[str, dict] = {}
        try:
            seg_res = self._scheduler.get(
                gv.SourceRole.PDBE_MAPPING_SEGMENT,
                PDBE_SEGMENT_URL.format(pdb_id=entry))
            seg_ref = seg_res.raw_ref
            if seg_res.ok and seg_res.status == 200:
                seg_stage = {"kind": "ok", "failure_code": None, "failure_detail": ""}
                raw = seg_res.content
                for chain in chains:
                    row_key = _row_key(entry, chain)
                    records[row_key] = self._process_row(
                        entry, chain, row_key, raw, seg_ref, seg_stage)
            elif seg_res.status == 404:
                # e423 repair: 404 classified directly on the status
                seg_stage = {"kind": "unavailable",
                             "failure_code": CensusFailureCode.SEGMENT_UNAVAILABLE.value,
                             "failure_detail": _clip("status=404")}
            elif not seg_res.ok:
                seg_stage = {"kind": "hard", "failure_code": seg_res.failure_code,
                             "failure_detail": seg_res.failure_detail}
            else:
                seg_stage = {"kind": "unavailable",
                             "failure_code": CensusFailureCode.HTTP_UNEXPECTED_STATUS.value,
                             "failure_detail": _clip(f"status={seg_res.status}")}
            for chain in chains:  # entry-level failure paths -> typed rows
                row_key = _row_key(entry, chain)
                if row_key not in records:
                    records[row_key] = self._entry_failure_record(
                        entry, chain, row_key, seg_ref, seg_stage)
        except Exception as exc:  # never lose a row to an internal error
            detail = _clip(f"{type(exc).__name__}: {exc}")
            for chain in chains:
                row_key = _row_key(entry, chain)
                records[row_key] = self._failure_record(
                    entry, chain, row_key, None, [ROLE_SEGMENT],
                    CensusFailureCode.INTERNAL_ERROR.value, detail)
        with self._records_lock:
            self._records.update(records)
            self._completed_since_flush += len(records)
        self._maybe_flush()

    def _process_row(self, entry: str, chain: str, row_key: str, raw: bytes,
                     seg_ref: dict | None, seg_stage: dict) -> dict:
        roles = [ROLE_SEGMENT]
        raw_refs: dict[str, dict] = {}
        if seg_ref is not None:
            raw_refs[ROLE_SEGMENT] = seg_ref
        resolution = resolve_chain(raw, entry, chain, self._project_root)
        accession = resolution["accession"]
        entity_id = resolution["entity_id"]
        primary_out: dict | None = None
        afdb_out: dict | None = None
        if resolution["ok"]:
            key = (entry, entity_id)
            with self._records_lock:
                primary_out = self._primary_known.get(key)
            if primary_out is None:
                primary_out = fetch_primary(self._scheduler, entry, entity_id,
                                            self._project_root)
                with self._records_lock:
                    self._primary_known[key] = primary_out
            roles.append(ROLE_PRIMARY)
            if primary_out.get("raw_ref"):
                raw_refs[ROLE_PRIMARY] = primary_out["raw_ref"]
            if primary_out["ok"]:
                afdb_out = self._afdb.get(
                    accession,
                    lambda acc=accession: fetch_afdb(self._scheduler, acc, self._project_root))
                roles.append(ROLE_AFDB)
                if afdb_out.get("raw_ref"):
                    raw_refs[ROLE_AFDB] = afdb_out["raw_ref"]
        return self._assemble_row(entry, chain, row_key, roles, raw_refs,
                                  resolution, primary_out, afdb_out)

    def _entry_failure_record(self, entry: str, chain: str, row_key: str,
                              seg_ref: dict | None, seg_stage: dict) -> dict:
        raw_refs = {ROLE_SEGMENT: seg_ref} if seg_ref is not None else {}
        resolution = {"ok": False, "kind": seg_stage["kind"], "accession": None,
                      "entity_id": None, "failure_code": seg_stage["failure_code"],
                      "failure_detail": seg_stage["failure_detail"]}
        return self._assemble_row(entry, chain, row_key, [ROLE_SEGMENT], raw_refs,
                                  resolution, None, None)

    def _failure_record(self, entry: str, chain: str, row_key: str,
                        seg_ref: dict | None, roles: list[str],
                        code: str, detail: str) -> dict:
        raw_refs = {ROLE_SEGMENT: seg_ref} if seg_ref is not None else {}
        resolution = {"ok": False, "kind": "hard", "accession": None, "entity_id": None,
                      "failure_code": code, "failure_detail": detail}
        return self._assemble_row(entry, chain, row_key, list(roles), raw_refs,
                                  resolution, None, None)

    def _assemble_row(self, entry: str, chain: str, row_key: str, roles: list[str],
                      raw_refs: dict, resolution: dict,
                      primary_out: dict | None, afdb_out: dict | None) -> dict:
        status = "ok"
        failure_code = None
        failure_detail = None
        exclusion_reason = None
        accessions: list[str] = []
        entity_id = resolution.get("entity_id")
        afdb_identity = None
        pdbe_mapping_ok = False
        afdb_ok = False

        res_kind = resolution["kind"]
        if res_kind in ("hard", "grammar"):
            status = "typed_failure"
            failure_code = resolution["failure_code"]
            failure_detail = resolution["failure_detail"]
        elif res_kind in ("unavailable", "unresolved", "ambiguous", "threshold"):
            exclusion_reason = resolution["failure_code"]
        elif res_kind == "ok":
            accessions = [resolution["accession"]]
            pdbe_mapping_ok = True
        else:  # pragma: no cover
            status = "typed_failure"
            failure_code = CensusFailureCode.INTERNAL_ERROR.value

        resolved_accession = resolution.get("accession")
        if pdbe_mapping_ok and primary_out is not None:
            if primary_out["kind"] in ("hard", "grammar"):
                status = "typed_failure"
                failure_code = primary_out["failure_code"]
                failure_detail = primary_out["failure_detail"]
                pdbe_mapping_ok = False
            elif primary_out["kind"] in ("unavailable", "absent"):
                exclusion_reason = primary_out["failure_code"]
                pdbe_mapping_ok = False
            elif resolved_accession is None or resolved_accession not in primary_out["accessions"]:
                exclusion_reason = CensusFailureCode.PRIMARY_ACCESSION_MISMATCH.value
                pdbe_mapping_ok = False

        if pdbe_mapping_ok and afdb_out is not None:
            if afdb_out["kind"] == "ok":
                afdb_identity = dict(afdb_out["identity"] or {})
                afdb_ok = True
            elif afdb_out["kind"] == "grammar":
                status = "typed_failure"
                failure_code = afdb_out["failure_code"]
                failure_detail = afdb_out["failure_detail"]
            elif afdb_out["kind"] == "hard":
                status = "typed_failure"
                failure_code = afdb_out["failure_code"]
                failure_detail = afdb_out["failure_detail"]
            else:  # absent / unavailable: availability outcome
                exclusion_reason = afdb_out["failure_code"]

        inclusion = "included" if (pdbe_mapping_ok and afdb_ok) else "excluded"
        return {
            "row_key": row_key,
            "entry": entry,
            "chain": chain,
            "source_roles_touched": sorted(set(roles)),
            "status": status,
            "raw_refs": {role: dict(ref) for role, ref in sorted(raw_refs.items())},
            "accessions": accessions,
            "entity_id": entity_id,
            "availability": {"pdbe_mapping": bool(pdbe_mapping_ok),
                             "afdb_prediction": bool(afdb_ok)},
            "inclusion_status": inclusion,
            "exclusion_reason": exclusion_reason,
            "failure_code": failure_code,
            "failure_detail": failure_detail,
            "afdb_identity": afdb_identity,
        }

    # -- entry point --------------------------------------------------------
    def run(self, rows: list[str], pins: dict, *, max_rows: int | None = None,
            schedule_pins: dict | None = None) -> dict:
        self._base_pins = dict(pins)
        self._schedule_pins = None if schedule_pins is None else dict(schedule_pins)
        selected = sorted(set(rows))
        if max_rows is not None:
            selected = selected[:max_rows]
        completed = self._load_progress()
        by_entry: dict[str, list[str]] = {}
        for row in selected:
            entry, chain = row.split("_", 1)
            by_entry.setdefault(entry, []).append(chain)
        pending = []
        for entry in sorted(by_entry):
            chains = sorted(by_entry[entry])
            if all(_row_key(entry, chain) in completed for chain in chains):
                continue
            pending.append((entry, chains))
        with concurrent.futures.ThreadPoolExecutor(max_workers=self._concurrency) as pool:
            futures = [pool.submit(self._process_entry, entry, chains)
                       for entry, chains in pending]
            for future in concurrent.futures.as_completed(futures):
                future.result()
        with self._records_lock:
            self._flush(final=True)
        summary = dict(self._manifest_obj(_utc_now())["summary"])
        resumed = sum(1 for row in selected
                      if _row_key(*row.split("_", 1)) in completed)
        summary.update({
            "flush_count": self._flush_count,
            "resumed_rows": resumed,
            "manifest_path": os.path.abspath(self._manifest_path),
            "raw_dir": os.path.abspath(self._raw_dir),
        })
        return summary


def _row_key(entry: str, chain: str) -> str:
    """Opaque candidate key: sha256 over the canonical census row literal."""
    return _sha256_hex(f"{entry}_{chain}".encode("ascii"))


# --------------------------------------------------------------------------
# CLI.
# --------------------------------------------------------------------------
def _schedule_pins(min_delay_s: float, jitter_s: float, backoff_initial_s: float,
                   max_attempts: int, timeout_s: float, concurrency: int) -> dict:
    return {
        "concurrency_max": int(concurrency),
        "inter_request_delay_ms_min": min_delay_s * 1000.0,
        "inter_request_delay_jitter_ms": jitter_s * 1000.0,
        "backoff_initial_s": backoff_initial_s,
        "backoff_factor": BACKOFF_FACTOR,
        "max_attempts": int(max_attempts),
        "request_timeout_s": timeout_s,
        "retry_after_cap_s": RETRY_AFTER_CAP_S,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="e420 gate-(e) phase (iii) availability census (predictor-blind).")
    parser.add_argument("--manifest", required=True, help="output manifest JSON path")
    parser.add_argument("--raw-dir", required=True, help="raw byte cache dir (hash filenames)")
    parser.add_argument("--resume", action="store_true",
                        help="skip rows completed in the progress file")
    parser.add_argument("--max-rows", type=int, default=None,
                        help="bound the census to the first N sorted rows")
    parser.add_argument("--project-root", default=ROOT)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--workers", type=int, default=CONCURRENCY_MAX)
    parser.add_argument("--min-delay-ms", type=float, default=INTER_REQUEST_DELAY_MS_MIN)
    parser.add_argument("--jitter-ms", type=float, default=INTER_REQUEST_JITTER_MS)
    parser.add_argument("--backoff-initial-s", type=float, default=BACKOFF_INITIAL_S)
    parser.add_argument("--max-attempts", type=int, default=MAX_ATTEMPTS)
    parser.add_argument("--registered-schedule-override", action="store_true",
                        help="required when deviating from the registered schedule (audit flag)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    min_delay_s = args.min_delay_ms / 1000.0
    jitter_s = args.jitter_ms / 1000.0
    registered = (min_delay_s == INTER_REQUEST_DELAY_MS_MIN / 1000.0
                  and jitter_s == INTER_REQUEST_JITTER_MS / 1000.0
                  and args.backoff_initial_s == BACKOFF_INITIAL_S
                  and args.max_attempts == MAX_ATTEMPTS
                  and args.workers <= CONCURRENCY_MAX)
    if not registered and not args.registered_schedule_override:
        print("refusing: schedule deviates from the registered bound without "
              "--registered-schedule-override", file=sys.stderr)
        return 2
    try:
        rows, pins = load_census_rows(args.project_root)
        runner = CensusRunner(
            manifest_path=args.manifest, raw_dir=args.raw_dir,
            project_root=args.project_root, resume=args.resume,
            concurrency=args.workers, min_delay_s=min_delay_s, jitter_s=jitter_s,
            backoff_initial_s=args.backoff_initial_s, max_attempts=args.max_attempts,
            seed=args.seed)
        schedule = _schedule_pins(min_delay_s, jitter_s, args.backoff_initial_s,
                                  args.max_attempts, REQUEST_TIMEOUT_S, args.workers)
        summary = runner.run(rows, pins, max_rows=args.max_rows, schedule_pins=schedule)
    except CensusError as exc:
        print(json.dumps({"typed_error": exc.code, "detail": exc.detail},
                         sort_keys=True, separators=(",", ":")), file=sys.stderr)
        return 3
    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
