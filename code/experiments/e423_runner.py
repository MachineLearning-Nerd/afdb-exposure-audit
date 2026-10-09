"""e423 dispatch-machinery repair: the e422 prospective batch runner.

ALLOCATED 2026-09-10 (w1:pG under delegated authority, DECISIONS same
date): a live dev-only shakedown of the frozen e422 runner against a
NON-registered historical window exposed two dispatch-machinery defects
— the e422 runner is content-SHA-bound at the gate-(a) freeze, so per
the freeze discipline this repair is a NEW E-number.  The e422 FROZEN
registration text is UNCHANGED and remains the registration of record;
e423 changes ONLY the dispatch implementation:

  1. the census request_options carried an INVALID "scoring_strategy":
     "none" (RCSB v2 enum rejects it, HTTP 400 on every page);
  2. a non-200 census response was marked ok by the scheduler
     fall-through and parsed as an EMPTY result set — a SILENT DROP
     (the shakedown produced a well-formed n_rows=0 checkpoint from an
     HTTP 400 error body).

Repairs here: the invalid request_options field is dropped (live
verified: HTTP 200); non-2xx responses are typed failures
(E420_HTTP_UNEXPECTED_STATUS in the shared scheduler); the census
response body is strictly validated (dict, integer total_count >= 0,
result_set list of identifier dicts) with the new typed code
E423_CENSUS_RESPONSE_MALFORMED — an empty batch window (total_count 0)
remains a legitimate outcome only on a validated 200.  Artifact paths,
schemas, allowlist lane roles, and the e422 study identity are
UNCHANGED (the e422 evaldriver consumes these checkpoints as-is); the
checkpoint records runner provenance "e423".

Executes one registered weekly batch end-to-end per the e422 FROZEN
registration (git-blob 6e1b9ce128ac3a1be4678463bdc777c25b50ca8d,
content sha256 dc7a3627…, freeze commit 9dd2380b) and e421-P §§2-3/§7:

  census receipt -> per-batch FREEZE receipt (BEFORE any label) -> label
  construction -> batch checkpoint (canonical-JSON SHA-256) -> [OTS
  anchor is a separate registered step] -> ledger append.

Lanes (docs/e421_source_allowlist.json, carried by e422 §2):
  - census enumeration: RCSB search v2 POST (metadata-only; the one
    allowlisted POST) with the registered predicate (X-RAY DIFFRACTION,
    resolution <= 2.5 A, initial release in the batch window);
  - mappings: GET https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry}
    (whole-entry UniProt mapping; ALL protein-chain triples processed —
    the exploratory yield basis);
  - PDBe entry files: GET https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif
    (reference structure);
  - AFDB metadata: GET https://alphafold.ebi.ac.uk/api/prediction/{acc};
  - AFDB child: GET exactly the manifest record's v6_cif_url.

Chain rule (registered here; the registration text is silent — flagged
for the gate-(c) disposition): every (entry, chain, accession) triple
from the mappings lane contributes rows when >= MIN_ALIGNED (30)
residues align; a UniProt position appears at most once per triple.
Label construction mirrors e421_pipeline.process_entry exactly
(Kabsch AF->ref, lDDT of sup vs ref, B-factor z over the whole mapped
chain) so the four channels {plddt, lddt, d_kabsch, b_factor_z} are the
registered e421-P §4 payload.

Schedule: e420 RequestScheduler gating (<=2 workers, >=250 ms + 100 ms
jitter inter-start, 2 s doubling backoff, max 3 attempts); the RCSB POST
is gated through the same lock via E423Scheduler.post.

Every failure is typed and tallied in the batch checkpoint — no silent
drops.  The presence manifest (experiments/e422_manifest.py) is loaded
from (or created at) the registered path and checkpointed per batch.
"""

from __future__ import annotations

import datetime
import json
import os
import sys
import time
from enum import Enum
from typing import Any, Callable, Mapping

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e420_census as census  # noqa: E402
import e421_pipeline as pipe  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402


# --------------------------------------------------------------------------- #
# Registered constants.
# --------------------------------------------------------------------------- #
RUNNER_SCHEMA = "e422-batch-checkpoint-v1"
FREEZE_SCHEMA = "e422-batch-freeze-v1"
MIN_ALIGNED = 30
CENSUS_PAGE_ROWS = 500
CENSUS_MAX_PAGES = 25
RESOLUTION_MAX_A = 2.5
EXPERIMENTAL_METHOD = "X-RAY DIFFRACTION"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
MAPPINGS_URL = "https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{entry}"
ENTRY_FILES_URL = "https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"
METADATA_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


class RunnerRole(str, Enum):
    """e422 source roles (allowlist lanes, e422-prefixed)."""

    CENSUS = "E422_CENSUS_SEARCH_V1"
    MAPPINGS = "E422_PDBE_MAPPINGS_V2"
    ENTRY_FILES = "E422_PDBE_ENTRY_FILES_V1"
    AFDB_METADATA = "E422_AFDB_METADATA_V1"
    AFDB_CHILD = "E422_AFDB_CHILD_V1"


class RunnerFailureCode(str, Enum):
    """Runner-level typed failures (manifest/source codes come from
    e422_manifest / census and pass through verbatim)."""

    CENSUS_QUERY_FAILED = "E422_CENSUS_QUERY_FAILED"
    CENSUS_RESPONSE_MALFORMED = "E423_CENSUS_RESPONSE_MALFORMED"
    CENSUS_RESULT_CAP = "E422_CENSUS_RESULT_CAP"
    WINDOW_NOT_CLOSED = "E423_WINDOW_NOT_CLOSED"
    CHECKPOINT_EXISTS = "E423_CHECKPOINT_EXISTS"
    T0_MISMATCH = "E423_T0_MISMATCH"
    CENSUS_DUPLICATE_ENTRIES = "E423_CENSUS_DUPLICATE_ENTRIES"
    MAPPING_UNAVAILABLE = "E422_MAPPING_UNAVAILABLE"
    MAPPING_PARSE_FAILED = "E422_MAPPING_PARSE_FAILED"
    NO_PROTEIN_MAPPING = "E422_NO_PROTEIN_MAPPING"
    METADATA_FETCH_FAILED = "E422_METADATA_FETCH_FAILED"
    CHILD_FETCH_FAILED = "E422_CHILD_FETCH_FAILED"
    REF_CIF_FETCH_FAILED = "E422_REF_CIF_FETCH_FAILED"
    AF_INSUFFICIENT_CA = "E422_AF_INSUFFICIENT_CA"
    REF_INSUFFICIENT_CA = "E422_REF_INSUFFICIENT_CA"
    INSUFFICIENT_ALIGNED = "E422_INSUFFICIENT_ALIGNED"
    KABSCH_FAILED = "E422_KABSCH_FAILED"


# --------------------------------------------------------------------------- #
# Batch window + census predicate (pure functions).
# --------------------------------------------------------------------------- #
# The frozen gate-(a) T0 bind (docs/e422_registration.md, content sha
# dc7a3627).  Must equal e423_dispatch.DEFAULT_T0 (pinned by test).
REGISTERED_T0 = "2026-09-11T00:00:00Z"


def batch_window_utc(t0_utc: str, k: int) -> tuple[str, str]:
    """Batch k window (start, close] as UTC ISO strings.

    ACCELERATION AMENDMENT (DECISIONS 2026-09-10, user-approved; the
    frozen §6 text is never edited — this amends the REGISTERED cadence
    only): registered batch 1 is a 48-HOUR window (T0, T0+2d]; batches
    2-12 remain 7-day and shift by batch 1's early close, so for the
    registered T0, close(k) = T0 + 2 + 7(k-1) days and the schedule
    stays gap-free: start(k) = close(k-1).  Nominal day-28/56/84
    evaluation labels therefore land at actual days 23/51/79
    (driver roles are keyed to batch INDEX, not dates — disclosed).
    Any NON-registered t0 (dev shakedowns, historical rehearsals) keeps
    the pure 7-day arithmetic: batch k = (T0+7(k-1)d, T0+7kd].
    """
    if k < 1:
        raise census.CensusError(RunnerFailureCode.CENSUS_QUERY_FAILED.value,
                                 f"batch index {k} < 1")
    t0 = manifest_mod.parse_t0(t0_utc)
    if t0_utc == REGISTERED_T0:
        start = (t0 if k == 1
                 else t0 + datetime.timedelta(days=2 + 7 * (k - 2)))
        close = t0 + datetime.timedelta(days=2 + 7 * (k - 1))
    else:
        start = t0 + datetime.timedelta(days=7 * (k - 1))
        close = t0 + datetime.timedelta(days=7 * k)
    fmt = lambda d: d.strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: E731
    return fmt(start), fmt(close)


def build_census_query(start_iso: str, close_iso: str, start_row: int = 0) -> dict:
    """The registered RCSB search v2 predicate (metadata-only)."""
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": EXPERIMENTAL_METHOD,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": RESOLUTION_MAX_A,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "greater",
                        "value": start_iso,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "less_or_equal",
                        "value": close_iso,
                    },
                },
            ],
        },
        "return_type": "entry",
        # e423 repair: NO scoring_strategy key.  The e422 runner sent
        # "scoring_strategy": "none", which the RCSB v2 enum rejects
        # (HTTP 400 on every page — found by the live shakedown).  The
        # field is scoring-only and irrelevant to this metadata predicate;
        # omitting it is live-verified (HTTP 200, 2026-09-10).
        "request_options": {
            "paginate": {"start": start_row, "rows": CENSUS_PAGE_ROWS},
            "results_content_type": ["experimental"],
        },
    }


# --------------------------------------------------------------------------- #
# Scheduler: RequestScheduler + gated POST (the one allowlisted POST).
# --------------------------------------------------------------------------- #
class E423Scheduler(census.RequestScheduler):
    """Adds a start-gated POST with identical receipt/backoff discipline.

    e423 repair: a non-2xx status that is not in the 429/5xx retry class
    is a TYPED failure (E420_HTTP_UNEXPECTED_STATUS) — the e422
    fall-through marked any other status (400/404/...) ok and its error
    body flowed into response parsing as if it were data (the silent
    empty-checkpoint defect the shakedown caught)."""

    @staticmethod
    def _settle_non_retryable(result: census.RequestResult,
                              status: int, content: bytes) -> None:
        result.failure_code = census.CensusFailureCode.HTTP_UNEXPECTED_STATUS.value
        result.failure_detail = census._clip(
            f"status={status} (non-retryable); body head: "
            + content[:200].decode("utf-8", errors="replace"))

    def post(self, role: RunnerRole, url: str,
             body: Mapping[str, Any]) -> census.RequestResult:
        payload = json.dumps(dict(body), sort_keys=True).encode("utf-8")
        result = census.RequestResult(role=role.value, url=url)
        for attempt in range(1, self._max_attempts + 1):
            self._enter_gate()
            try:
                resp = self._transport.post_json(url, payload, self._timeout_s)
            except Exception as exc:  # noqa: BLE001 - typed, never silent
                receipt = {"url": url, "method": "POST", "status": None,
                           "headers": {}, "bytes_len": 0, "sha256": None,
                           "utc": census._utc_now(), "role": role.value,
                           "attempt": attempt, "raw_name": None,
                           "retry_after_s": None, "next_backoff_s": None,
                           "error": type(exc).__name__}
                self._add_receipt(receipt)
                result.receipts.append(receipt)
                result.failure_code = census.CensusFailureCode.TRANSPORT_ERROR.value
                result.failure_detail = census._clip(f"{type(exc).__name__}: {exc}")
                return result
            content = resp.content
            sha = census._sha256_hex(content)
            raw_name = self._save_raw(content)
            receipt = {"url": url, "method": "POST", "status": int(resp.status),
                       "headers": dict(resp.headers), "bytes_len": len(content),
                       "sha256": sha, "utc": census._utc_now(), "role": role.value,
                       "attempt": attempt, "raw_name": raw_name,
                       "body_sha256": census._sha256_hex(payload),
                       "retry_after_s": None, "next_backoff_s": None, "error": None}
            self._add_receipt(receipt)
            result.receipts.append(receipt)
            result.status = int(resp.status)
            result.sha256 = sha
            result.bytes_len = len(content)
            result.content = content
            result.raw_name = raw_name
            if resp.status == 429 or 500 <= resp.status < 600:
                retry_after = census._retry_after_seconds(resp.headers)
                if retry_after is not None and retry_after > census.RETRY_AFTER_CAP_S:
                    retry_after = census.RETRY_AFTER_CAP_S
                receipt["retry_after_s"] = retry_after
                if attempt >= self._max_attempts:
                    result.failure_code = (
                        census.CensusFailureCode.HTTP_429_EXHAUSTED.value
                        if resp.status == 429
                        else census.CensusFailureCode.HTTP_5XX_EXHAUSTED.value)
                    result.failure_detail = census._clip(
                        f"status={resp.status} after {attempt} attempts")
                    return result
                delay = self._backoff_initial_s * (self._backoff_factor ** (attempt - 1))
                if retry_after is not None:
                    delay = max(delay, retry_after)
                receipt["next_backoff_s"] = delay
                self._sleeper(delay)  # backoff outside the start gate
                continue
            if not (200 <= resp.status < 300):
                self._settle_non_retryable(result, int(resp.status), content)
                return result
            result.ok = True
            return result
        return result  # pragma: no cover - loop always returns


# --------------------------------------------------------------------------- #
# Mapping parsing.
# --------------------------------------------------------------------------- #
def parse_mappings(payload: bytes, entry: str) -> list[dict]:
    """Flatten the PDBe whole-entry UniProt mapping to chain triples.

    Returns [{chain_id, entity_id, accession, res_map}] where res_map maps
    author residue number -> UniProt position (the e421 offset rule).
    Raises MAPPING_PARSE_FAILED on an unparseable payload.

    e423 repair (live shakedown 2026-09-10): mapping items are FLAT on the
    live PDBe v2 payload — ``chain_id`` at top level (with struct_asym_id
    alongside), exactly as the field-proven e420 resolve_chain /
    e421_pipeline parsers read them.  The e422 parser assumed a nested
    ``{"chain": {"chain_id": ...}}`` that only ever existed in its own
    test fixture (KeyError 'chain' on the first live receipt).  Every
    structural surprise is the TYPED MAPPING_PARSE_FAILED — never an
    untyped crash mid-batch.

    e423 hardening (detached review MAJOR-2, 2026-09-10): the typed
    envelope now covers the WHOLE walk — a list/null UniProt block, a
    non-dict accession value, a non-list mappings list, null/non-numeric
    residue bounds, and the res_map arithmetic itself are all inside it
    (they previously sat on untyped paths outside the per-item try).
    """
    try:
        doc = json.loads(payload.decode("utf-8"))
        uniprot = doc[entry.lower()]["UniProt"]
    except Exception as exc:  # noqa: BLE001 - typed
        raise census.CensusError(
            RunnerFailureCode.MAPPING_PARSE_FAILED.value,
            f"{entry}: {type(exc).__name__}",
        ) from exc
    if not isinstance(uniprot, dict):
        raise census.CensusError(
            RunnerFailureCode.MAPPING_PARSE_FAILED.value,
            f"{entry}: UniProt block is {type(uniprot).__name__}, not an object")
    triples: list[dict] = []
    for accession in sorted(uniprot):
        try:
            acc_block = uniprot[accession]
            mappings = (acc_block.get("mappings", [])
                        if isinstance(acc_block, dict) else None)
            if not isinstance(mappings, list):
                raise TypeError(
                    "accession value is not an object with a mappings list")
            for mapping in mappings:
                chain_id = mapping["chain_id"]
                entity_id = mapping["entity_id"]
                unp_s = mapping["unp_start"]
                unp_e = mapping["unp_end"]
                s = (mapping["start"]["author_residue_number"]
                     or mapping["start"]["residue_number"])
                e = (mapping["end"]["author_residue_number"]
                     or mapping["end"]["residue_number"])
                if any(not isinstance(v, int) or isinstance(v, bool)
                       for v in (unp_s, unp_e, s, e)):
                    raise TypeError("non-integer residue bounds")
                res_map = {s + off: unp_s + off
                           for off in range(max(0, e - s + 1))
                           if unp_s + off <= unp_e}
                if res_map:
                    triples.append({"chain_id": chain_id, "entity_id": entity_id,
                                    "accession": accession, "res_map": res_map})
        except (KeyError, TypeError) as exc:  # noqa: BLE001 - typed
            raise census.CensusError(
                RunnerFailureCode.MAPPING_PARSE_FAILED.value,
                f"{entry}/{accession}: malformed mapping item "
                f"({type(exc).__name__}: {exc})") from exc
    return triples


# --------------------------------------------------------------------------- #
# Label construction for one triple (mirrors e421_pipeline.process_entry).
# --------------------------------------------------------------------------- #
def label_triple(triple: dict, entry: str, af_cif_text: str, ref_cif_text: str,
                 af_sha: str, ref_sha: str) -> dict:
    """Rows + typed status for one (entry, chain, accession) triple."""
    out = {"entry": entry, "chain": triple["chain_id"],
           "accession": triple["accession"], "rows": [], "status": None,
           "af_sha256": af_sha, "ref_sha256": ref_sha}
    af_atoms = pipe.parse_cif_ca(af_cif_text)
    if len(af_atoms) < 3:
        out["status"] = RunnerFailureCode.AF_INSUFFICIENT_CA.value
        return out
    af_by_pos = {a["residue_number"]: a for a in af_atoms}
    ref_atoms = pipe.parse_cif_ca(ref_cif_text, chain_filter=triple["chain_id"])
    if len(ref_atoms) < MIN_ALIGNED:
        out["status"] = RunnerFailureCode.REF_INSUFFICIENT_CA.value
        return out
    ref_ca, af_ca, unp_positions, ref_idx = [], [], [], []
    seen_pos: set[int] = set()
    for i, atom in enumerate(ref_atoms):
        unp = triple["res_map"].get(atom["residue_number"])
        if unp is None or unp in seen_pos:
            continue
        af_atom = af_by_pos.get(unp)
        if af_atom is None:
            continue
        seen_pos.add(unp)
        ref_ca.append([atom["x"], atom["y"], atom["z"]])
        af_ca.append([af_atom["x"], af_atom["y"], af_atom["z"]])
        unp_positions.append(unp)
        ref_idx.append(i)
    if len(ref_ca) < MIN_ALIGNED:
        out["status"] = RunnerFailureCode.INSUFFICIENT_ALIGNED.value
        return out
    sup, dists = pipe.kabsch_superpose(af_ca, ref_ca)
    if sup is None:
        out["status"] = RunnerFailureCode.KABSCH_FAILED.value
        return out
    lddt, _nbr = pipe.lddt_calc(sup, np.array(ref_ca))
    # F channel: z-score over ALL mapped-chain CA b-factors (the e421 rule),
    # indexed back to the aligned subset via ref_idx.
    z, _med, _iqr = pipe.bfactor_z([a["b_factor"] for a in ref_atoms])
    af_bf = {a["residue_number"]: a["b_factor"] for a in af_atoms}
    rows = []
    for j, unp in enumerate(unp_positions):
        plddt = af_bf.get(unp)
        l = None if np.isnan(lddt[j]) else float(lddt[j])
        rows.append({
            "entry": entry, "chain": triple["chain_id"],
            "accession": triple["accession"], "uniprot_pos": int(unp),
            "plddt": round(float(plddt), 4) if plddt is not None else None,
            "lddt": round(l, 4) if l is not None else None,
            "d_kabsch": round(float(dists[j]), 4),
            "b_factor_z": round(float(z[ref_idx[j]]), 4),
        })
    out["rows"] = rows
    out["status"] = "ok"
    return out


# --------------------------------------------------------------------------- #
# Batch execution.
# --------------------------------------------------------------------------- #
def run_batch(
    t0_utc: str,
    batch_index: int,
    *,
    transport: Any,
    raw_dir: str,
    out_dir: str,
    manifest_path: str | None = None,
    max_entries: int | None = None,
    sleeper: Callable[[float], None] = time.sleep,
    now_iso: str | None = None,
    parse_mappings_fn: Callable[[bytes, str], list] | None = None,
    label_triple_fn: Callable[..., dict] | None = None,
    labeller_id: str | None = None,
) -> dict:
    """Execute batch k; returns the batch-checkpoint artifact (and writes
    the freeze receipt + checkpoint files under out_dir).

    e427 (docs/e427_registration.md): the mapping parser and triple labeller
    are injectable.  The defaults are this module's original functions (used
    by the e427 replay self-check); the registered dispatch passes the e427
    label-space repair and records ``labeller_id`` in the checkpoint."""
    parse_mappings_fn = parse_mappings_fn or parse_mappings
    label_triple_fn = label_triple_fn or label_triple
    start_iso, close_iso = batch_window_utc(t0_utc, batch_index)
    # Registered discipline, mechanically enforced (DECISIONS 2026-09-10,
    # pre-T0 hardening): a dispatch over an unclosed window is NOT the
    # registered population — the query would silently admit a partial
    # week.  Refuse before any filesystem or network side effect.
    # now == close is allowed: the (start, close] interval is complete at
    # that instant.
    now_dt = (manifest_mod.parse_t0(now_iso) if now_iso is not None
              else datetime.datetime.now(datetime.timezone.utc))
    if now_dt < manifest_mod.parse_t0(close_iso):
        raise census.CensusError(
            RunnerFailureCode.WINDOW_NOT_CLOSED.value,
            f"batch {batch_index}: window ({start_iso}, {close_iso}] is not "
            f"closed (now {now_dt.isoformat()}); refusing pre-close dispatch "
            "— the registered population is only defined over a closed "
            "window (runbook)")
    # Double-dispatch safety, mechanically enforced (DECISIONS
    # 2026-09-10): a checkpoint is an anchored dispatch artifact.  A
    # second run over the same batch would silently REPLACE it with a
    # rebuilt freeze/label construction (internally consistent, so the
    # eval-side digest chain would still load) — exactly the
    # indistinguishable-from-legitimate class the window guard exists
    # for.  Re-dispatch is a deliberate action: dated DECISIONS entry +
    # manual removal of the prior checkpoint, never an overwrite.
    checkpoint_path = os.path.join(
        out_dir, f"e422_batch{batch_index:02d}_checkpoint.json")
    if os.path.exists(checkpoint_path):
        raise census.CensusError(
            RunnerFailureCode.CHECKPOINT_EXISTS.value,
            f"batch {batch_index}: checkpoint already exists at "
            f"{checkpoint_path}; refusing to overwrite an anchored "
            "dispatch artifact — re-dispatch requires a dated DECISIONS "
            "entry and manual removal of the prior checkpoint (runbook: "
            "one dispatch per batch)")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    # Manifest integrity BEFORE any network phase (detached review
    # 2026-09-10, wf_d8bff144): an existing manifest is digest-verified
    # and T0-bound here, so a tampered manifest or a T0 that disagrees
    # with the accumulated chain fails before the census fetches rather
    # than after them.  The instance is reused for admission below.
    if manifest_path and os.path.exists(manifest_path):
        manifest = manifest_mod.ManifestBuilder.load(manifest_path)
        if manifest.t0_utc != t0_utc:
            raise census.CensusError(
                RunnerFailureCode.T0_MISMATCH.value,
                f"batch {batch_index}: manifest t0_utc "
                f"{manifest.t0_utc!r} != run t0_utc {t0_utc!r} — refusing "
                "to accumulate a chain across two different registrations "
                "(detached review 2026-09-10)")
    else:
        manifest = None
    scheduler = E423Scheduler(transport, raw_dir, sleeper=sleeper)

    # ---- 1. census enumeration (receipted POST, paginated) ----------------
    entries: list[str] = []
    page = 0
    while True:
        page += 1
        if page > CENSUS_MAX_PAGES:
            raise census.CensusError(
                RunnerFailureCode.CENSUS_RESULT_CAP.value,
                f"batch {batch_index}: exceeded {CENSUS_MAX_PAGES} x {CENSUS_PAGE_ROWS} rows",
            )
        query = build_census_query(start_iso, close_iso,
                                   start_row=(page - 1) * CENSUS_PAGE_ROWS)
        result = scheduler.post(RunnerRole.CENSUS, SEARCH_URL, query)
        if not result.ok or result.status != 200:
            raise census.CensusError(
                RunnerFailureCode.CENSUS_QUERY_FAILED.value,
                f"batch {batch_index} page {page}: "
                f"{result.failure_code} status={result.status}",
            )
        # e423 repair: strict response validation.  A malformed body is a
        # TYPED failure — the e422 runner defaulted missing total_count to
        # 0 and treated an error document as an empty (legitimate) window.
        try:
            doc = json.loads(result.content.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise census.CensusError(
                RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
                f"batch {batch_index} page {page}: undecodable census body "
                f"({type(exc).__name__})") from exc
        if not isinstance(doc, dict) or "total_count" not in doc \
                or not isinstance(doc.get("result_set"), list):
            raise census.CensusError(
                RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
                f"batch {batch_index} page {page}: census body is not a "
                "well-formed result document")
        raw_total = doc["total_count"]
        # e423 hardening (detached review MINOR-4, 2026-09-10): NO coercion —
        # int("7"), int(7.9), and int(True) all passed the old int() gate;
        # total_count must be a JSON integer, full stop.
        if isinstance(raw_total, bool) or not isinstance(raw_total, int):
            raise census.CensusError(
                RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
                f"batch {batch_index} page {page}: non-integer total_count "
                f"({type(raw_total).__name__})")
        total = raw_total
        if total < 0:
            raise census.CensusError(
                RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
                f"batch {batch_index} page {page}: negative total_count")
        ids = []
        for d in doc["result_set"]:
            if not isinstance(d, dict) or not isinstance(d.get("identifier"), str):
                raise census.CensusError(
                    RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
                    f"batch {batch_index} page {page}: result_set item has "
                    "no string identifier")
            ids.append(d["identifier"].lower())
        entries.extend(ids)
        if len(entries) >= total or not ids:
            break
    if len(entries) < total:
        # e423 hardening (review MINOR-4): a page stream that ends short of
        # total_count is a TRUNCATED census — a typed failure, never a
        # silently reduced registered population.
        raise census.CensusError(
            RunnerFailureCode.CENSUS_RESPONSE_MALFORMED.value,
            f"batch {batch_index}: census ended at {len(entries)} of "
            f"total_count {total} (truncated pagination)")
    if len(set(entries)) != len(entries):
        # e423 hardening (detached review wf_d8bff144, 2026-09-10): a
        # paginated census that repeats an identifier is a malformed
        # page stream — silent duplicates would inflate the census
        # population count (and the freeze receipt's entries tally)
        # while every digest still verifies.
        raise census.CensusError(
            RunnerFailureCode.CENSUS_DUPLICATE_ENTRIES.value,
            f"batch {batch_index}: census returned duplicate "
            f"identifiers ({len(entries)} entries, "
            f"{len(set(entries))} unique) — refusing a malformed page "
            "stream rather than silently deduplicating the registered "
            "population")
    if max_entries is not None:
        entries = entries[:max_entries]

    # ---- 2. mappings per entry (receipted GET) ----------------------------
    triples_by_entry: dict[str, list[dict]] = {}
    mapping_failures: list[dict] = []
    for entry in entries:
        result = scheduler.get(RunnerRole.MAPPINGS, MAPPINGS_URL.format(entry=entry))
        if not result.ok:
            mapping_failures.append({"entry": entry,
                                     "code": RunnerFailureCode.MAPPING_UNAVAILABLE.value,
                                     "detail": result.failure_code})
            continue
        try:
            triples = parse_mappings_fn(result.content, entry)
        except census.CensusError as exc:
            mapping_failures.append({"entry": entry, "code": exc.code,
                                     "detail": exc.detail})
            continue
        if not triples:
            mapping_failures.append({"entry": entry,
                                     "code": RunnerFailureCode.NO_PROTEIN_MAPPING.value,
                                     "detail": "no UniProt mappings"})
            continue
        triples_by_entry[entry] = triples

    # ---- 3. manifest admission per accession (receipted metadata GET) -----
    if manifest is None:  # early pre-network verification found no file
        manifest = manifest_mod.ManifestBuilder(t0_utc)
    admitted_new: list[str] = []
    manifest_failures: list[dict] = []
    needed = sorted({t["accession"] for ts in triples_by_entry.values() for t in ts})
    for accession in needed:
        if accession in manifest.records:
            continue
        result = scheduler.get(RunnerRole.AFDB_METADATA,
                               METADATA_URL.format(acc=accession))
        if result.status == 404:
            manifest_failures.append(
                {"accession": accession,
                 "code": manifest_mod.ManifestFailureCode.ABSENT_FROM_SNAPSHOT.value,
                 "detail": "metadata 404"})
            continue
        if not result.ok or result.status != 200:
            manifest_failures.append(
                {"accession": accession,
                 "code": RunnerFailureCode.METADATA_FETCH_FAILED.value,
                 "detail": result.failure_code or f"status={result.status}"})
            continue
        try:
            models = json.loads(result.content.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            manifest_failures.append(
                {"accession": accession,
                 "code": manifest_mod.ManifestFailureCode.MANIFEST_SCHEMA_INVALID.value,
                 "detail": f"undecodable metadata body ({type(exc).__name__})"})
            continue
        if isinstance(models, dict):
            models = [models]
        # e423 hardening (detached review MINOR-3, 2026-09-10): the FROZEN
        # manifest module assumes a non-empty array of JSON objects —
        # pre-validate here so a scalar/null/item-typed payload reaches it
        # as a TYPED schema failure, never an untyped crash.
        if not (isinstance(models, list) and models
                and all(isinstance(m, dict) for m in models)):
            manifest_failures.append(
                {"accession": accession,
                 "code": manifest_mod.ManifestFailureCode.MANIFEST_SCHEMA_INVALID.value,
                 "detail": "metadata payload is not a JSON object or "
                           "non-empty array of objects"})
            continue
        try:
            manifest.admit(models)
            admitted_new.append(accession)
        except manifest_mod.ManifestError as exc:
            manifest_failures.append({"accession": accession, "code": exc.code,
                                      "detail": exc.detail})
        except (TypeError, AttributeError, KeyError, IndexError,
                ValueError) as exc:
            # frozen-module validation gaps (e.g. int(None) latestVersion)
            # surface as TYPED schema failures, never untyped crashes
            manifest_failures.append(
                {"accession": accession,
                 "code": manifest_mod.ManifestFailureCode.MANIFEST_SCHEMA_INVALID.value,
                 "detail": f"{type(exc).__name__}: {exc}"})

    # ---- 4. FREEZE receipt (BEFORE any label construction) ---------------
    freeze_utc = now_iso or census._utc_now()
    freeze_receipt = {
        "schema": FREEZE_SCHEMA,
        "runner": "e423 (dispatch-machinery repair; DECISIONS 2026-09-10)",
        "e_number": "e422",
        "t0_utc": t0_utc,
        "batch_index": batch_index,
        "window": {"start": start_iso, "close": close_iso},
        "frozen_utc": freeze_utc,
        "entries": len(triples_by_entry),
        "manifest_accessions": len(manifest.records),
        "manifest_checkpoint_sha256": manifest.checkpoint(),
        "request_count": len(scheduler.receipts),
        "note": "freeze precedes label construction (e421-P section 7)",
    }
    freeze_sha = manifest_mod.canonical_sha256(freeze_receipt)
    freeze_receipt["sha256"] = freeze_sha
    census._atomic_write_json(
        os.path.join(out_dir, f"e422_batch{batch_index:02d}_freeze.json"),
        freeze_receipt)

    # ---- 5. label construction per entry then triple (children receipted) -
    rows: list[dict] = []
    triple_failures: list[dict] = []
    af_text_by_acc: dict[str, str] = {}
    af_sha_by_acc: dict[str, str] = {}
    ref_sha_by_entry: dict[str, str] = {}
    verified: set[str] = set()

    def fail(entry: str, t: dict, code: str, detail: str) -> None:
        triple_failures.append({"entry": entry, "chain": t["chain_id"],
                                "accession": t["accession"], "code": code,
                                "detail": detail})

    for entry in sorted(triples_by_entry):
        ref_result = scheduler.get(RunnerRole.ENTRY_FILES,
                                   ENTRY_FILES_URL.format(entry=entry))
        if not ref_result.ok:
            for t in triples_by_entry[entry]:
                fail(entry, t, RunnerFailureCode.REF_CIF_FETCH_FAILED.value,
                     f"status={ref_result.status}")
            continue
        ref_sha_by_entry[entry] = ref_result.sha256
        try:
            ref_text = ref_result.content.decode("utf-8")
        except UnicodeDecodeError:
            ref_text = ref_result.content.decode("utf-8", errors="replace")
        for t in sorted(triples_by_entry[entry],
                        key=lambda x: (x["chain_id"], x["accession"])):
            acc = t["accession"]
            if acc not in manifest.records:
                fail(entry, t,
                     manifest_mod.ManifestFailureCode.ABSENT_FROM_SNAPSHOT.value,
                     "not in manifest")
                continue
            if acc not in af_text_by_acc:
                # registered: verification metadata fetch + child fetch at
                # the accession's first labeling use (batch-time §2.2 check).
                ver = scheduler.get(RunnerRole.AFDB_METADATA,
                                    METADATA_URL.format(acc=acc))
                if ver.ok and ver.status == 200:
                    try:
                        model = json.loads(ver.content.decode("utf-8"))
                    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
                        fail(entry, t,
                             RunnerFailureCode.METADATA_FETCH_FAILED.value,
                             f"undecodable verification body "
                             f"({type(exc).__name__})")
                        continue
                    if isinstance(model, dict):
                        model = [model]
                    # e423 hardening (review MINOR-3): verify_against_manifest
                    # assumes a JSON object — a scalar/list/null verification
                    # payload is a TYPED lane failure, never an untyped crash.
                    if not (isinstance(model, list) and model
                            and isinstance(model[0], dict)):
                        fail(entry, t,
                             RunnerFailureCode.METADATA_FETCH_FAILED.value,
                             "verification payload is not a JSON object")
                        continue
                    try:
                        manifest_mod.verify_against_manifest(manifest.records[acc],
                                                             model[0])
                        verified.add(acc)
                    except manifest_mod.ManifestError as exc:
                        fail(entry, t, exc.code, exc.detail)
                        continue
                    except (TypeError, AttributeError, KeyError, IndexError,
                            ValueError) as exc:
                        fail(entry, t,
                             RunnerFailureCode.METADATA_FETCH_FAILED.value,
                             type(exc).__name__)
                        continue
                else:
                    fail(entry, t, RunnerFailureCode.METADATA_FETCH_FAILED.value,
                         ver.failure_code)
                    continue
                child = scheduler.get(RunnerRole.AFDB_CHILD,
                                      manifest.records[acc]["v6_cif_url"])
                if not child.ok:
                    fail(entry, t, RunnerFailureCode.CHILD_FETCH_FAILED.value,
                         child.failure_code)
                    continue
                af_sha_by_acc[acc] = child.sha256
                try:
                    af_text_by_acc[acc] = child.content.decode("utf-8")
                except UnicodeDecodeError:
                    af_text_by_acc[acc] = child.content.decode("utf-8",
                                                               errors="replace")
            labeled = label_triple_fn(t, entry, af_text_by_acc[acc], ref_text,
                                   af_sha_by_acc[acc], ref_sha_by_entry[entry])
            if labeled["status"] == "ok":
                rows.extend(labeled["rows"])
            else:
                fail(entry, t, labeled["status"], "")

    # ---- 6. batch checkpoint (canonical-JSON SHA-256) ---------------------
    checkpoint = {
        "schema": RUNNER_SCHEMA,
        "runner": "e423 (dispatch-machinery repair; DECISIONS 2026-09-10)",
        "e_number": "e422",
        "t0_utc": t0_utc,
        "batch_index": batch_index,
        "window": {"start": start_iso, "close": close_iso},
        "completed_utc": now_iso or census._utc_now(),
        "census_entries": len(entries),
        # e423 (review F6): an explicit bounding-override field — never an
        # implicit bound inferable only from census_entries.
        "max_entries": max_entries,
        "entries_with_mappings": len(triples_by_entry),
        "manifest_checkpoint_sha256": manifest.checkpoint(),
        "manifest_admitted_new": sorted(admitted_new),
        "freeze_receipt_sha256": freeze_sha,
        "n_rows": len(rows),
        "rows": rows,
        "exclusions": {
            "mapping_failures": mapping_failures,
            "manifest_failures": manifest_failures,
            "triple_failures": triple_failures,
        },
        "verified_accessions": sorted(verified),
        "receipts": scheduler.receipts,
    }
    if labeller_id is not None:
        checkpoint["labeller"] = labeller_id
    checkpoint["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in checkpoint.items() if k != "sha256"})
    path = checkpoint_path
    census._atomic_write_json(path, checkpoint)
    if manifest_path:
        manifest.save(manifest_path)
    return checkpoint


__all__ = [
    "RUNNER_SCHEMA",
    "FREEZE_SCHEMA",
    "RunnerRole",
    "RunnerFailureCode",
    "E423Scheduler",
    "batch_window_utc",
    "build_census_query",
    "parse_mappings",
    "label_triple",
    "run_batch",
]
