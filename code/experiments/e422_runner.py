"""e422 receipted prospective batch runner (gate-(e) dispatch implementation).

Executes one registered weekly batch end-to-end per the amended e422
registration (blob 9cdf7056, commit 9dec8a20) and e421-P §§2-3/§7:

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
is gated through the same lock via E422Scheduler.post.

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
    CENSUS_RESULT_CAP = "E422_CENSUS_RESULT_CAP"
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
def batch_window_utc(t0_utc: str, k: int) -> tuple[str, str]:
    """Batch k window: (T0 + 7(k-1) days, T0 + 7k days] as UTC ISO strings."""
    if k < 1:
        raise census.CensusError(RunnerFailureCode.CENSUS_QUERY_FAILED.value,
                                 f"batch index {k} < 1")
    t0 = manifest_mod.parse_t0(t0_utc)
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
        "request_options": {
            "paginate": {"start": start_row, "rows": CENSUS_PAGE_ROWS},
            "results_content_type": ["experimental"],
            "scoring_strategy": "none",
        },
    }


# --------------------------------------------------------------------------- #
# Scheduler: RequestScheduler + gated POST (the one allowlisted POST).
# --------------------------------------------------------------------------- #
class E422Scheduler(census.RequestScheduler):
    """Adds a start-gated POST with identical receipt/backoff discipline."""

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
    """
    try:
        doc = json.loads(payload.decode("utf-8"))
        uniprot = doc[entry.lower()]["UniProt"]
    except Exception as exc:  # noqa: BLE001 - typed
        raise census.CensusError(
            RunnerFailureCode.MAPPING_PARSE_FAILED.value,
            f"{entry}: {type(exc).__name__}",
        ) from exc
    triples: list[dict] = []
    for accession in sorted(uniprot):
        for mapping in uniprot[accession].get("mappings", []):
            chain_id = mapping["chain"]["chain_id"]
            entity_id = mapping["entity_id"]
            unp_s = mapping["unp_start"]
            unp_e = mapping["unp_end"]
            s = (mapping["start"]["author_residue_number"]
                 or mapping["start"]["residue_number"])
            e = (mapping["end"]["author_residue_number"]
                 or mapping["end"]["residue_number"])
            res_map = {s + off: unp_s + off for off in range(max(0, e - s + 1))
                       if unp_s + off <= unp_e}
            if res_map:
                triples.append({"chain_id": chain_id, "entity_id": entity_id,
                                "accession": accession, "res_map": res_map})
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
) -> dict:
    """Execute batch k; returns the batch-checkpoint artifact (and writes
    the freeze receipt + checkpoint files under out_dir)."""
    start_iso, close_iso = batch_window_utc(t0_utc, batch_index)
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    scheduler = E422Scheduler(transport, raw_dir, sleeper=sleeper)

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
        if not result.ok:
            raise census.CensusError(
                RunnerFailureCode.CENSUS_QUERY_FAILED.value,
                f"batch {batch_index} page {page}: {result.failure_code}",
            )
        doc = json.loads(result.content.decode("utf-8"))
        total = int(doc.get("total_count", 0))
        ids = [str(d["identifier"]).lower() for d in doc.get("result_set", [])]
        entries.extend(ids)
        if len(entries) >= total or not ids:
            break
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
            triples = parse_mappings(result.content, entry)
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
    if manifest_path and os.path.exists(manifest_path):
        manifest = manifest_mod.ManifestBuilder.load(manifest_path)
    else:
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
            if isinstance(models, dict):
                models = [models]
            manifest.admit(models)
            admitted_new.append(accession)
        except manifest_mod.ManifestError as exc:
            manifest_failures.append({"accession": accession, "code": exc.code,
                                      "detail": exc.detail})
        except json.JSONDecodeError as exc:
            manifest_failures.append(
                {"accession": accession,
                 "code": manifest_mod.ManifestFailureCode.MANIFEST_SCHEMA_INVALID.value,
                 "detail": str(exc)})

    # ---- 4. FREEZE receipt (BEFORE any label construction) ---------------
    freeze_utc = now_iso or census._utc_now()
    freeze_receipt = {
        "schema": FREEZE_SCHEMA,
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
                        if isinstance(model, dict):
                            model = [model]
                        manifest_mod.verify_against_manifest(manifest.records[acc],
                                                             model[0])
                        verified.add(acc)
                    except manifest_mod.ManifestError as exc:
                        fail(entry, t, exc.code, exc.detail)
                        continue
                    except (json.JSONDecodeError, KeyError, IndexError) as exc:
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
            labeled = label_triple(t, entry, af_text_by_acc[acc], ref_text,
                                   af_sha_by_acc[acc], ref_sha_by_entry[entry])
            if labeled["status"] == "ok":
                rows.extend(labeled["rows"])
            else:
                fail(entry, t, labeled["status"], "")

    # ---- 6. batch checkpoint (canonical-JSON SHA-256) ---------------------
    checkpoint = {
        "schema": RUNNER_SCHEMA,
        "e_number": "e422",
        "t0_utc": t0_utc,
        "batch_index": batch_index,
        "window": {"start": start_iso, "close": close_iso},
        "completed_utc": now_iso or census._utc_now(),
        "census_entries": len(entries),
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
    checkpoint["sha256"] = manifest_mod.canonical_sha256(
        {k: v for k, v in checkpoint.items() if k != "sha256"})
    path = os.path.join(out_dir, f"e422_batch{batch_index:02d}_checkpoint.json")
    census._atomic_write_json(path, checkpoint)
    if manifest_path:
        manifest.save(manifest_path)
    return checkpoint


__all__ = [
    "RUNNER_SCHEMA",
    "FREEZE_SCHEMA",
    "RunnerRole",
    "RunnerFailureCode",
    "E422Scheduler",
    "batch_window_utc",
    "build_census_query",
    "parse_mappings",
    "label_triple",
    "run_batch",
]
