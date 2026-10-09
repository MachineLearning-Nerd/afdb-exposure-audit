"""Offline tests for the e422 batch runner (no network; injected transport)."""

from __future__ import annotations

import json
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e420_census as census  # noqa: E402
import e421_pipeline as pipe  # noqa: E402
import e422_manifest as manifest_mod  # noqa: E402
import e422_runner as runner  # noqa: E402

T0 = "2026-09-10T00:00:00Z"


# --------------------------------------------------------------------------- #
# Fixtures: mmCIF text + PDBe mapping payloads.
# --------------------------------------------------------------------------- #
def cif_text(residues, b_start, coords):
    """Build a minimal mmCIF atom_site loop.

    residues: list of (chain, auth_seq); coords: list of (x, y, z) aligned
    with residues; b_start: b-factor = b_start + index.
    """
    lines = ["data_test", "loop_",
             "_atom_site.group_PDB", "_atom_site.label_atom_id",
             "_atom_site.label_comp_id", "_atom_site.label_asym_id",
             "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
             "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
             "_atom_site.B_iso_or_equiv", "_atom_site.label_alt_id",
             "_atom_site.pdbx_PDB_model_num"]
    for i, (chain, seq) in enumerate(residues):
        x, y, z = coords[i]
        lines.append(f"ATOM CA ALA {chain} {chain} {seq} {x:.3f} {y:.3f} {z:.3f} "
                     f"{b_start + i:.2f} . 1")
    return "\n".join(lines) + "\n"


def helix(n, scale=3.0):
    return [(scale * math.cos(0.4 * i), scale * math.sin(0.4 * i), 1.5 * i)
            for i in range(n)]


def mapping_payload(entry, pairs):
    """pairs: list of (accession, chain_id, n_res)."""
    uniprot = {}
    for accession, chain, n in pairs:
        uniprot[accession] = {"mappings": [{
            "chain": {"chain_id": chain, "author_asym_id": chain},
            "entity_id": 1,
            "unp_start": 1, "unp_end": n,
            "start": {"author_residue_number": 1, "residue_number": 1},
            "end": {"author_residue_number": n, "residue_number": n},
        }]}
    return {entry: {"UniProt": uniprot}}


def af_metadata(accession, created="2025-01-01", version=6):
    return {
        "uniprotAccession": accession,
        "modelCreatedDate": created,
        "latestVersion": version,
        "modelEntityId": f"AF-{accession}-F1",
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif",
    }


class FakeTransport:
    """Offline transport with call recording (get + post_json)."""

    def __init__(self, routes=None, post_routes=None):
        self.routes = routes or {}
        self.post_routes = post_routes or {}
        self.calls = []

    def _respond(self, routes, url, key):
        self.calls.append(key)
        outcome = routes[url]
        if callable(outcome):
            outcome = outcome()
        status, content = outcome
        if isinstance(content, str):
            content = content.encode("utf-8")
        return census.HttpResponse(status, {"Content-Type": "application/json"},
                                   content)

    def get(self, url, timeout_s):
        return self._respond(self.routes, url, ("GET", url))

    def post_json(self, url, payload, timeout_s):
        return self._respond(self.post_routes, url,
                             ("POST", url, payload.decode("utf-8")))

    def count_role_url(self, fragment):
        return sum(1 for c in self.calls if fragment in c[1])


def census_body(entries):
    return json.dumps({"total_count": len(entries),
                       "result_set": [{"identifier": e} for e in entries]})


def standard_routes(entries=("aaaa",), pairs=(("P12345", "A", 34),),
                    n_res=34, af_created="2025-01-01"):
    entry = entries[0]
    accession = pairs[0][0]
    ref = cif_text([(pairs[0][1], i + 1) for i in range(n_res)], 20.0, helix(n_res))
    af = cif_text([("A", i + 1) for i in range(n_res)], 50.0, helix(n_res))
    metadata = json.dumps(af_metadata(accession, created=af_created))
    child_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.cif"
    routes = {
        runner.MAPPINGS_URL.format(entry=entry): (
            200, json.dumps(mapping_payload(entry, list(pairs)))),
        runner.METADATA_URL.format(acc=accession): (200, metadata),
        runner.ENTRY_FILES_URL.format(entry=entry): (200, ref),
        child_url: (200, af),
    }
    post_routes = {
        runner.SEARCH_URL: (200, census_body(list(entries))),
    }
    return FakeTransport(routes, post_routes)


# --------------------------------------------------------------------------- #
# Pure-function tests.
# --------------------------------------------------------------------------- #
def test_batch_window_arithmetic():
    start, close = runner.batch_window_utc(T0, 1)
    assert start == "2026-09-10T00:00:00Z"
    assert close == "2026-09-17T00:00:00Z"
    start3, close3 = runner.batch_window_utc(T0, 3)
    assert start3 == "2026-09-24T00:00:00Z"
    assert close3 == "2026-10-01T00:00:00Z"
    with pytest.raises(census.CensusError):
        runner.batch_window_utc(T0, 0)


def test_batch_windows_are_disjoint_and_covered():
    edges = set()
    for k in range(1, 13):
        start, close = runner.batch_window_utc(T0, k)
        assert start not in edges
        edges.add(start)
        assert runner.batch_window_utc(T0, k + 1)[0] == close


def test_census_query_registered_predicate():
    q = runner.build_census_query("2026-09-10T00:00:00Z", "2026-09-17T00:00:00Z")
    nodes = q["query"]["nodes"]
    attrs = {(n["parameters"]["attribute"], n["parameters"]["operator"])
             for n in nodes}
    assert ("exptl.method", "exact_match") in attrs
    assert ("rcsb_entry_info.resolution_combined", "less_or_equal") in attrs
    assert ("rcsb_accession_info.initial_release_date", "greater") in attrs
    assert ("rcsb_accession_info.initial_release_date", "less_or_equal") in attrs
    assert q["return_type"] == "entry"
    assert q["request_options"]["paginate"]["rows"] == runner.CENSUS_PAGE_ROWS


def test_parse_mappings_triples_and_offsets():
    payload = json.dumps(mapping_payload("aaaa", [("P12345", "A", 10),
                                                  ("Q99999", "B", 5)])).encode()
    triples = runner.parse_mappings(payload, "aaaa")
    assert len(triples) == 2
    by_acc = {t["accession"]: t for t in triples}
    assert by_acc["P12345"]["chain_id"] == "A"
    assert by_acc["P12345"]["res_map"][1] == 1
    assert by_acc["P12345"]["res_map"][10] == 10
    assert by_acc["Q99999"]["res_map"][5] == 5


def test_parse_mappings_author_numbering_offset():
    # author numbering offset from UniProt start (e.g. expression tag):
    # author 5..14 maps to unp 1..10
    doc = {"aaaa": {"UniProt": {"P12345": {"mappings": [{
        "chain": {"chain_id": "A", "author_asym_id": "A"}, "entity_id": 1,
        "unp_start": 1, "unp_end": 10,
        "start": {"author_residue_number": 5, "residue_number": 5},
        "end": {"author_residue_number": 14, "residue_number": 14},
    }]}}}}
    triples = runner.parse_mappings(json.dumps(doc).encode(), "aaaa")
    assert triples[0]["res_map"][5] == 1
    assert triples[0]["res_map"][14] == 10
    assert 4 not in triples[0]["res_map"]


def test_parse_mappings_malformed_is_typed():
    with pytest.raises(census.CensusError) as excinfo:
        runner.parse_mappings(b"not json", "aaaa")
    assert excinfo.value.code == runner.RunnerFailureCode.MAPPING_PARSE_FAILED.value


# --------------------------------------------------------------------------- #
# label_triple: coordinate geometry.
# --------------------------------------------------------------------------- #
def test_label_triple_identical_coordinates():
    n = 34
    coords = helix(n)
    ref = cif_text([("A", i + 1) for i in range(n)], 20.0, coords)
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, coords)
    triple = {"chain_id": "A", "entity_id": 1, "accession": "P12345",
              "res_map": {i + 1: i + 1 for i in range(n)}}
    out = runner.label_triple(triple, "aaaa", af, ref, "afsha", "refsha")
    assert out["status"] == "ok"
    assert len(out["rows"]) == n
    row0 = out["rows"][0]
    assert row0["uniprot_pos"] == 1
    assert row0["d_kabsch"] == 0.0  # identical structures superpose exactly
    assert row0["lddt"] == 1.0
    assert row0["plddt"] == 50.0
    z, _med, _iqr = pipe.bfactor_z([20.0 + i for i in range(n)])
    assert row0["b_factor_z"] == round(float(z[0]), 4)
    assert out["af_sha256"] == "afsha" and out["ref_sha256"] == "refsha"


def test_label_triple_insufficient_ref_chain():
    n = 34
    ref = cif_text([("A", i + 1) for i in range(10)], 20.0, helix(10))
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, helix(n))
    triple = {"chain_id": "A", "entity_id": 1, "accession": "P12345",
              "res_map": {i + 1: i + 1 for i in range(n)}}
    out = runner.label_triple(triple, "aaaa", af, ref, "afsha", "refsha")
    assert out["status"] == runner.RunnerFailureCode.REF_INSUFFICIENT_CA.value


def test_label_triple_insufficient_aligned_when_map_misses():
    n = 34
    ref = cif_text([("A", i + 100) for i in range(n)], 20.0, helix(n))  # offset auth nums
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, helix(n))
    triple = {"chain_id": "A", "entity_id": 1, "accession": "P12345",
              "res_map": {i + 1: i + 1 for i in range(n)}}  # expects auth 1..34
    out = runner.label_triple(triple, "aaaa", af, ref, "afsha", "refsha")
    assert out["status"] == runner.RunnerFailureCode.INSUFFICIENT_ALIGNED.value


def test_label_triple_duplicate_positions_deduped():
    n = 34
    coords = helix(n)
    ref = cif_text([("A", i + 1) for i in range(n)], 20.0, coords)
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, coords)
    # two author numbers mapped to the same UniProt position
    res_map = {i + 1: i + 1 for i in range(n)}
    res_map[34] = 1
    triple = {"chain_id": "A", "entity_id": 1, "accession": "P12345",
              "res_map": res_map}
    out = runner.label_triple(triple, "aaaa", af, ref, "afsha", "refsha")
    assert out["status"] == "ok"
    positions = [r["uniprot_pos"] for r in out["rows"]]
    # author 34 collides with author 1 (both -> unp 1); dedup drops one atom
    assert len(positions) == len(set(positions)) == n - 1


# --------------------------------------------------------------------------- #
# Scheduler POST discipline.
# --------------------------------------------------------------------------- #
def test_post_receipts_and_backoff_then_success():
    calls = {"n": 0}

    def flaky(url, payload, timeout_s):
        calls["n"] += 1
        if calls["n"] == 1:
            return census.HttpResponse(500, {}, b"boom")
        return census.HttpResponse(200, {}, census_body(["aaaa"]).encode())

    transport = FakeTransport()
    transport.post_json = flaky
    scheduler = runner.E422Scheduler(transport, "/tmp/e422_test_raw", sleeper=lambda s: None)
    result = scheduler.post(runner.RunnerRole.CENSUS, runner.SEARCH_URL,
                            runner.build_census_query("a", "b"))
    assert result.ok and calls["n"] == 2
    assert len(scheduler.receipts) == 2
    assert all(r["method"] == "POST" for r in scheduler.receipts)
    assert scheduler.receipts[0]["status"] == 500
    assert scheduler.receipts[0]["next_backoff_s"] == pytest.approx(2.0)
    assert scheduler.receipts[1]["status"] == 200
    assert result.content == census_body(["aaaa"]).encode()


def test_post_429_exhaustion_is_typed():
    transport = FakeTransport()
    transport.post_json = lambda url, payload, timeout_s: census.HttpResponse(
        429, {}, b"slow down")
    scheduler = runner.E422Scheduler(transport, "/tmp/e422_test_raw", sleeper=lambda s: None)
    result = scheduler.post(runner.RunnerRole.CENSUS, runner.SEARCH_URL,
                            runner.build_census_query("a", "b"))
    assert not result.ok
    assert result.failure_code == census.CensusFailureCode.HTTP_429_EXHAUSTED.value
    assert len(result.receipts) == scheduler._max_attempts


# --------------------------------------------------------------------------- #
# run_batch end-to-end (offline).
# --------------------------------------------------------------------------- #
def _checkpoint_of(transport, tmp_path, t0=T0, batch=1, manifest_path=None,
                   now_iso="2026-09-14T00:00:00Z"):
    return runner.run_batch(
        t0, batch, transport=transport,
        raw_dir=str(tmp_path / "raw"), out_dir=str(tmp_path / "out"),
        manifest_path=manifest_path, sleeper=lambda s: None, now_iso=now_iso)


def test_run_batch_happy_path(tmp_path):
    transport = standard_routes()
    checkpoint = _checkpoint_of(transport, tmp_path)
    assert checkpoint["schema"] == runner.RUNNER_SCHEMA
    assert checkpoint["census_entries"] == 1
    assert checkpoint["entries_with_mappings"] == 1
    assert checkpoint["n_rows"] == 34
    assert checkpoint["exclusions"] == {"mapping_failures": [],
                                        "manifest_failures": [],
                                        "triple_failures": []}
    assert checkpoint["verified_accessions"] == ["P12345"]
    # checkpoint digest verifies
    body = {k: v for k, v in checkpoint.items() if k != "sha256"}
    assert manifest_mod.canonical_sha256(body) == checkpoint["sha256"]
    # freeze receipt written, precedes labels (fewer requests than final)
    freeze_path = tmp_path / "out" / "e422_batch01_freeze.json"
    assert freeze_path.exists()
    freeze = json.loads(freeze_path.read_text())
    assert freeze["manifest_checkpoint_sha256"] == checkpoint[
        "manifest_checkpoint_sha256"]
    assert freeze["request_count"] < len(checkpoint["receipts"])
    assert freeze["frozen_utc"] == "2026-09-14T00:00:00Z"
    # freeze digest verifies
    freeze_body = {k: v for k, v in freeze.items() if k != "sha256"}
    assert manifest_mod.canonical_sha256(freeze_body) == freeze["sha256"]
    # rows carry the four channels
    row = checkpoint["rows"][0]
    for key in ("entry", "chain", "accession", "uniprot_pos", "plddt", "lddt",
                "d_kabsch", "b_factor_z"):
        assert key in row
    # every lane role appears in receipts exactly as registered
    roles = {r["role"] for r in checkpoint["receipts"]}
    assert roles == {role.value for role in runner.RunnerRole}


def test_run_batch_writes_and_reloads_manifest(tmp_path):
    manifest_path = str(tmp_path / "manifest" / "presence.json")
    transport = standard_routes()
    _checkpoint_of(transport, tmp_path, manifest_path=manifest_path)
    assert os.path.exists(manifest_path)
    # reload + digest check passes
    manifest = manifest_mod.ManifestBuilder.load(manifest_path)
    assert set(manifest.records) == {"P12345"}
    # second batch with a fresh transport: no re-admission fetch
    transport2 = standard_routes()
    _checkpoint_of(transport2, tmp_path, batch=2, manifest_path=manifest_path,
                   now_iso="2026-09-21T00:00:00Z")
    # census + mappings + entry files still hit; metadata admission NOT re-hit
    # (only the batch-time verification fetch remains)
    admission_urls = [c for c in transport2.calls
                      if runner.METADATA_URL.format(acc="P12345") in c[1]]
    assert len(admission_urls) == 1  # verification fetch only


def test_run_batch_absent_accession_fail_closed(tmp_path):
    entry = "aaaa"
    ref = cif_text([("A", i + 1) for i in range(34)], 20.0, helix(34))
    af = cif_text([("A", i + 1) for i in range(34)], 50.0, helix(34))
    metadata = json.dumps(af_metadata("P12345"))
    child_url = "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif"
    transport = FakeTransport(
        {
            runner.MAPPINGS_URL.format(entry=entry): (
                200, json.dumps(mapping_payload(entry, [("P12345", "A", 34),
                                                         ("P00000", "B", 34)]))),
            runner.METADATA_URL.format(acc="P12345"): (200, metadata),
            runner.METADATA_URL.format(acc="P00000"): (404, "{}"),
            runner.ENTRY_FILES_URL.format(entry=entry): (200, ref),
            child_url: (200, af),
            "https://alphafold.ebi.ac.uk/files/AF-P00000-F1-model_v6.cif": (200, af),
        },
        {runner.SEARCH_URL: (200, census_body([entry]))},
    )
    checkpoint = _checkpoint_of(transport, tmp_path)
    codes = {f["code"] for f in checkpoint["exclusions"]["manifest_failures"]}
    assert "E422_ABSENT_FROM_SNAPSHOT" in codes
    triple_codes = {f["code"] for f in checkpoint["exclusions"]["triple_failures"]}
    assert "E422_ABSENT_FROM_SNAPSHOT" in triple_codes
    assert checkpoint["n_rows"] == 34  # P12345 chain A still labeled
    chains = {r["chain"] for r in checkpoint["rows"]}
    assert chains == {"A"}


def test_run_batch_pre_snapshot_release_fail_closed(tmp_path):
    entry = "aaaa"
    ref = cif_text([("A", i + 1) for i in range(34)], 20.0, helix(34))
    transport = FakeTransport(
        {
            runner.MAPPINGS_URL.format(entry=entry): (
                200, json.dumps(mapping_payload(entry, [("P12345", "A", 34)]))),
            runner.METADATA_URL.format(acc="P12345"): (
                200, json.dumps(af_metadata("P12345", created=T0[:10]))),
            runner.ENTRY_FILES_URL.format(entry=entry): (200, ref),
        },
        {runner.SEARCH_URL: (200, census_body([entry]))},
    )
    checkpoint = _checkpoint_of(transport, tmp_path)
    codes = {f["code"] for f in checkpoint["exclusions"]["manifest_failures"]}
    assert "E422_PRE_SNAPSHOT_RELEASE" in codes
    assert checkpoint["n_rows"] == 0
    # the frozen manifest is not persisted when nothing was admitted?  it is:
    # (empty-manifest save is still a registered checkpoint of the batch)
    assert "manifest_checkpoint_sha256" in checkpoint


def test_run_batch_mapping_failures_tallied(tmp_path):
    transport = standard_routes(entries=("aaaa", "bbbb"))
    transport.routes[runner.MAPPINGS_URL.format(entry="bbbb")] = (503, "nope")
    checkpoint = _checkpoint_of(transport, tmp_path)
    assert checkpoint["census_entries"] == 2
    assert checkpoint["entries_with_mappings"] == 1
    assert checkpoint["exclusions"]["mapping_failures"] == [
        {"entry": "bbbb", "code": "E422_MAPPING_UNAVAILABLE",
         "detail": "E420_HTTP_5XX_MAX_ATTEMPTS"}]


def test_run_batch_mapping_404_is_typed_unavailable(tmp_path):
    # AMENDED 2026-09-10 under e423 (DECISIONS same date): the shared
    # census scheduler now types any non-2xx as a transport-layer failure
    # (the e422 silent-ok fall-through was a defect — see e423_runner
    # docstring); a 404 mapping response surfaces as E422_MAPPING_UNAVAILABLE
    # instead of reaching the JSON parser as a misleading parse failure.
    # The frozen e422_runner.py module bytes are unchanged.
    transport = standard_routes(entries=("aaaa", "bbbb"))
    transport.routes[runner.MAPPINGS_URL.format(entry="bbbb")] = (404, "nope")
    checkpoint = _checkpoint_of(transport, tmp_path)
    assert checkpoint["entries_with_mappings"] == 1
    assert checkpoint["exclusions"]["mapping_failures"] == [
        {"entry": "bbbb", "code": "E422_MAPPING_UNAVAILABLE",
         "detail": "E420_HTTP_UNEXPECTED_STATUS"}]


def test_run_batch_census_cap_is_typed(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "CENSUS_MAX_PAGES", 2)

    class EndlessTransport(FakeTransport):
        def post_json(self, url, payload, timeout_s):
            body = json.dumps({"total_count": 10 ** 6,
                               "result_set": [{"identifier": f"aaaa{x}"} for x in range(500)]})
            return census.HttpResponse(200, {}, body.encode())

    transport = EndlessTransport()
    with pytest.raises(census.CensusError) as excinfo:
        _checkpoint_of(transport, tmp_path)
    assert excinfo.value.code == runner.RunnerFailureCode.CENSUS_RESULT_CAP.value


def test_run_batch_multi_chain_same_accession(tmp_path):
    # the registered chain rule: ALL mapped chains contribute rows
    entry = "aaaa"
    n = 34
    ref = cif_text([("A", i + 1) for i in range(n)] + [("B", i + 1) for i in range(n)],
                   20.0, helix(2 * n))
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, helix(n))
    metadata = json.dumps(af_metadata("P12345"))
    child_url = "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif"
    uniprot = {"P12345": {"mappings": [
        {"chain": {"chain_id": "A", "author_asym_id": "A"}, "entity_id": 1,
         "unp_start": 1, "unp_end": n,
         "start": {"author_residue_number": 1, "residue_number": 1},
         "end": {"author_residue_number": n, "residue_number": n}},
        {"chain": {"chain_id": "B", "author_asym_id": "B"}, "entity_id": 1,
         "unp_start": 1, "unp_end": n,
         "start": {"author_residue_number": 1, "residue_number": 1},
         "end": {"author_residue_number": n, "residue_number": n}},
    ]}}
    doc = {entry: {"UniProt": uniprot}}
    transport = FakeTransport(
        {
            runner.MAPPINGS_URL.format(entry=entry): (200, json.dumps(doc)),
            runner.METADATA_URL.format(acc="P12345"): (200, metadata),
            runner.ENTRY_FILES_URL.format(entry=entry): (200, ref),
            child_url: (200, af),
        },
        {runner.SEARCH_URL: (200, census_body([entry]))},
    )
    checkpoint = _checkpoint_of(transport, tmp_path)
    assert checkpoint["n_rows"] == 2 * n
    chains = sorted({r["chain"] for r in checkpoint["rows"]})
    assert chains == ["A", "B"]
    # both chains map to the same UniProt positions
    pos_a = sorted(r["uniprot_pos"] for r in checkpoint["rows"] if r["chain"] == "A")
    pos_b = sorted(r["uniprot_pos"] for r in checkpoint["rows"] if r["chain"] == "B")
    assert pos_a == pos_b


def test_run_batch_child_url_taken_from_manifest_not_predicted(tmp_path):
    # the child fetch must use the manifest record's v6_cif_url verbatim
    entry = "aaaa"
    n = 34
    ref = cif_text([("A", i + 1) for i in range(n)], 20.0, helix(n))
    af = cif_text([("A", i + 1) for i in range(n)], 50.0, helix(n))
    metadata = json.dumps(af_metadata("P12345"))
    child_url = "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif"
    transport = FakeTransport(
        {
            runner.MAPPINGS_URL.format(entry=entry): (
                200, json.dumps(mapping_payload(entry, [("P12345", "A", 34)]))),
            runner.METADATA_URL.format(acc="P12345"): (200, metadata),
            runner.ENTRY_FILES_URL.format(entry=entry): (200, ref),
            child_url: (200, af),
        },
        {runner.SEARCH_URL: (200, census_body([entry]))},
    )
    _checkpoint_of(transport, tmp_path)
    child_calls = [c for c in transport.calls if c[0] == "GET" and "files/AF-" in c[1]]
    assert child_calls == [("GET", child_url)]
