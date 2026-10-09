"""e427 label-mapping repair + program-pause guard (offline)."""

from __future__ import annotations

import hashlib
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "experiments"))

import e423_dispatch as dispatch  # noqa: E402
import e423_runner as old  # noqa: E402
import e427_label_repair as new  # noqa: E402
import e427_pause  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
AA = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAACVRESVPALL"
THREE = {v: k for k, v in {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}.items()}


def _walk(n=200, seed=7):
    import numpy as np
    rng = np.random.default_rng(seed)
    steps = rng.normal(size=(n, 3))
    steps = 3.8 * steps / np.linalg.norm(steps, axis=1, keepdims=True)
    return np.cumsum(steps, axis=0)


_TRACE = _walk()


def coord(p):
    """Deterministic, aperiodic CA random walk (3.8 A steps) keyed by UniProt
    position, so a register shift genuinely changes the geometry."""
    return tuple(float(v) for v in _TRACE[p - 1])


ATOM_TAGS = ["group_PDB", "id", "label_atom_id", "label_alt_id", "label_comp_id",
             "label_asym_id", "label_entity_id", "label_seq_id", "Cartn_x", "Cartn_y",
             "Cartn_z", "B_iso_or_equiv", "auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num"]


def cif(atoms, entity_seq=None):
    """atoms: (label_seq, auth_seq, aa1, xyz, b); entity_seq: [(num, aa1)]."""
    out = ["data_x", "#"]
    if entity_seq:
        out += ["loop_", "_entity_poly_seq.entity_id", "_entity_poly_seq.num",
                "_entity_poly_seq.mon_id", "_entity_poly_seq.hetero"]
        out += [f"1 {n} {THREE[a]} n" for n, a in entity_seq]
        out.append("#")
    out += ["loop_"] + [f"_atom_site.{t}" for t in ATOM_TAGS]
    for i, (lab, auth, aa, (x, y, z), b) in enumerate(atoms, 1):
        out.append(f"ATOM {i} CA . {THREE[aa]} A 1 {lab} {x:.3f} {y:.3f} {z:.3f} {b:.2f} {auth} A 1")
    out.append("#")
    return "\n".join(out) + "\n"


def af_cif(n=120):
    # AFDB F1: auth_seq_id == UniProt position; B column = pLDDT
    return cif([(p, p, AA[p - 1], coord(p), 95.0) for p in range(1, n + 1)])


def sifts(entry, segs):
    maps = [{"chain_id": "A", "struct_asym_id": "A", "entity_id": 1,
             "unp_start": us, "unp_end": ue,
             "start": {"residue_number": ls, "author_residue_number": a_s},
             "end": {"residue_number": le, "author_residue_number": a_e}}
            for (ls, le, us, ue, a_s, a_e) in segs]
    return json.dumps({entry: {"UniProt": {"P00001": {"mappings": maps}}}}).encode()


def ref_offset_case():
    """Reference covers UniProt 35..110 as label 1..76; label 1 is
    unobserved; author numbering = UniProt numbering (9qj6 pattern)."""
    atoms = [(lab, lab + 34, AA[lab + 33], coord(lab + 34), 20.0 + lab % 7)
             for lab in range(2, 77)]
    payload = sifts("1abc", [(1, 76, 35, 110, None, 110)])
    return payload, cif(atoms)


def test_old_runner_reproduces_register_shift():
    payload, ref = ref_offset_case()
    (t,) = old.parse_mappings(payload, "1abc")
    res = old.label_triple(t, "1abc", af_cif(), ref, "", "")
    # the frozen path keys label numbers (1..110) onto author atoms (36..110):
    # every residue lands 34 positions off
    assert res["status"] == "ok"
    assert {r["uniprot_pos"] for r in res["rows"]} == set(range(70, 111))
    assert sorted(r["lddt"] for r in res["rows"])[len(res["rows"]) // 2] < 0.5


def test_repair_maps_in_label_space():
    payload, ref = ref_offset_case()
    (t,) = new.parse_mappings_labelspace(payload, "1abc")
    assert t["map_method"] == "sifts_offset"
    assert t["res_map"][1] == 35 and t["res_map"][76] == 110
    res = new.label_triple_labelspace(t, "1abc", af_cif(), ref, "", "")
    assert res["status"] == "ok"
    assert res["n_identity_mismatch"] == 0
    assert [r["uniprot_pos"] for r in res["rows"]] == list(range(36, 111))
    assert all(r["lddt"] == 1.0 for r in res["rows"])
    assert all(r["d_kabsch"] < 1e-3 for r in res["rows"])


def test_identity_check_rejects_register_shift():
    payload, ref = ref_offset_case()
    (t,) = new.parse_mappings_labelspace(payload, "1abc")
    t = dict(t, res_map={k: v + 7 for k, v in t["res_map"].items()})  # forced shift
    res = new.label_triple_labelspace(t, "1abc", af_cif(), ref, "", "")
    assert res["status"] == new.E427FailureCode.SEQUENCE_IDENTITY_FAILED.value
    assert res["rows"] == []


def test_point_mutations_are_excluded_not_fatal():
    payload, _ = ref_offset_case()
    atoms = []
    for lab in range(2, 77):
        aa = AA[lab + 33]
        if lab in (10, 40):  # two engineered mutations
            aa = "W" if aa != "W" else "A"
        atoms.append((lab, lab + 34, aa, coord(lab + 34), 30.0))
    (t,) = new.parse_mappings_labelspace(payload, "1abc")
    res = new.label_triple_labelspace(t, "1abc", af_cif(), cif(atoms), "", "")
    assert res["status"] == "ok" and res["n_identity_mismatch"] == 2
    assert len(res["rows"]) == 73


def test_indel_segment_resolved_by_alignment():
    # label 1..70 covers UniProt 31..103 with a 3-residue internal deletion
    # (UniProt 61..63 absent from the construct)
    unp = [p for p in range(31, 104) if p not in (61, 62, 63)]
    entity = [(i + 1, AA[p - 1]) for i, p in enumerate(unp)]
    atoms = [(i + 1, 500 + i, AA[p - 1], coord(p), 25.0) for i, p in enumerate(unp)]
    payload = sifts("2def", [(1, 70, 31, 103, 500, 569)])
    (t,) = new.parse_mappings_labelspace(payload, "2def")
    assert t["map_method"] == "needs_alignment" and t["res_map"] is None
    res = new.label_triple_labelspace(t, "2def", af_cif(), cif(atoms, entity), "", "")
    assert res["status"] == "ok" and res["map_method"] == "needs_alignment"
    assert [r["uniprot_pos"] for r in res["rows"]] == unp
    assert all(r["lddt"] == 1.0 for r in res["rows"])


def test_missing_struct_asym_is_typed():
    payload = json.dumps({"3ghi": {"UniProt": {"P00001": {"mappings": [{
        "chain_id": "A", "entity_id": 1, "unp_start": 1, "unp_end": 40,
        "start": {"residue_number": 1, "author_residue_number": 1},
        "end": {"residue_number": 40, "author_residue_number": 40}}]}}}}).encode()
    with pytest.raises(Exception) as exc:
        new.parse_mappings_labelspace(payload, "3ghi")
    assert getattr(exc.value, "code", "") == old.RunnerFailureCode.MAPPING_PARSE_FAILED.value


# --------------------------------------------------------------------------- #
# Real cached bytes: 9qj6 / K7PQ54 (the audited case).  Skipped without cache.
# --------------------------------------------------------------------------- #
def _raw_for(url):
    ck = os.path.join(ROOT, "results/e422/batches/e422_batch04_checkpoint.json")
    if not os.path.exists(ck):
        return None
    for r in json.load(open(ck))["receipts"]:
        if r.get("url") == url and r.get("raw_name"):
            p = os.path.join(ROOT, "results/e422/raw", r["raw_name"])
            if os.path.exists(p):
                data = open(p, "rb").read()
                assert hashlib.sha256(data).hexdigest() == r["raw_name"]
                return data
    return None


def test_9qj6_real_bytes_repaired():
    mp = _raw_for("https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/9qj6")
    ref = _raw_for("https://www.ebi.ac.uk/pdbe/entry-files/9qj6.cif")
    af = None
    for b in (2, 3, 4):
        ck = os.path.join(ROOT, f"results/e422/batches/e422_batch{b:02d}_checkpoint.json")
        for r in json.load(open(ck))["receipts"] if os.path.exists(ck) else []:
            if "AF-K7PQ54-F1" in r.get("url", "") and r["url"].endswith(".cif") and r.get("raw_name"):
                af = open(os.path.join(ROOT, "results/e422/raw", r["raw_name"]), "rb").read()
    if mp is None or ref is None or af is None:
        pytest.skip("e422 raw cache for 9qj6/K7PQ54 not present")
    ref_t, af_t = ref.decode(), af.decode()
    for chain in ("A", "B"):
        o = next(t for t in old.parse_mappings(mp, "9qj6") if t["chain_id"] == chain)
        n = next(t for t in new.parse_mappings_labelspace(mp, "9qj6") if t["chain_id"] == chain)
        ro = old.label_triple(o, "9qj6", af_t, ref_t, "", "")
        rn = new.label_triple_labelspace(n, "9qj6", af_t, ref_t, "", "")
        med = lambda rows: sorted(r["lddt"] for r in rows)[len(rows) // 2]  # noqa: E731
        assert med(ro["rows"]) < 0.3            # the ledger's artifact
        assert rn["status"] == "ok" and rn["n_identity_mismatch"] == 0
        assert med(rn["rows"]) > 0.98           # audit F1: 0.993 / 0.995
        assert min(r["uniprot_pos"] for r in rn["rows"]) in (36, 37)


# --------------------------------------------------------------------------- #
# Pause guard.
# --------------------------------------------------------------------------- #
def _pause(tmp_path, monkeypatch, content='{"reason": "test", "paused_utc": "x"}'):
    p = tmp_path / "PAUSE.json"
    p.write_text(content)
    monkeypatch.setattr(e427_pause, "PAUSE_PATH", str(p))


def test_dispatch_refuses_while_paused(tmp_path, monkeypatch, capsys):
    _pause(tmp_path, monkeypatch)
    rc = dispatch.main(["--batch", "5", "--raw-dir", str(tmp_path / "r"),
                        "--out-dir", str(tmp_path / "o")])
    assert rc == 2
    assert "TYPED_FAILURE E427_PROGRAM_PAUSED" in capsys.readouterr().err
    assert not (tmp_path / "r").exists() and not (tmp_path / "o").exists()


def test_dispatch_dry_run_allowed_while_paused(tmp_path, monkeypatch):
    _pause(tmp_path, monkeypatch)
    assert dispatch.main(["--batch", "5", "--raw-dir", "r", "--out-dir", "o", "--dry-run"]) == 0


def test_malformed_pause_file_still_pauses(tmp_path, monkeypatch):
    _pause(tmp_path, monkeypatch, content="{not json")
    assert e427_pause.pause_reason() is not None


def test_evaluation_refuses_while_paused(tmp_path, monkeypatch, capsys):
    import e425_execreceipt
    _pause(tmp_path, monkeypatch)
    rc = e425_execreceipt.main(["--batch", "8", "--out-dir", str(tmp_path / "o"),
                                "--result-dir", str(tmp_path / "e")])
    assert rc == 2
    assert "TYPED_FAILURE E427_PROGRAM_PAUSED" in capsys.readouterr().err
    assert not (tmp_path / "e").exists()


# --------------------------------------------------------------------------- #
# Rewire (docs/e427_registration.md §4.2) and replay transport (§4.1).
# --------------------------------------------------------------------------- #
FROZEN = os.path.abspath(os.path.join(ROOT, "results", "e422", "batches"))


def test_dispatch_refuses_frozen_ledger_dir(capsys):
    rc = dispatch.main(["--batch", "5", "--raw-dir", "r", "--out-dir", FROZEN])
    assert rc == 2
    assert "E427_FROZEN_LEDGER_DIR" in capsys.readouterr().err


def test_dispatch_passes_e427_labeller(monkeypatch, tmp_path):
    seen = {}

    def fake_run_batch(*a, **kw):
        seen.update(kw)
        return {"n_rows": 0, "sha256": "x", "exclusions": {"mapping_failures": [],
                "manifest_failures": [], "triple_failures": []}, "census_entries": 0}
    monkeypatch.setattr(dispatch, "run_batch", fake_run_batch)
    dispatch.main(["--batch", "5", "--raw-dir", str(tmp_path / "r"),
                   "--out-dir", str(tmp_path / "o")])
    assert seen["parse_mappings_fn"] is new.parse_mappings_labelspace
    assert seen["label_triple_fn"] is new.label_triple_labelspace
    assert seen["labeller_id"] == dispatch.E427_LABELLER_ID


def test_evaluation_refuses_frozen_paths_and_defaults_to_e427(capsys):
    import e425_execreceipt
    rc = e425_execreceipt.main(["--batch", "8", "--out-dir", FROZEN])
    assert rc == 2 and "E427_FROZEN_LEDGER_DIR" in capsys.readouterr().err
    src = open(os.path.join(ROOT, "experiments", "e425_execreceipt.py")).read()
    assert 'default="results/e427/batches"' in src
    assert 'default="results/e427/evals"' in src


def test_replay_transport_order_and_hash(tmp_path):
    import e427_replay_relabel as rp
    body = b'{"ok": 1}'
    name = hashlib.sha256(body).hexdigest()
    (tmp_path / name).write_bytes(body)
    receipts = [{"method": "GET", "url": "u", "status": 429, "raw_name": None, "headers": {}},
                {"method": "GET", "url": "u", "status": 200, "raw_name": name, "headers": {}}]
    tr = rp.ReplayTransport(receipts, str(tmp_path))
    assert tr.get("u", 1).status == 429
    r = tr.get("u", 1)
    assert r.status == 200 and r.content == body and tr.leftover() == 0
    with pytest.raises(rp.ReplayError):
        tr.get("u", 1)
    (tmp_path / name).write_bytes(b"tampered")
    tr2 = rp.ReplayTransport(receipts[1:], str(tmp_path))
    with pytest.raises(rp.ReplayError):
        tr2.get("u", 1)
