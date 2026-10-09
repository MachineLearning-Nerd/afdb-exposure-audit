"""Check the numbers quoted in the paper against the shipped data, offline.

Every check is printed with its class:

  RECOMPUTED  the value is computed here from row-level or request-level
              inputs shipped in data/derived (census ledgers, e420 receipts,
              e421 census rows, cached RCSB/UniProt/PDBe responses, the AFDB
              template table and the post-2022 AFDB model table), with code
              written for this script;
  RECONCILED  a stored output is compared with this script's recomputation,
              or a stored summary is compared with its own per-record fields
              (catches an output whose summary and records disagree);
  RE-READ     the value can only be read from a stored output, because the
              raw third-party files it was computed from are not shipped
              (9qj6 mmCIF/SIFTS for Figure 1, the raw SIFTS/mmCIF files
              behind V7, the raw files behind the independent label audit,
              and the biotite cross-check and excluded-residue counts in
              audits/label_crosscheck_biotite, whose per-segment row counts
              and medians are then reconciled with the shipped ledger).

The script imports nothing from code/experiments, code/verification or
code/reproduce_p5.py: the V1, V2/V2c, V5, V6, post-2022 summary, P5 and
K7PQ54 computations are re-derived here from the registration texts.
Percentages are computed from computed numerators. Values are compared with
the paper at the paper's printed precision and with stored outputs at 1e-12
(5e-5 where a stored value was computed from raw files and the ledger holds
4-decimal values). Dates, design parameters and the timeline in the paper's
appendix are not checked unless a stored receipt records them.

Blinding: pLDDT is read from the preregistered census ledgers only (a) as
the registered P5 eligibility filter (0 <= pLDDT < 100.01) and (b) as the
per-chain median pLDDT of the Table 2 screen and of the 9qj6 chain, both
reported in the paper. No pLDDT-binned or pLDDT-conditional coverage,
Mondrian or CQR quantity is computed for the census. pLDDT-binned error
rates are computed only for the exploratory post-2022 (e421) census, as in
the paper.

Usage: python -I code/paper_numbers.py [<derived_root>]   (default: data/derived)
Exit status 0 iff every expected check ran and passed; 1 if a check failed;
2 if an input is missing or unreadable, a section raised, or the number of
checks differs from EXPECTED_CHECKS.
"""
import gzip
import hashlib
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

EXPECTED_CHECKS = 236
ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/derived")
AUD = ROOT.resolve().parent.parent / "audits"   # <release>/audits next to <release>/data/derived
VER = "paper/arxiv/verification"
CUTOFF = "2018-04-30"
RESULTS = []          # (class, ok)
ERRORS = []           # section-level errors (missing input, exception)
MISSING = object()
S = {}                # state shared between sections


class MissingInput(Exception):
    pass


def load(rel, base=None):
    p = (base or ROOT) / rel
    if not p.is_file():
        raise MissingInput(f"missing input: {p}")
    try:
        if p.suffix == ".gz":
            return json.loads(gzip.open(p).read())
        return json.loads(p.read_text())
    except (OSError, ValueError) as e:
        raise MissingInput(f"unreadable input {p}: {e}") from e


def get(obj, *path):
    for k in path:
        if isinstance(obj, dict) and k in obj:
            obj = obj[k]
        elif isinstance(obj, list) and isinstance(k, int) and -len(obj) <= k < len(obj):
            obj = obj[k]
        else:
            return MISSING
    return obj


def close(a, b, tol):
    """Deep equality with an absolute tolerance on floats; MISSING never matches."""
    if a is MISSING or b is MISSING:
        return False
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if isinstance(a, float) or isinstance(b, float):
            return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol
        return a == b
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k], tol) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))
    return a == b


def check(cls, where, what, paper, value, source, tol=None):
    """Compare a number printed in the paper with the value obtained here."""
    if value is MISSING:
        ok = False
    elif tol is None:
        ok = paper == value
    else:
        ok = close(paper, value, tol)
    RESULTS.append((cls, ok))
    shown = "MISSING" if value is MISSING else value
    print(f"{'OK  ' if ok else 'DIFF'} | {cls:<10} | {where:<6} | {what:<62} | paper {str(paper):<12} "
          f"| got {str(shown)[:28]:<28} | {source}")


def reconcile(where, what, stored, recomputed, source, tol=1e-12):
    """Compare a stored output (or summary) with the recomputation (or its own records)."""
    ok = close(stored, recomputed, tol)
    RESULTS.append(("RECONCILED", ok))
    print(f"{'OK  ' if ok else 'DIFF'} | RECONCILED | {where:<6} | {what:<62} | {'agree' if ok else 'DISAGREE':<55} "
          f"| {source}")
    if not ok:
        def brief(x):
            return "MISSING" if x is MISSING else json.dumps(x, default=lambda o: "MISSING")[:300]
        print(f"       stored:     {brief(stored)}\n       recomputed: {brief(recomputed)}")


def pct(a, b, nd=1):
    return round(100.0 * a / b, nd)


def fin(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def canonical_sha256(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode("utf-8")).hexdigest()


# =========================================================================== ATLAS
def atlas():
    """ATLAS availability census and V1 channel recount (Sec. 4, 5)"""
    man = load("results/e420/source_availability_manifest.json")
    rows = man["rows"]
    n_rows, n_ent = len(rows), len({r["entry"] for r in rows})
    n_acc = len({a for r in rows for a in r["accessions"]})
    check("RECOMPUTED", "Sec4", "ATLAS entry-chain records", 1938, n_rows, "e420 manifest rows")
    check("RECOMPUTED", "Sec4", "ATLAS PDB entries", 1735, n_ent, "e420 manifest rows")
    check("RECOMPUTED", "Sec4", "ATLAS UniProt accessions", 943, n_acc, "e420 manifest rows")
    reconcile("Sec4", "manifest summary equals counts from its rows and requests", man.get("summary", MISSING),
              {"rows_attempted": n_rows, "unique_entries": n_ent, "unique_accessions": n_acc,
               "rows_succeeded": sum(r["status"] == "ok" for r in rows),
               "rows_typed_failed": sum(r["status"] != "ok" for r in rows),
               "request_count": len(man["requests"])}, "results/e420/source_availability_manifest.json")

    rec = load("results/e420/label_run_receipt.json")
    ovs = load("results/e420/overflow_set.json")
    ovt = load("results/e420/overflow_temporal_check.json")
    tf = rec["typed_failures"]
    admitted = {t["accession"] for t in tf}
    check("RECOMPUTED", "Sec4", "accessions with AFDB model + SIFTS mapping (600 + 269)", 869,
          len(admitted) + len(set(ovs["overflow_accessions"])), "label_run_receipt + overflow_set")
    # The admitted pool is enumerated through its typed failures. That is complete only because
    # no admitted row reached the label ledger (every admitted pair failed the strict temporal
    # rule). V1 relies on the same fact; check it instead of assuming it.
    reconcile("Sec5", "receipt: ledger_rows == 0 and every admitted accession is typed",
              [get(rec, "summary", "ledger_rows"), get(rec, "summary", "accessions_admitted")],
              [0, len(admitted)], "results/e420/label_run_receipt.json")
    reconcile("Sec5", "overflow set equals the overflow temporal check",
              sorted(set(ovs["overflow_accessions"])),
              sorted({r["accession"] for r in ovt["temporal_fail"] + ovt["temporal_pass"]}),
              "overflow_set.json vs overflow_temporal_check.json")
    pat = re.compile(r"modelCreatedDate=(\d{4}-\d{2}-\d{2})\S*\s*!<\s*release=(\d{4}-\d{2}-\d{2})")
    pairs, undated = {}, []
    for t in tf:
        m = pat.search(t.get("detail", ""))
        if m:
            pairs[(t["accession"], t["entry"])] = (m.group(1), m.group(2))
        else:
            undated.append(t["accession"])
    for r in ovt["temporal_fail"] + ovt["temporal_pass"]:
        pairs[(r["accession"], r["entry"])] = (r["model_created"][:10], r["release"][:10])
    by = defaultdict(list)
    for (a, _), v in pairs.items():
        by[a].append(v)
    accs = len(by)
    acc_t = sum(all(rel <= CUTOFF for _, rel in v) for v in by.values())
    acc_p = sum(all(rel <= mc for mc, rel in v) for v in by.values())
    ent_t = sum(rel <= CUTOFF for _, rel in pairs.values())
    ent_p = sum(rel <= mc for mc, rel in pairs.values())
    src = "label_run_receipt + overflow_temporal_check"
    check("RECOMPUTED", "Sec5", "undecodable AFDB record dropped", ["Q582G4"], undated, "label_run_receipt")
    check("RECOMPUTED", "Sec5", "accessions", 868, accs, src)
    check("RECOMPUTED", "Sec5", "accession-entry pairs", 935, len(pairs), src)
    check("RECOMPUTED", "Sec5", "training-eligible accessions (all refs <= 2018-04-30)", 791, acc_t, src)
    check("RECOMPUTED", "Sec5", "training-eligible accessions (%)", 91.1, pct(acc_t, accs), src)
    check("RECOMPUTED", "Sec5", "training-eligible pairs", 851, ent_t, src)
    check("RECOMPUTED", "Sec5", "training-eligible pairs (%)", 91.0, pct(ent_t, len(pairs)), src)
    check("RECOMPUTED", "Sec5", "accessions with >= 1 later reference", 77, accs - acc_t, src)
    check("RECOMPUTED", "Sec5", "accessions flagged by the AFDB creation-date rule", 863, acc_p, src)
    check("RECOMPUTED", "Sec5", "  ... as % of accessions", 99.4, pct(acc_p, accs), src)
    years = sorted(int(rel[:4]) for _, rel in pairs.values())
    check("RECOMPUTED", "Sec5", "dated references", 935, len(years), src)
    check("RECOMPUTED", "Sec5", "first / last release year", "1988-2023", f"{years[0]}-{years[-1]}", src)
    check("RECOMPUTED", "Sec5", "median release year", 2008, int(statistics.median(years)), src)
    check("RECOMPUTED", "Abs", "ATLAS training-eligible share, rounded (abstract, conclusion)", "91%",
          f"{round(100 * acc_t / accs)}%", src)
    check("RECOMPUTED", "Sec4", "ATLAS web resource queried on", "2026-09-04", man.get("started_utc", "")[:10],
          "e420 manifest started_utc")
    v1 = load(f"{VER}/out/v1_channel_recount.json")
    u = v1.get("union", {})
    names = ("accessions", "entries", "accessions_all_refs_pre_training_cutoff", "entries_training_channel",
             "accessions_all_refs_pre_model_date", "entries_template_channel", "accessions_with_post_cutoff_ref")
    reconcile("Sec5", "out/v1 union counts and year histogram equal recomputation",
              [[u.get(k, MISSING) for k in names], v1.get("release_year_histogram", MISSING)],
              [[accs, len(pairs), acc_t, ent_t, acc_p, ent_p, accs - acc_t],
               dict(sorted(Counter(str(y) for y in years).items()))], "out/v1_channel_recount.json")


# =========================================================================== census ledgers
def ledgers():
    """Preregistered census ledgers: defect scope, 9qj6, label checks (Table 2, Sec. 4, 7.1, 7.2, Fig. 2)"""
    digest_ok = []

    def ledger(sub):
        rows = []
        for b in (1, 2, 3, 4):
            ck = load(f"results/{sub}/batches/e422_batch0{b}_checkpoint.json")
            digest_ok.append(canonical_sha256({k: v for k, v in ck.items() if k != "sha256"}) == ck.get("sha256")
                             and ck.get("n_rows") == len(ck["rows"]))
            rows += ck["rows"]
        return rows

    orig, corr = ledger("e422"), ledger("e427")
    S["orig"], S["corr"] = orig, corr
    reconcile("Tab2", "8 checkpoints: stored canonical SHA-256 and n_rows match content", [True] * 8, digest_ok,
              "results/{e422,e427}/batches")

    def segs(rows):
        s = defaultdict(lambda: ([], []))
        for r in rows:
            if fin(r.get("lddt")):
                k = (r["entry"], r["chain"], r["accession"])
                s[k][0].append(r["lddt"])
                if fin(r.get("plddt")):
                    s[k][1].append(r["plddt"])
        return {k: (statistics.median(v[0]), statistics.median(v[1]) if v[1] else None, len(v[0]))
                for k, v in s.items()}

    so, sc = segs(orig), segs(corr)
    for tag, rows, s, n_rows, n_seg, n_acc in (("original", orig, so, 117905, 488, 138),
                                               ("corrected", corr, sc, 123454, 495, 141)):
        src = f"results/{'e422' if tag == 'original' else 'e427'}/batches"
        check("RECOMPUTED", "Tab2", f"{tag}: residue rows", n_rows, len(rows), src)
        check("RECOMPUTED", "Tab2", f"{tag}: chain segments", n_seg, len(s), src)
        check("RECOMPUTED", "Tab2", f"{tag}: accessions", n_acc, len({r['accession'] for r in rows}), src)
    check("RECOMPUTED", "Sec4", "original accession-entry pairs", 272,
          len({(r["accession"], r["entry"]) for r in orig}), "results/e422/batches")

    def low_share(rows):
        y = [r["lddt"] for r in rows if fin(r.get("lddt"))]
        return sum(v < 0.60 for v in y) / len(y)
    lo_o, lo_c = low_share(orig), low_share(corr)
    check("RECOMPUTED", "Tab2", "original: residues with lDDT < 0.60 (%)", 31.9, round(100 * lo_o, 1), "ledger")
    check("RECOMPUTED", "Tab2", "corrected: residues with lDDT < 0.60 (%)", 1.1, round(100 * lo_c, 1), "ledger")

    def hcla(s):  # chain screen: >= 30 residues, median pLDDT >= 90, median lDDT < 0.4
        return [k for k, (ml, mp, n) in s.items() if n >= 30 and mp is not None and mp >= 90 and ml < 0.4]
    check("RECOMPUTED", "Tab2", "original: high-confidence low-accuracy chains", 116, len(hcla(so)),
          "ledger (per-chain median pLDDT)")
    check("RECOMPUTED", "Tab2", "  ... spanning accessions", 32, len({k[2] for k in hcla(so)}), "ledger")
    check("RECOMPUTED", "Tab2", "corrected: high-confidence low-accuracy chains", 0, len(hcla(sc)), "ledger")
    k7o = [v[0] for k, v in so.items() if k[2] == "K7PQ54"]
    k7c = [v[0] for k, v in sc.items() if k[2] == "K7PQ54"]
    check("RECOMPUTED", "Tab2", "K7PQ54 chains with median lDDT < 0.4 (original)", "29/38",
          f"{sum(x < 0.4 for x in k7o)}/{len(k7o)}", "ledger")
    check("RECOMPUTED", "Tab2", "K7PQ54 chains with median lDDT < 0.4 (corrected)", "0/38",
          f"{sum(x < 0.4 for x in k7c)}/{len(k7c)}", "ledger")
    both = set(so) & set(sc)
    changed = [k for k in both if so[k][0] != sc[k][0]]
    big = {k for k in both if abs(so[k][0] - sc[k][0]) > 0.2}
    S["damaged"] = big
    check("RECOMPUTED", "Sec7.2", "segments whose median lDDT changed", 188, len(changed), "ledgers")
    check("RECOMPUTED", "Sec7.2", "segments changed by > 0.2", 151, len(big), "ledgers")
    check("RECOMPUTED", "Sec7.2", "  ... across accessions", 40, len({k[2] for k in big}), "ledgers")

    def disagreeing_entries(s):
        by = defaultdict(list)
        for (e, c, a), v in s.items():
            by[(e, a)].append(v[0])
        return {e for (e, a), v in by.items() if len(v) > 1 and max(v) - min(v) >= 0.4}
    check("RECOMPUTED", "Sec7.2", "entries whose same-protein chains disagree by >= 0.4 (orig)", 13,
          len(disagreeing_entries(so)), "ledger")
    check("RECOMPUTED", "Sec7.2", "  ... after correction", 0, len(disagreeing_entries(sc)), "ledger")
    check("RECOMPUTED", "Sec7.2", "net change in rows", 5549, len(corr) - len(orig), "ledgers")
    key = lambda r: (r["entry"], r["chain"], r["accession"], r["uniprot_pos"])  # noqa: E731
    ko, kc = {key(r) for r in orig}, {key(r) for r in corr}
    check("RECOMPUTED", "Sec7.2", "residue positions only in corrected", 6465, len(kc - ko), "ledgers")
    check("RECOMPUTED", "Sec7.2", "residue positions only in original", 916, len(ko - kc), "ledgers")
    check("RECOMPUTED", "Sec7.2", "accessions added", 3,
          len({r["accession"] for r in corr} - {r["accession"] for r in orig}), "ledgers")
    check("RECOMPUTED", "Fig2", "segments only in the corrected ledger", 7, len(set(sc) - set(so)), "ledgers")
    fn = load(f"{VER}/out/figures_rev3_numbers.json")
    f2 = fn.get("f2_scope", {})
    mine = {"segments_in_both": len(both), "segments_original_only": len(set(so) - set(sc)),
            "segments_corrected_only": len(set(sc) - set(so)), "segments_changed_any": len(changed),
            "segments_abs_delta_gt_0.2": len(big), "accessions_abs_delta_gt_0.2": len({k[2] for k in big}),
            "accessions_in_both": len({k[2] for k in both}),
            "residues_original_finite_lddt": sum(fin(r.get("lddt")) for r in orig),
            "residues_corrected_finite_lddt": sum(fin(r.get("lddt")) for r in corr),
            "rows_original_total": len(orig), "rows_corrected_total": len(corr),
            "share_lddt_lt_0.60_original": lo_o, "share_lddt_lt_0.60_corrected": lo_c}
    reconcile("Fig2", "figures_rev3_numbers f2_scope equals recomputation", {k: f2.get(k, MISSING) for k in mine},
              mine, "out/figures_rev3_numbers.json")
    check("RECOMPUTED", "Sec7.2", "all 151 segments changed > 0.2 moved upward", True,
          all(sc[k][0] > so[k][0] for k in big), "ledgers")
    check("RECOMPUTED", "Fig2", "their original medians (min-max)", "0.17-0.70",
          f"{min(so[k][0] for k in big):.2f}-{max(so[k][0] for k in big):.2f}", "ledgers")
    check("RECOMPUTED", "Fig2", "their corrected medians at least", 0.87,
          math.floor(100 * min(sc[k][0] for k in big)) / 100, "ledgers")
    small = [k for k in changed if k not in big]
    check("RECOMPUTED", "Sec7.2", "other changed segments: count, max |change|", "37, 0.056",
          f"{len(small)}, {max(abs(so[k][0] - sc[k][0]) for k in small):.3f}", "ledgers")

    # census design as executed: windows, empty first batch, T0, end date
    cks = [load(f"results/e422/batches/e422_batch0{b}_checkpoint.json") for b in (1, 2, 3, 4)]
    from datetime import datetime  # noqa: PLC0415
    ts = lambda x: datetime.fromisoformat(x.replace("Z", "+00:00"))  # noqa: E731
    hours = [(ts(c["window"]["close"]) - ts(c["window"]["start"])).total_seconds() / 3600 for c in cks]
    contiguous = all(cks[i]["window"]["close"] == cks[i + 1]["window"]["start"] for i in range(3))
    check("RECOMPUTED", "Sec4", "batch windows (h), contiguous; batch 1 rows", "48/168/168/168 True 0",
          f"{'/'.join(str(round(h)) for h in hours)} {contiguous} {len(cks[0]['rows'])}", "results/e422/batches")
    check("RECOMPUTED", "Sec4", "T0 = first window start; last window closes", "2026-09-11 / 2026-10-04",
          f"{cks[0]['t0_utc'][:10] if cks[0]['t0_utc'] == cks[0]['window']['start'] else 'T0 != start'} / "
          f"{cks[3]['window']['close'][:10]}", "results/e422/batches")
    check("RECOMPUTED", "Sec5", "census released > 8 years after the 2018-04-30 cutoff", True,
          (ts(cks[0]["window"]["start"]) - ts("2018-04-30T23:59:59Z")).days > 8 * 365.25, "T0 vs cutoff")
    # 9zxa (Q92802, 8 chains): qualifying but never labelled. Original run: no protein mapping, so its
    # reference file was never fetched; corrected relabel: the reference fetch was absent (status None) for all chains
    def zxa(sub):
        ex = load(f"results/{sub}/batches/e422_batch03_checkpoint.json").get("exclusions", {})
        return sorted({(f.get("accession", "-"), f["code"].replace("E422_", "")) for k in ex for f in ex[k]
                       if f.get("entry") == "9zxa"}), sum(f.get("entry") == "9zxa" for k in ex for f in ex[k])
    z = [zxa("e422"), zxa("e427")]
    check("RECOMPUTED", "Sec7.2", "9zxa absent from both ledgers; why (original; corrected)",
          "absent; NO_PROTEIN_MAPPING; Q92802 x8 REF_CIF_FETCH_FAILED",
          f"{'absent' if not any(r['entry'] == '9zxa' for r in orig + corr) else 'PRESENT'}; "
          f"{'+'.join(c_ for _, c_ in z[0][0])}; " + "+".join(f"{a} x{z[1][1]} {c_}" for a, c_ in z[1][0]),
          "results/{e422,e427}/batches/e422_batch03 exclusions")
    rr = load("results/e427/replay_report.json")
    want = {str(b): [load(f"results/e422/batches/e422_batch0{b}_checkpoint.json")["sha256"],
                     load(f"results/e427/batches/e422_batch0{b}_checkpoint.json")["sha256"], True] for b in (2, 3, 4)}
    reconcile("Sec7.2", "replay report: original relabel reproduced every stored batch (digests, self-checks)",
              [rr.get("verdict"), {b: [get(d, "original_sha256"), get(d, "corrected_sha256"),
                                       all(get(d, "self_check").values()) if isinstance(get(d, "self_check"), dict)
                                       else MISSING] for b, d in rr.get("batches", {}).items()}],
              ["OK", want], "results/e427/replay_report.json")
    check("RECOMPUTED", "AppA", "replay and relabel date", "2026-10-07", rr.get("replay_utc", "")[:10],
          "results/e427/replay_report.json")

    # 9qj6 (Sec. 7.1, Fig. 1): what the ledgers determine is recomputed; the rest needs raw files
    f1 = fn.get("f1_mechanism_9qj6", {})
    aud = load(f"{VER}/audit/k7pq54_9qj6_audit.json")
    q = {(tag, ch): [r for r in rows if r["entry"] == "9qj6" and r["chain"] == ch]
         for tag, rows in (("o", orig), ("c", corr)) for ch in "AB"}
    first_o = min(r["uniprot_pos"] for r in q[("o", "A")])
    first_c = min(r["uniprot_pos"] for r in q[("c", "A")])
    check("RECOMPUTED", "Sec7.1", "register shift (first UniProt position, original - corrected)", 34,
          first_o - first_c, "ledgers")
    check("RECOMPUTED", "Sec7.1", "author residue 36 assigned UniProt (original first position)", 70, first_o,
          "results/e422/batches")
    reconcile("Sec7.1", "Figure 1 numbers agree with the ledgers (34, 70, 36, 507)",
              [f1.get(k, MISSING) for k in ("register_shift", "faulty_uniprot_for_first_atom",
                                            "first_observed_correct_uniprot", "sifts_unp_end")],
              [first_o - first_c, first_o, first_c, max(r["uniprot_pos"] for r in q[("c", "A")])],
              "out/figures_rev3_numbers.json vs ledgers")
    check("RE-READ", "Sec7.1", "SIFTS segment start (label 1 -> UniProt)", 35, f1.get("sifts_unp_start", MISSING),
          "out/figures_rev3_numbers.json (raw SIFTS not shipped)")
    check("RE-READ", "Sec7.1", "UniProt length (segment end)", 507, f1.get("sifts_unp_end", MISSING),
          "out/figures_rev3_numbers.json")
    med = {k: statistics.median(r["lddt"] for r in v) for k, v in q.items()}
    check("RECOMPUTED", "Sec7.1", "shifted median lDDT chain A", 0.198, round(med[("o", "A")], 3), "results/e422")
    check("RECOMPUTED", "Sec7.1", "shifted median lDDT chain B", 0.199, round(med[("o", "B")], 3), "results/e422")
    check("RECOMPUTED", "Sec7.1", "median pLDDT of shifted chain A", 98.4,
          round(statistics.median(r["plddt"] for r in q[("o", "A")]), 1), "results/e422 (chain median)")
    em = aud.get("e422_checkpoint_matches_summary", {})
    reconcile("Sec7.1", "9qj6 audit: shifted-chain medians and row counts equal the ledger",
              [get(em, "A", "median_lddt_ca"), get(em, "B", "median_lddt_ca"), get(em, "A", "n_rows"),
               get(em, "B", "n_rows")],
              [med[("o", "A")], med[("o", "B")], len(q[("o", "A")]), len(q[("o", "B")])],
              "audit/k7pq54_9qj6_audit.json", tol=5e-5)
    ch = aud.get("chains", {})
    check("RE-READ", "Sec7.1", "fixed: aligned residues A / B (raw files)", "472/471",
          f"{get(ch, 'A', 'n_aligned')}/{get(ch, 'B', 'n_aligned')}", "audit/k7pq54_9qj6_audit.json")
    mm = [get(ch, c, "n_sequence_mismatches") for c in "AB"]
    check("RE-READ", "Sec7.1", "fixed: sequence mismatches A + B (raw files)", 0,
          sum(mm) if all(fin(x) for x in mm) else MISSING, "audit json")
    check("RECOMPUTED", "Sec7.1", "fixed: median lDDT A (corrected ledger)", 0.993, round(med[("c", "A")], 3),
          "results/e427")
    check("RECOMPUTED", "Sec7.1", "fixed: median lDDT B (corrected ledger)", 0.995, round(med[("c", "B")], 3),
          "results/e427")
    reconcile("Sec7.1", "9qj6 audit (raw files) agrees with the corrected ledger (n, median)",
              [get(ch, "A", "n_aligned"), get(ch, "B", "n_aligned"),
               get(ch, "A", "median_lddt_ca"), get(ch, "B", "median_lddt_ca")],
              [len(q[("c", "A")]), len(q[("c", "B")]), med[("c", "A")], med[("c", "B")]],
              "audit json vs results/e427", tol=5e-5)

    # independent label audit (needs raw files): re-read, reconciled with the corrected ledger
    ils = load("results/e427/audit/independent_label_sample.json").get("independent_segment_sample", {})
    check("RE-READ", "Sec7.2", "independent audit: sampled segments reproduced exactly", "36 True",
          f"{ils.get('n_segments')} {ils.get('all_rows_match')}", "results/e427/audit/independent_label_sample.json")
    check("RE-READ", "Sec7.2", "  ... of which targeted / random", "26/10",
          f"{ils.get('targeted_segments')}/{ils.get('random_unchanged_segments')}", "same")
    rs = Counter("targeted" if "target" in s.get("reason", "") else "random" if "sample" in s.get("reason", "")
                 else "other" for s in ils.get("segments", []))
    reconcile("Sec7.2", "independent audit: targeted/random counts equal per-segment reasons",
              [ils.get("targeted_segments"), ils.get("random_unchanged_segments")], [rs["targeted"], rs["random"]],
              "independent_label_sample.json")
    cnt = Counter((r["entry"], r["chain"], r["accession"]) for r in corr)
    segl = ils.get("segments", [])
    n_led = [cnt[(s.get("entry"), s.get("chain"), s.get("accession"))] for s in segl]
    reconcile("Sec7.2", "independent audit: per-segment flags and row counts vs corrected ledger",
              [len(segl), [s.get("row_values_match") for s in segl], [s.get("corrected_checkpoint_rows") for s in segl],
               [s.get("independent_rows") for s in segl]],
              [ils.get("n_segments", MISSING), [True] * len(segl), n_led, n_led], "independent_label_sample.json")


# =========================================================================== V7
def trigger():
    """Trigger prevalence (Sec. 7.3, V7): raw SIFTS/mmCIF not shipped, so re-read and reconciled"""
    v7 = load(f"{VER}/out/v7_trigger_prevalence.json")
    core = v7.get("core_label_free", {})
    sg = core.get("segments", {})
    trig = v7.get("triggered_segments", [])
    pa = v7.get("per_accession", {})
    det = [s for s in trig if get(s, "offset_check", "status") == "DETERMINED"]
    nz = [s for s in det if s["offset_check"]["offset"] != 0]
    mapped, ntrig, nnz = sg.get("mapped", MISSING), sg.get("start_author_null_TRIGGER", MISSING), \
        sg.get("offset_nonzero", MISSING)
    check("RE-READ", "Sec7.3", "SIFTS segments mapped to roster proteins", 492, mapped, "out/v7")
    check("RE-READ", "Sec7.3", "segments starting on an unobserved residue", 304, ntrig, "out/v7")
    check("RE-READ", "Sec7.3", "  ... share (%)", 62, round(100 * ntrig / mapped) if fin(mapped) and fin(ntrig)
          else MISSING, "out/v7")
    check("RE-READ", "Sec7.3", "triggered segments with non-zero offset", 147, nnz, "out/v7")
    check("RE-READ", "Sec7.3", "  ... share of mapped (%)", 30, round(100 * nnz / mapped) if fin(mapped) and fin(nnz)
          else MISSING, "out/v7")
    check("RE-READ", "Sec7.3", "  ... proteins", 39, get(core, "accessions", "with_nonzero_offset_trigger"), "out/v7")
    offs = sorted(int(x) for x in sg.get("offset_distribution", {}) if int(x) != 0)
    check("RE-READ", "Sec7.3", "offset range", "-31..201", f"{offs[0]}..{offs[-1]}" if offs else MISSING, "out/v7")
    reconcile("Sec7.3", "V7 summary counts equal its per-segment records",
              [ntrig, sg.get("offset_distribution", MISSING), sg.get("offset_zero", MISSING), nnz,
               get(sg, "offset_status_counts", "DETERMINED"), get(core, "accessions", "with_nonzero_offset_trigger"),
               get(core, "accessions", "with_any_trigger")],
              [len(trig), {str(k): v for k, v in sorted(Counter(s["offset_check"]["offset"] for s in det).items())},
               len(det) - len(nz), len(nz), len(det), len({s["accession"] for s in nz}),
               len({s["accession"] for s in trig})], "out/v7_trigger_prevalence.json")
    reconcile("Sec7.3", "V7 summary counts equal its per-accession records",
              [mapped, ntrig, nnz, get(core, "accessions", "with_nonzero_offset_trigger"),
               core.get("roster_accessions", MISSING)],
              [sum(v["segments"] for v in pa.values()), sum(v["triggered"] for v in pa.values()),
               sum(v["triggered_offset_nonzero"] for v in pa.values()),
               sum(v["triggered_offset_nonzero"] > 0 for v in pa.values()),
               len({r["accession"] for r in S["orig"]})], "out/v7 per_accession")
    reconcile("Sec7.3", "V7 per-accession non-zero counts equal its per-segment records",
              {a: v["triggered_offset_nonzero"] for a, v in pa.items() if v["triggered_offset_nonzero"]},
              dict(Counter(s["accession"] for s in nz)), "out/v7")
    # link to the defect scope, against this script's own damage set from the ledgers
    damaged = S["damaged"]
    k = lambda s: (s["entry"], s["chain_id"], s["accession"])  # noqa: E731
    nz_in = sum(k(s) in damaged for s in nz)
    zero_only = len((damaged & {k(s) for s in det if s["offset_check"]["offset"] == 0}) - {k(s) for s in nz})
    check("RECONCILED", "Sec7.3", "non-zero-offset segments among the 151 changed > 0.2", 147, nz_in,
          "out/v7 segments x recomputed ledger damage set")
    check("RECONCILED", "Sec7.3", "other changed segments (zero offset, offset varies)", 4, zero_only, "same")
    lk = v7.get("link_to_fscope_NOT_IN_A4_V7_TEXT", {})
    reconcile("Sec7.3", "stored V7 link section equals recomputation",
              [lk.get(x, MISSING) for x in ("segments_abs_delta_gt_0.2", "triggered_nonzero_offset_segments_in_damaged_set",
                                            "damaged_keys_with_only_zero_offset_triggers",
                                            "damaged_keys_without_any_trigger")],
              [len(damaged), nz_in, zero_only, len(damaged - {k(s) for s in trig})], "out/v7")


# =========================================================================== V2 / V2c
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"


class Cache:
    """Read-only view of a request-keyed HTTP cache (file name = SHA-256 of 'METHOD\\nURL\\nBODY')."""

    def __init__(self, rel):
        self.dir = ROOT / rel
        if not self.dir.is_dir():
            raise MissingInput(f"missing cache directory: {self.dir}")
        self.log = []
        self.keys = []

    def get(self, method, url, body=None):
        k = hashlib.sha256(f"{method}\n{url}\n{body or ''}".encode()).hexdigest()
        p = self.dir / f"{k}.json"
        if not p.is_file():
            self.log.append((url.split("/")[2], "NOT_CACHED"))
            return None, None
        rec = json.loads(p.read_text())
        if hashlib.sha256((rec["body"] or "").encode()).hexdigest() != rec["sha256"] or rec["url"] != url:
            raise MissingInput(f"cache record {p} fails its own SHA-256/URL check")
        self.log.append((url.split("/")[2], rec["status"]))
        self.keys.append((str(self.dir / k), url.split("/")[2], rec["status"]))
        return rec["status"], rec["body"]


def hits(cache, nodes):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "entry", "request_options": {"return_all_hits": True, "results_verbosity": "compact"}}
    st, body = cache.get("POST", SEARCH, json.dumps(q, sort_keys=True))
    if st == 204:
        return 0
    if st != 200:
        return None
    return len(json.loads(body).get("result_set", []))


def exposure(cache, acc, seq, date):
    """Registered V2 rule: strongest of same accession, >= 95 %, >= 30 % identity released <= date."""
    dnode = {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_accession_info.initial_release_date", "operator": "less_or_equal",
        "value": f"{date}T23:59:59Z"}}
    anode = {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}}
    n = {"same_accession": hits(cache, [anode, dnode])}
    for lab, ident in (("seq95", 0.95), ("seq30", 0.30)):
        n[lab] = None if seq is None else hits(cache, [{"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 0.1, "identity_cutoff": ident, "sequence_type": "protein", "value": seq}}, dnode])
    if any(v is None for v in n.values()):
        n["class"] = "QUERY_FAILED"
    else:
        n["class"] = ("SAME_ACCESSION" if n["same_accession"] else "SEQ95" if n["seq95"]
                      else "SEQ30" if n["seq30"] else "NOVEL")
    return n


def classify(cache, accs, mc_of):
    out = {}
    for acc in accs:
        st, body = cache.get("GET", f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        seq = "".join(x.strip() for x in body.splitlines()[1:]) if st == 200 and body.startswith(">") else None
        mc = mc_of(acc)
        tr = exposure(cache, acc, seq, CUTOFF)
        if tr["class"] == "SAME_ACCESSION" and mc and mc >= CUTOFF:
            tb = "SAME_ACCESSION"  # a later date can only add hits
        else:
            tb = exposure(cache, acc, seq, mc)["class"] if mc else None
        out[acc] = {"training": tr, "template_bound": tb, "seq_len": len(seq) if seq else None}
    return out


def homology():
    """Pre-cutoff homology exposure (Table 1, Sec. 4, 6): V2/V2c replayed from the shipped caches"""
    census = load("results/e421/census_results.json")
    qual = load("results/e421/qualifying_entries.json")
    e422_mc = {r["accession"]: r["model_created_date"][:10]
               for r in load("results/e422/manifest_checkpoint.json")["records"]}
    e421_mc = {}
    for q in qual:
        e421_mc.setdefault(q["afdb_accession"], q["model_created_date"][:10])
    mc_of = lambda a: e422_mc.get(a) or e421_mc.get(a)  # noqa: E731  (registered lookup order)
    e421_ok = {v["accession"] for v in census.values() if v["status"] == "ok"}
    e422 = {r["accession"] for r in S["orig"]}
    universe = sorted(e421_ok | e422)
    c2 = Cache("data/paper_verification/http")
    cls2 = classify(c2, universe, mc_of)
    ids = sorted({q["pdb_id"] for q in qual})  # V2b release-date lookups, 50 ids per GraphQL request
    for i in range(0, len(ids), 50):
        qq = '{ entries(entry_ids: %s) { rcsb_id rcsb_accession_info { initial_release_date } } }' % json.dumps(
            [c.upper() for c in ids[i:i + 50]])
        c2.get("POST", GRAPHQL, json.dumps({"query": qq}))
    relabel = load("results/e427/e421_relabel/census_results.json")
    new_a = sorted({r["accession"] for r in S["corr"]} - set(universe))
    new_b = sorted({v["accession"] for v in relabel.values() if v["status"] == "ok"} - set(universe))
    c2c = Cache(f"{VER}/cache/v2c/http")
    cls2c = {"a_e427_e422_ledger_new": classify(c2c, new_a, mc_of),
             "b_e427_e421_relabel_new": classify(c2c, new_b, mc_of)}
    cls = {a: r["training"]["class"] for a, r in cls2.items()}
    for tab in cls2c.values():
        for a, r in tab.items():
            cls.setdefault(a, r["training"]["class"])
    S["class"] = cls

    rec = {a: r for a, r in cls2.items() if a in e422}
    tc = Counter(r["training"]["class"] for r in rec.values())
    src = "RCSB/UniProt cache replay"
    check("RECOMPUTED", "Tab1", "census proteins", 138, len(rec), "results/e422/batches")
    for c, n, p in (("SAME_ACCESSION", 78, 56.5), ("SEQ95", 9, 6.5), ("SEQ30", 41, 29.7), ("NOVEL", 10, 7.2)):
        check("RECOMPUTED", "Tab1", f"homology class {c}", f"{n} ({p}%)", f"{tc[c]} ({pct(tc[c], len(rec))}%)", src)
    close_n = tc["SAME_ACCESSION"] + tc["SEQ95"]
    check("RECOMPUTED", "Sec6", "same accession or >= 95% (abstract: 87, 63%)", "87 (63.0%)",
          f"{close_n} ({pct(close_n, len(rec))}%)", src)
    check("RECOMPUTED", "Sec6", "relative at >= 30%", "128 (92.8%)",
          f"{len(rec) - tc['NOVEL']} ({pct(len(rec) - tc['NOVEL'], len(rec))}%)", src)
    check("RECOMPUTED", "Sec6", "unresolved homology queries (census)", 0, tc["QUERY_FAILED"], src)
    check("RECOMPUTED", "Abs", "rounded shares: close 63%, none at 30% 7%, relative at 30% 93%", "63/7/93",
          f"{round(100 * close_n / len(rec))}/{round(100 * tc['NOVEL'] / len(rec))}/"
          f"{round(100 * (len(rec) - tc['NOVEL']) / len(rec))}", src)
    check("RECOMPUTED", "Sec4", "post-2022 census proteins in the homology search (64 + 22)", 86,
          len(e421_ok) + len(new_b), "e421 census + relabel")
    check("RECOMPUTED", "Sec7.1", "K7PQ54 UniProt sequence length", 507, cls2["K7PQ54"]["seq_len"], "UniProt cache")
    for acc, n in (("P00698", 730), ("P00918", 700), ("P61769", 693)):
        check("RECOMPUTED", "Sec6", f"{acc} same-accession entries before cutoff", n,
              cls2[acc]["training"]["same_accession"], src)

    v2 = load(f"{VER}/out/v2_prior_exposure.json")
    v2c = load(f"{VER}/out/v2c_extension.json")

    def view(r):
        t = r.get("training", {})
        tb = r.get("template_bound")
        return [t.get("class"), t.get("same_accession"), t.get("seq95"), t.get("seq30"),
                tb.get("class") if isinstance(tb, dict) else tb, r.get("seq_len")]
    reconcile("Tab1", "out/v2 per-accession classes and hit counts equal the replay",
              {a: view(r) for a, r in v2.get("accessions", {}).items()}, {a: view(r) for a, r in cls2.items()},
              "out/v2_prior_exposure.json")
    reconcile("Sec4", "out/v2c per-accession classes and hit counts equal the replay",
              {t: {a: view(r) for a, r in (get(v2c, "results", t, "accessions") or {}).items()}
               if isinstance(get(v2c, "results", t, "accessions"), dict) else MISSING for t in cls2c},
              {t: {a: view(r) for a, r in tab.items()} for t, tab in cls2c.items()}, "out/v2c_extension.json")
    check("RECOMPUTED", "Sec4", "accessions queried (both censuses)", 224, len(universe) + len(new_a) + len(new_b),
          "ledgers + e421 census + relabel")
    st = Counter(c2.log + c2c.log)
    search_keys = Counter(k for k, h, s_ in c2.keys + c2c.keys if h == "search.rcsb.org" and s_ == 200)
    check("RECOMPUTED", "Sec4", "HTTP 200 lookups that repeat an identical query", 2,
          sum(n - 1 for n in search_keys.values()), src)
    check("RECOMPUTED", "Sec4", "RCSB Search queries returning 200", 594, st[("search.rcsb.org", 200)], src)
    check("RECOMPUTED", "Sec4", "RCSB Search queries returning 204", 408, st[("search.rcsb.org", 204)], src)
    check("RECOMPUTED", "Sec4", "failed or uncached requests", 0,
          sum(n for (h, s), n in st.items() if s not in (200, 204)), src)
    stored = Counter((r["url"].split("/")[2], r["status"]) for r in v2.get("receipts", []) + v2c.get("receipts", []))
    reconcile("Sec4", "stored receipts: (host, status) multiset equals the replay",
              {f"{h} {s}": n for (h, s), n in stored.items()},
              {f"{h} {s}": n for (h, s), n in st.items()}, "out/v2 + out/v2c receipts")


# =========================================================================== V5
def templates():
    """AFDB template channel (Table 1, Sec. 4, 6): V5 replayed from the template table + PDBe/RCSB cache"""
    tab = load(f"{VER}/v5_afdb_template_table.json").get("accessions", {})
    roster = sorted({r["accession"] for r in S["orig"]})
    receipts = {}
    for b in (1, 2, 3, 4):
        for r in load(f"results/e422/batches/e422_batch0{b}_checkpoint.json")["receipts"]:
            receipts.setdefault(r["url"], r)
    reconcile("Tab1", "template table: one AFDB file per roster protein, SHA-256 = census receipt",
              {a: [get(tab, a, "sha256"), get(tab, a, "raw_name")] for a in sorted(set(tab) | set(roster))},
              {a: [get(receipts, get(tab, a, "url"), "sha256")] * 2 for a in roster},
              "v5_afdb_template_table.json vs results/e422/batches receipts")
    cache = Cache(f"{VER}/cache/v5/http")
    look = {}
    for pid in sorted({t["pdb_id"] for a in roster for t in tab[a]["templates"]}):
        st, body = cache.get("GET", f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pid.lower()}")
        keys = sorted(json.loads(body).get(pid.lower(), {}).get("UniProt") or {}) if st == 200 else None
        st2, body2 = cache.get("GET", f"https://data.rcsb.org/rest/v1/core/entry/{pid.upper()}")
        rel = None
        if st2 == 200:
            rel = (json.loads(body2).get("rcsb_accession_info", {}).get("initial_release_date") or "")[:10] or None
        look[pid] = (keys, rel)

    def maps(a, pid):
        keys = look[pid][0]
        return keys is not None and any(k == a or k.split("-")[0] == a for k in keys)
    cls = {}
    for a in roster:
        ts = tab[a]["templates"]
        failed = any(look[t["pdb_id"]][0] is None for t in ts)
        cls[a] = ("NONE" if not ts else "SAME_ACC" if any(maps(a, t["pdb_id"]) for t in ts)
                  else "UNRESOLVED" if failed else "OTHER")
    t = Counter(cls.values())
    src = "template table + PDBe/RCSB cache"
    for c, n, p in (("SAME_ACC", 45, 32.6), ("OTHER", 88, 63.8), ("UNRESOLVED", 5, 3.6)):
        check("RECOMPUTED", "Tab1", f"template class {c}", f"{n} ({p}%)", f"{t[c]} ({pct(t[c], len(roster))}%)", src)
    check("RECOMPUTED", "Abs", "same-protein template share, rounded", "33%", f"{round(100 * t['SAME_ACC'] / len(roster))}%",
          src)
    nt = [len(tab[a]["templates"]) for a in roster]
    check("RECOMPUTED", "Sec6", "templates per model (min / max)", "4/4", f"{min(nt)}/{max(nt)}", "template table")
    dated = {p: r for p, (_, r) in look.items() if r}
    check("RECOMPUTED", "Sec6", "earliest / latest template release", "1988 / 2021-01",
          f"{min(dated.values())[:4]} / {max(dated.values())[:7]}", "RCSB cache")
    check("RECOMPUTED", "Sec6", "templates released after 2021-02-15", 0,
          sum(r > "2021-02-15" for r in dated.values()), "RCSB cache")
    mg = Counter(m for a in roster for m in tab[a]["model_group_name"])
    check("RECOMPUTED", "Sec4", "model software", "AlphaFold Monomer v2.0 model: 138",
          ", ".join(f"{k}: {v}" for k, v in mg.items()), "template table (AFDB files)")
    hcls = S["class"]
    xt = defaultdict(Counter)
    for a in roster:
        xt[cls[a]][hcls[a]] += 1
    check("RECOMPUTED", "Sec6", "same-accession structure AND same-protein template", 43,
          xt["SAME_ACC"]["SAME_ACCESSION"], "V2 replay x V5 replay")
    check("RECOMPUTED", "Sec6", "same-accession structure, other templates only", 35, xt["OTHER"]["SAME_ACCESSION"],
          "same")
    extra = sorted(a for a in roster if cls[a] == "SAME_ACC" and hcls[a] != "SAME_ACCESSION")
    check("RECOMPUTED", "Sec6", "same-protein template without pre-cutoff same-accession", 2, len(extra), "same")
    yrs = sorted(min(look[x["pdb_id"]][1] for x in tab[a]["templates"] if maps(a, x["pdb_id"]))[:4] for a in extra)
    check("RECOMPUTED", "Sec6", "  ... their same-protein templates released in", "2019, 2020", ", ".join(yrs), "same")
    check("RECOMPUTED", "Sec6", "  ... Q9F0J8 has no relative at 30%", "NOVEL",
          hcls["Q9F0J8"] if "Q9F0J8" in extra else "not in set", "same")
    v5 = load(f"{VER}/out/v5_templates.json").get("results", {})
    per = v5.get("per_accession", {})
    reconcile("Tab1", "out/v5 per-accession template classes equal the replay",
              {a: get(per, a, "template_class") for a in roster}, cls, "out/v5_templates.json per_accession")
    reconcile("Tab1", "out/v5 class counts and cross-tab equal the replay",
              [v5.get("template_class_counts", MISSING), v5.get("crosstab_template_class_x_v2_training_class", MISSING)],
              [dict(t), {k: dict(v) for k, v in xt.items()}], "out/v5_templates.json")
    tl = lambda ts: [[x.get("pdb_id"), x.get("template_auth_asym_id")] for x in ts]  # noqa: E731
    reconcile("Tab1", "out/v5 template lists, failed lookups, release range equal table + cache",
              [{a: tl(get(per, a, "templates")) if isinstance(get(per, a, "templates"), list) else MISSING
                for a in roster}, v5.get("sifts_failed_pdb_ids", MISSING),
               v5.get("earliest_template_release_overall", MISSING), v5.get("latest_template_release_overall", MISSING),
               v5.get("model_group_names", MISSING)],
              [{a: tl(tab[a]["templates"]) for a in roster}, sorted(p for p, (k, _) in look.items() if k is None),
               min(dated.values()), max(dated.values()), dict(mg)], "out/v5_templates.json")


# =========================================================================== V6
def effect():
    """Exploratory exposure effect on the corrected post-2022 census (Sec. 6.3, V6)"""
    census = load("results/e427/e421_relabel/census_results.json")
    rows = defaultdict(lambda: ([], []))
    for eid, e in sorted(census.items()):
        if e["status"] != "ok":
            continue
        for lab in e["labels"]:
            if fin(lab.get("lddt")) and fin(lab.get("plddt")):
                rows[e["accession"]][0].append(lab["lddt"])
                rows[e["accession"]][1].append(lab["plddt"])
    per = {}
    for a, (l, p) in sorted(rows.items()):
        l, p = np.asarray(l, float), np.asarray(p, float)
        per[a] = [float(np.median(l)), float(np.mean(l < 0.60)), float(np.mean((p >= 90.0) & (l < 0.60))),
                  float(np.median(p))]
    cls = S["class"]
    gC = sorted(a for a in per if cls.get(a) in ("SAME_ACCESSION", "SEQ95"))
    gN = sorted(a for a in per if cls.get(a) in ("SEQ30", "NOVEL"))
    check("RECOMPUTED", "Sec6", "proteins with / without a close relative", "40/45", f"{len(gC)}/{len(gN)}",
          "relabel census + V2/V2c replay")
    XC, XN = np.array([per[a] for a in gC]), np.array([per[a] for a in gN])
    rng = np.random.default_rng(0)  # registered: 10,000 resamples within groups, seed 0, percentile CI
    iC = rng.integers(0, len(gC), size=(10_000, len(gC)))
    iN = rng.integers(0, len(gN), size=(10_000, len(gN)))
    boot = XC[iC].mean(axis=1) - XN[iN].mean(axis=1)
    eff = {}
    for j, q in enumerate(("median_lddt", "frac_lddt_lt_060", "frac_confident_error", "median_plddt")):
        lo, hi = np.percentile(boot[:, j], [2.5, 97.5])
        eff[q] = {"mean_CLOSE": float(XC[:, j].mean()), "mean_NOT_CLOSE": float(XN[:, j].mean()),
                  "diff_CLOSE_minus_NOT_CLOSE": float(XC[:, j].mean() - XN[:, j].mean()),
                  "ci95_percentile_bootstrap": [float(lo), float(hi)]}
    src = "relabel census, bootstrap over proteins"
    m = eff["median_lddt"]
    check("RECOMPUTED", "Sec6", "median lDDT CLOSE vs NOT_CLOSE", "0.975/0.959",
          f"{m['mean_CLOSE']:.3f}/{m['mean_NOT_CLOSE']:.3f}", src)
    check("RECOMPUTED", "Sec6", "  difference [95% CI]", "0.016 [0.003, 0.032]",
          f"{m['diff_CLOSE_minus_NOT_CLOSE']:.3f} [{m['ci95_percentile_bootstrap'][0]:.3f}, "
          f"{m['ci95_percentile_bootstrap'][1]:.3f}]", src)
    check("RECOMPUTED", "Sec6", "lDDT gap on the 0-100 scale vs pLDDT gap", "1.6 vs 3.8",
          f"{100 * eff['median_lddt']['diff_CLOSE_minus_NOT_CLOSE']:.1f} vs "
          f"{eff['median_plddt']['diff_CLOSE_minus_NOT_CLOSE']:.1f}", src)
    m = eff["median_plddt"]
    check("RECOMPUTED", "Sec6", "median pLDDT CLOSE vs NOT_CLOSE", "96.1/92.3",
          f"{m['mean_CLOSE']:.1f}/{m['mean_NOT_CLOSE']:.1f}", src)
    check("RECOMPUTED", "Sec6", "  difference [95% CI]", "3.8 [1.7, 6.0]",
          f"{m['diff_CLOSE_minus_NOT_CLOSE']:.1f} [{m['ci95_percentile_bootstrap'][0]:.1f}, "
          f"{m['ci95_percentile_bootstrap'][1]:.1f}]", src)
    m = eff["frac_confident_error"]
    check("RECOMPUTED", "Sec6", "confident-error share (%) CLOSE vs NOT_CLOSE", "0.12/0.13",
          f"{100 * m['mean_CLOSE']:.2f}/{100 * m['mean_NOT_CLOSE']:.2f}", src)
    check("RECOMPUTED", "Sec6", "  difference in pp [95% CI]", "-0.01 [-0.18, 0.17]",
          f"{100 * m['diff_CLOSE_minus_NOT_CLOSE']:.2f} [{100 * m['ci95_percentile_bootstrap'][0]:.2f}, "
          f"{100 * m['ci95_percentile_bootstrap'][1]:.2f}]", src)
    # AFDB model types of the post-2022 census (derived table of the October 2026 model files)
    mt = load("results/e427/e421_relabel/afdb_model_table.json").get("files", {})
    rcp = load("results/e427/e421_relabel/receipts.json")
    got = {r["url"]: r["sha256"] for r in rcp if r.get("status") == 200 and re.match(
        r"^https://alphafold\.ebi\.ac\.uk/files/AF-[^/]+\.cif$", r["url"])}
    reconcile("Sec4", "post-2022 model table: one row per downloaded AFDB file, SHA-256 = receipt",
              {u: get(v, "sha256") for u, v in mt.items()}, got, "afdb_model_table.json vs relabel receipts")
    kind = Counter("ColabFold" if any("ColabFold" in (g or "") for g in v["model_group_name"]) else
                   "AF2" if v["model_group_name"] == ["AlphaFold Monomer v2.0 model"] else "other" for v in mt.values())
    cf_acc = {a for v in mt.values() if any("ColabFold" in (g or "") for g in v["model_group_name"])
              for a in v["accessions_in_census"]}
    check("RECOMPUTED", "Sec4", "post-2022 model files: total / ColabFold v1.5.2 / AF Monomer v2.0", "93/6/87",
          f"{len(mt)}/{sum('ColabFold Monomer v1.5.2 model' in v['model_group_name'] for v in mt.values())}/"
          f"{kind['AF2']}", "afdb_model_table.json")
    check("RECOMPUTED", "Sec4", "ColabFold proteins among the 85 in the exposure comparison", 0,
          len(cf_acc & set(per)), "afdb_model_table.json x relabel census")
    v6 = load(f"{VER}/out/v6_exposure_effect.json").get("results", {})
    qn = ("median_lddt", "frac_lddt_lt_060", "frac_confident_error", "median_plddt")
    reconcile("Sec6", "out/v6 group sizes and per-protein quantities equal recomputation",
              [v6.get("group_sizes", MISSING), {a: [get(r, x) for x in qn] for a, r in v6.get("per_accession", {}).items()}],
              [{"CLOSE": len(gC), "NOT_CLOSE": len(gN), "excluded": len(per) - len(gC) - len(gN)}, per],
              "out/v6_exposure_effect.json")
    reconcile("Sec6", "out/v6 effects (means, differences, bootstrap CIs) equal recomputation",
              {q: {k: get(v6, "effects", q, k) for k in eff[q]} for q in eff}, eff, "out/v6_exposure_effect.json")


# =========================================================================== post-2022
def post2022():
    """Post-2022 census before and after the repair (Sec. 4, 7.3, 7.4; exploratory, not blinded)"""
    old = load("results/e421/census_results.json")
    old_labs = load("results/e421/labeled_residues_with_plddt.json")
    new = load("results/e427/e421_relabel/census_results.json")

    def summ(results, labs):
        ok = [r for r in results.values() if r["status"] == "ok"]
        grad = []
        for lo, hi in ((0, 50), (50, 70), (70, 80), (80, 90), (90, 100.01)):
            sel = [x["lddt"] for x in labs if fin(x.get("lddt")) and lo <= x["plddt"] < hi]
            grad.append({"bin": [lo, hi], "n": len(sel),
                         "err_rate_lddt_lt_0.60": sum(v < 0.60 for v in sel) / len(sel) if sel else None})
        f = [x for x in labs if fin(x.get("lddt"))]
        return {"entries": len(results), "status_counts": dict(Counter(r["status"] for r in results.values())),
                "ok_entries": len(ok), "ok_accessions": len({r["accession"] for r in ok}), "residue_rows": len(labs),
                "gradient": grad, "n_spearman": len(f),
                "spearman_plddt_lddt": float(spearmanr([x["plddt"] for x in f], [x["lddt"] for x in f]).statistic)}
    new_labs = [x for r in new.values() if r["status"] == "ok" for x in r["labels"]]
    o, c = summ(old, old_labs), summ(new, new_labs)
    reconcile("Sec7.4", "original e421 label list equals the labels of its ok entries", len(old_labs),
              sum(len(r.get("labels", [])) for r in old.values() if r["status"] == "ok"), "results/e421")
    tr = Counter(f"{old.get(k, {}).get('status')} -> {v['status']}" for k, v in new.items())
    src = "results/e421 + results/e427/e421_relabel census rows"
    check("RECOMPUTED", "Sec4", "post-2022 census entries", 168, o["entries"], "results/e421/census_results.json")
    check("RECOMPUTED", "Sec7.4", "original: rejected (too few aligned residues)", 69,
          o["status_counts"].get("insufficient_aligned"), "results/e421/census_results.json")
    check("RECOMPUTED", "Sec7.4", "original: model could not be downloaded", 3, o["status_counts"].get("af_fetch_fail"),
          "results/e421/census_results.json")
    check("RECOMPUTED", "Sec7.4", "of those 69, now labelled", 45, tr["insufficient_aligned -> ok"], src)
    check("RECOMPUTED", "Sec7.4", "labelled entries original -> corrected", "96 -> 139",
          f"{o['ok_entries']} -> {c['ok_entries']}", src)
    check("RECOMPUTED", "Sec7.4", "residue rows original -> corrected", "24457 -> 37531",
          f"{o['residue_rows']} -> {c['residue_rows']}", src)
    check("RECOMPUTED", "Sec7.4", "Spearman pLDDT-lDDT", "0.617 -> 0.643",
          f"{o['spearman_plddt_lddt']:.3f} -> {c['spearman_plddt_lddt']:.3f}", src)
    g = lambda s, lo: next(b for b in s["gradient"] if b["bin"][0] == lo)  # noqa: E731
    check("RECOMPUTED", "Sec7.4", "pLDDT>=90: share lDDT<0.60 (%) (n)", "0.58 (20777) -> 0.29 (31032)",
          f"{100 * g(o, 90)['err_rate_lddt_lt_0.60']:.2f} ({g(o, 90)['n']}) -> "
          f"{100 * g(c, 90)['err_rate_lddt_lt_0.60']:.2f} ({g(c, 90)['n']})", src)
    check("RECOMPUTED", "Sec7.4", "pLDDT<50: share lDDT<0.60 (%) (n)", "58 (177) -> 75 (314)",
          f"{100 * g(o, 0)['err_rate_lddt_lt_0.60']:.0f} ({g(o, 0)['n']}) -> "
          f"{100 * g(c, 0)['err_rate_lddt_lt_0.60']:.0f} ({g(c, 0)['n']})", src)
    ys = [x["lddt"] for x in new_labs if fin(x.get("lddt"))]
    check("RECOMPUTED", "Sec7.3", "corrected post-2022: residues below 0.60 lDDT (%)", 2.6,
          round(100 * sum(y < 0.60 for y in ys) / len(ys), 1), src)
    check("RECOMPUTED", "Sec4", "labelled entries / proteins: corrected; original", "139/85; 96/64",
          f"{c['ok_entries']}/{c['ok_accessions']}; {o['ok_entries']}/{o['ok_accessions']}", src)
    oldok = {k for k, v in old.items() if v["status"] == "ok"}
    newok = {k for k, v in new.items() if v["status"] == "ok"}
    check("RECOMPUTED", "Sec7.4", "previously labelled entries now rejected", "10pa, 10ps",
          ", ".join(sorted(oldok - newok)), src)
    ce = Counter(x["entry"] for x in old_labs if fin(x.get("lddt")) and x["plddt"] >= 90 and x["lddt"] < 0.60)
    check("RECOMPUTED", "Sec7.4", "original confident errors; held by 10ke", "121; 22",
          f"{sum(ce.values())}; {ce['10ke']}", "results/e421/labeled_residues_with_plddt.json")
    # decomposition: recomputed old entries vs the 45 new ones
    old_new = summ({k: new[k] for k in oldok & newok}, [x for k in sorted(oldok & newok) for x in new[k]["labels"]])
    added = summ({k: new[k] for k in newok - oldok}, [x for k in sorted(newok - oldok) for x in new[k]["labels"]])
    rise = c["spearman_plddt_lddt"] - o["spearman_plddt_lddt"]
    frac_rho = (old_new["spearman_plddt_lddt"] - o["spearman_plddt_lddt"]) / rise
    ce_share = lambda s_: g(s_, 90)["err_rate_lddt_lt_0.60"]  # noqa: E731
    frac_ce = (ce_share(o) - ce_share(old_new)) / (ce_share(o) - ce_share(c))
    check("RECOMPUTED", "Sec7.4", "share of the Spearman rise from relabelled old entries (> 0.5)", True,
          frac_rho > 0.5, f"{frac_rho:.2f} of the rise")
    check("RECOMPUTED", "Sec7.4", "share of the confident-error fall from old entries (> 0.5)", True,
          frac_ce > 0.5, f"{frac_ce:.2f} of the fall (pLDDT >= 90 share)")
    check("RECOMPUTED", "Sec7.4", "added residues from the 45 new entries (> 90% of net gain)", True,
          added["residue_rows"] / (c["residue_rows"] - o["residue_rows"]) > 0.9,
          f"{added['residue_rows']} of {c['residue_rows'] - o['residue_rows']}")
    check("RECOMPUTED", "Sec7.4", "pLDDT<50 rise not from old entries (old entries' share <= original)", True,
          g(old_new, 0)["err_rate_lddt_lt_0.60"] <= g(o, 0)["err_rate_lddt_lt_0.60"] < g(added, 0)["err_rate_lddt_lt_0.60"],
          f"old {g(old_new, 0)['err_rate_lddt_lt_0.60']:.2f}, new {g(added, 0)['err_rate_lddt_lt_0.60']:.2f}")
    rdates = {r["utc"][:7] for r in load("results/e427/e421_relabel/receipts.json")}
    check("RECOMPUTED", "Sec7.4", "post-2022 source files fetched again in", "2026-10", ", ".join(sorted(rdates)),
          "results/e427/e421_relabel/receipts.json")
    sm = load("results/e427/e421_relabel/summary.json")
    keys = ("entries", "status_counts", "ok_entries", "ok_accessions", "residue_rows", "gradient", "n_spearman",
            "spearman_plddt_lddt")
    reconcile("Sec7.4", "summary.json (original, corrected, transitions) equals recomputation",
              [{k: get(sm, "original", k) for k in keys}, {k: get(sm, "corrected", k) for k in keys},
               sm.get("status_transitions", MISSING)], [o, c, dict(tr)], "results/e427/e421_relabel/summary.json")


# =========================================================================== P5
def p5_eval(rows, open_interval=False, miss_level=None):
    """Registered P5 (e422 registration 3 / 3A.1): pooled split conformal over 200 accession-level splits."""
    acc, y = [], []
    for r in rows:
        if all(fin(r.get(k)) for k in ("lddt", "plddt", "d_kabsch", "b_factor_z")) and 0.0 <= r["plddt"] < 100.01:
            acc.append(r["accession"])
            y.append(float(r["lddt"]))
    y = np.asarray(y)
    names = sorted(set(acc))
    idx = {n: i for i, n in enumerate(names)}
    a = np.fromiter((idx[x] for x in acc), dtype=np.int64, count=len(acc))
    n = len(names)
    cov = {0.90: [], 0.95: []}
    miss = np.zeros(n)
    for seed in range(200):
        perm = np.random.default_rng(seed).permutation(n)
        t = n // 2
        c = (n - t) // 2
        role = np.empty(n, dtype=np.int8)
        role[perm[:t]], role[perm[t:t + c]], role[perm[t + c:]] = 0, 1, 2
        rr = role[a]
        centre = float(np.median(y[rr == 0]))
        sc = np.sort(np.abs(y[rr == 1] - centre))
        yt = y[rr == 2]
        for lv in (0.90, 0.95):
            k = math.ceil((sc.size + 1) * (1.0 - (1.0 - lv)))  # alpha_c = 1 - level, registered float form
            if k > sc.size:
                raise ValueError("calibration set too small for the order statistic")
            q = float(sc[k - 1])
            lo, hi = centre - q, centre + q
            inside = ((lo < yt) & (yt < hi)) if open_interval else ((lo <= yt) & (yt <= hi))
            cov[lv].append(float(inside.mean()))
            if miss_level == lv:  # V4's miss definition: |y - centre| > q, test proteins only
                np.add.at(miss, a[rr == 2], np.abs(yt - centre) > q)
    out = {}
    for lv, v in cov.items():
        v = np.asarray(v)
        mcse = float(v.std(ddof=1) / math.sqrt(v.size))
        thr = lv - 3.0 * mcse
        out[lv] = {"mean_cov": float(v.mean()), "mcse": mcse, "threshold": thr,
                   "outcome": "HOLD" if float(v.mean()) >= thr else "FAIL", "r_eff": int(v.size)}
    return out, int(y.size), {nm: float(miss[i]) for i, nm in enumerate(names)}


def coverage():
    """Registered coverage prediction P5 (Table 3, Sec. 7, 7.5): own re-implementation"""
    stored = load("p5_batch04_outcomes.json")
    res = {}
    for run, rows in (("original_e422", S["orig"]), ("corrected_e427", S["corr"])):
        res[run], n_elig, miss = p5_eval(rows, miss_level=0.90 if run == "original_e422" else None)
        if run == "original_e422":
            k7 = sum(r["accession"] == "K7PQ54" and all(fin(r.get(x)) for x in ("lddt", "plddt", "d_kabsch",
                                                                                 "b_factor_z"))
                     and 0 <= r["plddt"] < 100.01 for r in rows)
            check("RECOMPUTED", "Sec7.2", "K7PQ54 share of original eligible rows (%)", 14, round(100 * k7 / n_elig),
                  "results/e422/batches")
            top = max(miss, key=miss.get)
            share = miss[top] / sum(miss.values())
            check("RECOMPUTED", "Sec7.2", "K7PQ54 share of 0.90-level test misses (about half)", 0.50,
                  round(share, 2) if top == "K7PQ54" else f"top protein {top}", "same, 200 splits", tol=0.01)
            v4 = load(f"{VER}/out/v4_k7pq54.json")
            reconcile("Sec7.2", "out/v4: K7PQ54 rows, eligible rows, top miss share, top protein",
                      [v4.get("k7pq54_rows", MISSING), v4.get("n_rows", MISSING), get(v4, "miss_share_top", "1"),
                       get(v4, "top10", 0, "accession")], [k7, n_elig, share, top], "out/v4_k7pq54.json")
        reconcile("Tab3", f"{run}: stored n_rows_cumulative equals the ledger",
                  get(stored, run, "n_rows_cumulative"), len(rows), "p5_batch04_outcomes.json")
    paper = {("original_e422", 0.90): (0.8580, 0.8760, "FAIL"), ("original_e422", 0.95): (0.9152, 0.9328, "FAIL"),
             ("corrected_e427", 0.90): (0.8945, 0.8914, "HOLD"), ("corrected_e427", 0.95): (0.9456, 0.9449, "HOLD")}
    for (run, lv), (pm, pt, po) in paper.items():
        r = res[run][lv]
        lab = f"{run.split('_')[0]} P5-{int(round(lv * 100))}"
        check("RECOMPUTED", "Tab3", f"{lab} mean coverage", pm, round(r["mean_cov"], 4), "ledger, own P5")
        check("RECOMPUTED", "Tab3", f"{lab} threshold (nominal - 3 MCSE)", pt, round(r["threshold"], 4), "same")
        check("RECOMPUTED", "Tab3", f"{lab} outcome", po, r["outcome"], "same")
        o5 = get(stored, run, f"p5_{int(round(lv * 100))}", "outcome")
        reconcile("Tab3", f"{lab}: registered evaluator mean, MCSE, threshold, outcome (1e-12)",
                  [get(o5, "mean_cov"), get(o5, "mcse"), get(o5, "threshold"), get(o5, "outcome"), get(o5, "r_eff")],
                  [r["mean_cov"], r["mcse"], r["threshold"], r["outcome"], r["r_eff"]], "p5_batch04_outcomes.json")
    o90, o95 = res["original_e422"][0.90], res["original_e422"][0.95]
    c90, c95 = res["corrected_e427"][0.90], res["corrected_e427"][0.95]
    check("RECOMPUTED", "Tab3", "original margins at 0.90 / 0.95", "-0.0180/-0.0176",
          f"{o90['mean_cov'] - o90['threshold']:.4f}/{o95['mean_cov'] - o95['threshold']:.4f}", "own P5")
    check("RECOMPUTED", "Sec7.5", "original MCSE at 0.90 / 0.95", "0.0080/0.0057", f"{o90['mcse']:.4f}/{o95['mcse']:.4f}",
          "own P5")
    check("RECOMPUTED", "Sec7.5", "original MCSE more than twice the corrected (both levels)", True,
          o90["mcse"] > 2 * c90["mcse"] and o95["mcse"] > 2 * c95["mcse"],
          f"ratios {o90['mcse'] / c90['mcse']:.2f}, {o95['mcse'] / c95['mcse']:.2f}")
    check("RECOMPUTED", "Sec7.5", "0.90 margin about one MCSE; 0.95 margin below half an MCSE", "1 / True",
          f"{round((c90['mean_cov'] - c90['threshold']) / c90['mcse'])} / "
          f"{(c95['mean_cov'] - c95['threshold']) < 0.5 * c95['mcse']}", "own P5")
    check("RECOMPUTED", "Sec7.5", "margin at 0.90", 0.0031, round(c90["mean_cov"] - c90["threshold"], 4), "own P5")
    check("RECOMPUTED", "Sec7.5", "margin at 0.95", 0.0007, round(c95["mean_cov"] - c95["threshold"], 4), "own P5")
    check("RECOMPUTED", "Sec7.5", "MCSE at 0.95", 0.0017, round(c95["mcse"], 4), "own P5")
    check("RECOMPUTED", "Sec7.5", "corrected rows / accessions", "123454/141",
          f"{len(S['corr'])}/{len({r['accession'] for r in S['corr']})}", "results/e427/batches")
    op, _, _ = p5_eval(S["corr"], open_interval=True)
    # the two reimplementations cited in Sec. 7.5 (shipped outputs) against the registered evaluator
    st_ = {(run, lv): get(stored, run, f"p5_{int(round(lv * 100))}", "outcome") for run in res for lv in (0.90, 0.95)}
    pick = lambda d: [get(d, "mean_cov") if "mean_cov" in d else get(d, "mean"), get(d, "mcse"), get(d, "threshold"),
                      get(d, "outcome")] if isinstance(d, dict) else MISSING  # noqa: E731
    rec = load("results/e427/audit/phaseB_subagent/recompute_p5.json").get("summary", {})
    sen = load("results/e427/audit/phaseB_subagent/sensitivity_p5.json")
    ind = load("census_correction/independent_p5_audit.json", AUD).get("p5_recomputed", {})
    reconcile("Sec7.5", "reimplementation 1 (phase B) equals the evaluator; also reproduces the original failure",
              [pick(rec.get("0.9", {})), pick(rec.get("0.95", {})), pick(get(sen, "original_ledger_closed", "0.9")),
               pick(get(sen, "original_ledger_closed", "0.95"))],
              [pick(st_[("corrected_e427", 0.90)]), pick(st_[("corrected_e427", 0.95)]),
               pick(st_[("original_e422", 0.90)]), pick(st_[("original_e422", 0.95)])],
              "results/e427/audit/phaseB_subagent/{recompute,sensitivity}_p5.json")
    reconcile("Sec7.5", "reimplementation 2 (audit) equals the evaluator",
              [pick(ind.get("P5-90", {})), pick(ind.get("P5-95", {}))],
              [pick(st_[("corrected_e427", 0.90)]), pick(st_[("corrected_e427", 0.95)])],
              "audits/census_correction/independent_p5_audit.json")
    reconcile("Sec7.5", "open-interval sensitivity equals the shipped phase-B output",
              [pick(get(sen, "corrected_ledger_open_interval", "0.9")), pick(get(sen, "corrected_ledger_open_interval", "0.95"))],
              [pick(op[0.90]), pick(op[0.95])], "sensitivity_p5.json")
    rc_ = {s_: load(f"results/{s_}/evals/e422_batch04_execution_receipt.json") for s_ in ("e422", "e427")}
    check("RECOMPUTED", "AppA", "original evaluation ran; corrected re-evaluation completed", "2026-10-04..05; 2026-10-09",
          f"{rc_['e422'].get('started_utc', '')[:10]}..{rc_['e422'].get('completed_utc', '')[8:10]}; "
          f"{rc_['e427'].get('completed_utc', '')[:10]}", "results/{e422,e427}/evals execution receipts")
    check("RECOMPUTED", "Sec7.5", "boundary counted as uncovered: corrected 0.95 mean", 0.9453,
          round(op[0.95]["mean_cov"], 4), "own P5, open interval")
    check("RECOMPUTED", "Sec7.5", "  ... against threshold", 0.9449, round(op[0.95]["threshold"], 4), "same")
    pre = load("results/e427_repair/repaired_rows_batches01-04.json.gz")
    pr, n_pre, _ = p5_eval(pre)
    src = "results/e427_repair/repaired_rows_batches01-04.json.gz"
    check("RECOMPUTED", "Sec7.5", "preliminary (unregistered) repair: eligible residues", 122024, n_pre, src)
    check("RECOMPUTED", "Sec7.5", "preliminary repair P5 means (0.90 / 0.95)", "0.898/0.947",
          f"{pr[0.90]['mean_cov']:.3f}/{pr[0.95]['mean_cov']:.3f}", src)


def crosscheck():
    """Second label check (biotite), excluded residues, chance identity under a shift (Sec. 4, 7.2, 7.3)"""
    fl = load("label_crosscheck_biotite/full_ledger_out.json", AUD)
    a21 = load("label_crosscheck_biotite/analyse_21xg_out.json", AUD)
    ns = load("label_crosscheck_biotite/nonstd_out.json", AUD)
    okg = [r for r in fl if r.get("study_status") == "ok" and "error" not in r]
    exact = [r for r in okg if r.get("n_diff_gt_1e-4") == 0 and r.get("only_ledger") == 0 and r.get("only_biotite") == 0]
    src = "audits/label_crosscheck_biotite (needs raw files and PDBe updated mmCIF)"
    check("RE-READ", "Sec7.2", "biotite check: segments reproduced / compared", "494/495", f"{len(exact)}/{len(okg)}", src)
    mx = max(r["max_abs_diff_segment"] for r in exact) if exact else float("inf")
    check("RE-READ", "Sec7.2", "  ... every lDDT value within 5e-5", True, mx <= 5e-5 + 1e-12, f"max {mx:.3g}; {src}")
    x = next((r for r in okg if (r["entry"], r["chain"]) == ("21xg", "A")), {})
    check("RE-READ", "Sec7.2", "21xg A: values moved; largest move", "18; 0.054",
          f"{a21.get('n_lddt_values_differing_gt_1e-4')}; {a21.get('max_abs_lddt_difference_unrounded', 0):.3f}", src)
    check("RE-READ", "Sec7.2", "21xg A: median change; residues below 0.60 (ledger, biotite)", "0.003; 5, 5",
          f"{abs(a21.get('median_difference', 0)):.3f}; {a21.get('n_below_0.60_ledger')}, {a21.get('n_below_0.60_biotite')}",
          src)
    # reconcile the cross-check with the shipped corrected ledger
    by = defaultdict(list)
    for r in S["corr"]:
        by[(r["entry"].lower(), r["chain"], r["accession"])].append(r["lddt"])
    reconcile("Sec7.2", "biotite output: every compared segment's row count and median equal the ledger",
              {"/".join((r["entry"], r["chain"], r["acc"])): [r.get("ledger_n"), r.get("ledger_median")] for r in okg},
              {"/".join(k): [len(v), float(np.median(v))] for k, v in by.items()}, "full_ledger_out.json vs results/e427",
              tol=1e-12)
    led21 = by[("21xg", "A", "Q4KH59")]
    reconcile("Sec7.2", "21xg A detail agrees with the full run and the ledger",
              [a21.get("n_lddt_values_differing_gt_1e-4"), a21.get("max_abs_lddt_difference_unrounded"),
               a21.get("median_lddt_biotite"), a21.get("median_lddt_ledger"), a21.get("n_below_0.60_ledger")],
              [x.get("n_diff_gt_1e-4"), x.get("max_abs_diff_segment"), x.get("bt_median_segment"),
               float(np.median(led21)), sum(v < 0.60 for v in led21)], "analyse_21xg_out.json")
    check("RE-READ", "Sec4", "excluded residues: altloc other than A / modified other than MSE", "16/57",
          f"{ns.get('excluded_altloc_not_A_residues')}/{ns.get('excluded_modified_residues')}", src)
    pe = ns.get("per_entry", {})
    reconcile("Sec4", "excluded-residue totals equal their per-entry and per-component records",
              [ns.get("excluded_altloc_not_A_residues"), ns.get("excluded_modified_residues"),
               ns.get("excluded_modified_residues")],
              [sum(v["altloc_not_A"] for v in pe.values()), sum(v["modified"] for v in pe.values()),
               sum(ns.get("modified_by_component", {}).values())], "nonstd_out.json")
    # chance identity under a register shift, for the four chains used in the label-code fault tests
    # (code audit of 2026-10-09): shifts of 1 in each direction, 3, and the chain's true author offset
    seqs = {}
    for rel in ("data/paper_verification/http", f"{VER}/cache/v2c/http"):
        cache = Cache(rel)
        for acc in ("K7PQ54", "A0A4Q0WMG2", "P07342", "P00426"):
            if acc not in seqs:
                st_, body = cache.get("GET", f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
                if st_ == 200:
                    seqs[acc] = "".join(z.strip() for z in body.splitlines()[1:])
    fr = []
    for e, ch, acc, off in (("9qj6", "A", "K7PQ54", 34), ("9wvj", "A", "A0A4Q0WMG2", 14),
                            ("9y8y", "A", "P07342", 10), ("36cv", "E", "P00426", 0)):
        sq = seqs[acc]
        pos = [r["uniprot_pos"] for r in S["corr"] if (r["entry"], r["chain"], r["accession"]) == (e, ch, acc)]
        for k in {1, -1, -3, off, -off} - {0}:
            pr = [q for q in pos if 1 <= q + k <= len(sq)]
            fr.append(sum(sq[q - 1] == sq[q + k - 1] for q in pr) / len(pr))
    check("RECOMPUTED", "Sec7.3", "pairs matching by chance under a register shift (%)", "5-7",
          f"{round(100 * min(fr))}-{round(100 * max(fr))}", f"corrected ledger + UniProt cache, {len(fr)} shifts")


def rev6():
    """Revision 6 additions: census flow, sample composition, template chain names (Sec. 4, 6, 7.4)"""
    man = load("results/e420/source_availability_manifest.json")["rows"]
    check("RECOMPUTED", "Sec4", "ATLAS records not checked because of network errors", 4,
          sum((r.get("failure_code") or "") == "E420_TRANSPORT_ERROR" for r in man), "e420 manifest rows")
    check("RECOMPUTED", "Sec4", "ATLAS proteins in kept records (identity/coverage/AFDB rule)", 869,
          len({a for r in man if r["inclusion_status"] == "included" for a in r["accessions"]}), "e420 manifest rows")
    flow = {}
    for tag, sub in (("original", "e422"), ("corrected", "e427")):
        cks = [load(f"results/{sub}/batches/e422_batch0{b}_checkpoint.json") for b in (1, 2, 3, 4)]
        rows = [r for ck in cks for r in ck["rows"]]
        flow[tag] = (sum(ck["census_entries"] for ck in cks), sum(ck["entries_with_mappings"] for ck in cks),
                     len({r["entry"] for r in rows}), len({(r["accession"], r["entry"]) for r in rows}))
    src = "results/{e422,e427}/batches checkpoints"
    check("RECOMPUTED", "Tab2", "census entries enumerated (original/corrected)", "409/409",
          f"{flow['original'][0]}/{flow['corrected'][0]}", src)
    check("RECOMPUTED", "Tab2", "entries with a SIFTS mapping (original/corrected)", "364/365",
          f"{flow['original'][1]}/{flow['corrected'][1]}", src)
    check("RECOMPUTED", "Tab2", "entries with a scored chain segment (original/corrected)", "218/220",
          f"{flow['original'][2]}/{flow['corrected'][2]}", src)
    check("RECOMPUTED", "Tab2", "accession-entry pairs, corrected", 275, flow["corrected"][3], src)
    bt = load("label_crosscheck_biotite/full_ledger_out.json", AUD)
    ok = [r for r in bt if r.get("study_status") == "ok"]
    check("RE-READ", "Sec4", "scored chain segments made of one SIFTS segment", "495/495",
          f"{sum(r.get('n_segments') == 1 for r in ok)}/{len(ok)}", "audits/label_crosscheck_biotite")
    census = load("results/e427/e421_relabel/census_results.json")
    accs = {e["accession"] for e in census.values() if e["status"] == "ok"}
    cls = Counter(S["class"].get(a) for a in accs)
    check("RECOMPUTED", "Sec6", "proteins without a close relative: 30-95% relative / none", "38/7",
          f"{cls['SEQ30']}/{cls['NOVEL']}", "relabel census + V2/V2c replay")
    old_labs = load("results/e421/labeled_residues_with_plddt.json")
    new_labs = [x for e in census.values() if e["status"] == "ok" for x in e["labels"]]
    ce = lambda labs: sum(fin(x.get("lddt")) and x["plddt"] >= 90 and x["lddt"] < 0.60 for x in labs)  # noqa: E731
    check("RECOMPUTED", "Sec7.4", "confident errors (pLDDT>=90, lDDT<0.60) original -> corrected", "121 -> 90",
          f"{ce(old_labs)} -> {ce(new_labs)}", "results/e421 + results/e427/e421_relabel")
    v5 = get(load(f"{VER}/out/v5_templates.json"), "results", "template_class_counts_chain_level_NOT_PRESPECIFIED")
    check("RE-READ", "Sec6", "same-protein templates confirmed by chain name", 44,
          v5.get("SAME_ACC", MISSING) if isinstance(v5, dict) else MISSING, f"{VER}/out/v5_templates.json")
    st_, body = Cache(f"{VER}/cache/v5/http").get("GET", "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/6swu")
    m = json.loads(body).get("6swu", {}).get("UniProt", {}) if st_ == 200 else {}
    chains = {a: sorted({x.get("chain_id") for x in v.get("mappings", [])}) for a, v in m.items()}
    check("RECOMPUTED", "Sec6", "6SWU: accession and mapped chains", "Q5UE59: A,B,C,D,E,F",
          "; ".join(f"{a}: {','.join(c)}" for a, c in sorted(chains.items())) or MISSING, "PDBe SIFTS cache")
    old = load("results/e421/census_results.json")
    o_ok = {e["accession"] for e in old.values() if e["status"] == "ok"}
    check("RECOMPUTED", "Sec4", "post-2022 protein queried but not scored after correction", "Q2G0L4",
          ",".join(sorted(o_ok - accs)) or MISSING, "results/e421 + results/e427/e421_relabel")


def rev7():
    """Revision 7: homology exposure in two published evaluations (Sec. 6.4; registration e428)"""
    import csv
    t = load("results/e428/t_structures.json")
    a = load("results/e428/a_chains.json")
    st = load("results/e428/stats.json")
    cache = Cache("data/e428/http")
    close_ = ("SAME_ACCESSION", "SEQ95")
    src = "results/e428 per-chain classes"
    tc = Counter(r["class"] for r in t)
    classed = [r for r in t if r["class"] != "QUERY_FAILED"]
    check("RECOMPUTED", "Sec6.4", "Terwilliger: structures classed / with a close relative", "101/24",
          f"{len(classed)}/{sum(r['class'] in close_ for r in classed)}", src)
    check("RECOMPUTED", "Sec6.4", "Terwilliger: structures with no relative at >=30%", 49, tc["NOVEL"], src)
    failed = [r["pdb"] for r in t if r["class"] == "QUERY_FAILED"]
    st_, _ = cache.get("GET", "https://data.rcsb.org/rest/v1/core/entry/7DRH")
    check("RECOMPUTED", "Sec6.4", "Terwilliger: unclassed entry and its RCSB entry status", "7DRH: 404",
          f"{','.join(failed)}: {st_}", "results/e428 + e428 HTTP cache")
    reconcile("Sec6.4", "Terwilliger class counts: stats.json vs per-structure file",
              st["T"]["class_counts"], {k: tc.get(k, 0) for k in st["T"]["class_counts"]},
              "results/e428/stats.json")
    ac = Counter(r["class"] for r in a)
    n_close = sum(r["class"] in close_ for r in a)
    check("RECOMPUTED", "Sec6.4", "AlphaFlow test: chains with a close relative / total", "32/82",
          f"{n_close}/{len(a)}", src)
    check("RECOMPUTED", "Sec6.4", "AlphaFlow test: same protein / >=95% / >=30% / none", "28/4/25/25",
          f"{ac['SAME_ACCESSION']}/{ac['SEQ95']}/{ac['SEQ30']}/{ac['NOVEL']}", src)
    reconcile("Sec6.4", "AlphaFlow class counts: stats.json vs per-chain file",
              st["A"]["class_counts"], {k: ac.get(k, 0) for k in st["A"]["class_counts"]},
              "results/e428/stats.json")
    # Fusion constructs: chains with two UniProt accessions. The fusion partners are identified by name
    # (maltose-binding protein P0AEX9, T4 lysozyme D9IEF7); the other accession is the protein of interest.
    tags = {"P0AEX9", "D9IEF7"}
    fusions = [r for r in a if len(r["accessions"]) > 1]
    own = 0
    for r in fusions:
        n = 0
        for acc in (x for x in r["accessions"] if x not in tags):
            q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers"
                                 ".database_accession", "operator": "exact_match", "value": acc}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_accession_info.initial_release_date", "operator": "less_or_equal",
                    "value": f"{CUTOFF}T23:59:59Z"}}]},
                "return_type": "entry", "request_options": {"return_all_hits": True, "results_verbosity": "compact"}}
            s_, body = cache.get("POST", "https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q, sort_keys=True))
            if s_ not in (200, 204):
                raise MissingInput(f"e428 cache lacks the same-accession search for {acc}")
            n += 0 if s_ == 204 else len(json.loads(body).get("result_set", []))
        own += r["class"] in close_ and n > 0
    tag_only = sum(r["class"] in close_ for r in fusions) - own
    lo = n_close - tag_only
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: fusion chains / close only via the fusion partner", "3/2",
          f"{len(fusions)}/{tag_only}", "per-chain file + e428 HTTP cache")
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: close match via protein of interest / fusion-inclusive", "30 (37%)/32 (39%)",
          f"{lo} ({round(100 * lo / len(a))}%)/{n_close} ({round(100 * n_close / len(a))}%)", src)
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: test chains with train-split relative, same acc/>=95%/>=30%",
          "3/5/10", "/".join(str(sum(bool(r[k]) for r in a)) for k in
                              ("train_split_same_accession", "train_split_seq95", "train_split_seq30")), src)
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: >=95% train-split overlaps not via a fusion partner",
          "7buy_A,7e2s_A", ",".join(sorted(r["name"] for r in a if r["train_split_seq95"] and r not in fusions)), src)
    rows = {s: list(csv.DictReader((ROOT / f"data/e428/alphaflow/atlas_{s}.csv").open())) for s in ("train", "test")}
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: chains in the released train split file", 1266, len(rows["train"]),
          "data/e428/alphaflow/atlas_train.csv")
    check("RECOMPUTED", "Sec6.4", "AlphaFlow: test entries released on or before 2020-05-01", 40,
          sum(r["release_date"] <= "2020-05-01" for r in rows["test"]), "data/e428/alphaflow/atlas_test.csv")
    p1 = st["T"]["P1_rmsd"]
    check("RE-READ", "Sec6.4", "Terwilliger: median r.m.s.d. close (n) / other (n)", "0.847 (24)/0.956 (77)",
          f"{p1['median_close']} ({p1['n_close']})/{p1['median_other']} ({p1['n_other']})", "results/e428/stats.json")
    check("RE-READ", "Sec6.4", "Terwilliger: difference, 95% interval, permutation p", "-0.11 (-0.41 to 0.53), 0.38",
          f"{p1['diff_median_close_minus_other']:.2f} ({p1['ci95'][0]:.2f} to {p1['ci95'][1]:.2f}), "
          f"{p1['perm_p_two_sided']:.2f}", "results/e428/stats.json")
    check("RE-READ", "Sec6.4", "Terwilliger: median r.m.s.d. with no relative found", 0.943,
          st["T"]["S4_median_rmsd_by_class"]["NOVEL"], "results/e428/stats.json")


def run_section(fn):
    print(f"\n# {fn.__doc__}")
    try:
        fn()
    except MissingInput as e:
        ERRORS.append(str(e))
        print(f"ERROR | {e}")
    except Exception as e:  # noqa: BLE001 - reported and fatal for the exit status, never silent
        ERRORS.append(f"{fn.__name__}: {type(e).__name__}: {e}")
        print(f"ERROR | {fn.__name__}: {type(e).__name__}: {e}")


def main():
    if not ROOT.is_dir():
        print(f"derived root not found: {ROOT}")
        return 2
    needs = {trigger: "orig", homology: "orig", templates: "class", effect: "class", coverage: "orig",
             crosscheck: "orig"}
    for fn in (atlas, ledgers, trigger, homology, templates, effect, post2022, coverage, crosscheck, rev6, rev7):
        if fn in needs and needs[fn] not in S:
            ERRORS.append(f"{fn.__name__}: not run (a prerequisite section failed)")
            print(f"\nERROR | {fn.__name__} not run: a prerequisite section failed")
            continue
        run_section(fn)
    bad = sum(not ok for _, ok in RESULTS)
    total = Counter(c for c, _ in RESULTS)
    good = Counter(c for c, ok in RESULTS if ok)
    print("\nby class: " + ", ".join(f"{c} {good[c]}/{total[c]}" for c in ("RECOMPUTED", "RECONCILED", "RE-READ")))
    print(f"{len(RESULTS) - bad}/{len(RESULTS)} checks OK (expected {EXPECTED_CHECKS} checks)")
    if ERRORS:
        print(f"FAILED: {len(ERRORS)} section error(s): " + "; ".join(ERRORS))
        return 2
    if len(RESULTS) != EXPECTED_CHECKS:
        print(f"FAILED: ran {len(RESULTS)} checks, expected {EXPECTED_CHECKS}")
        return 2
    if bad:
        print(f"FAILED: {bad} check(s) differ")
        return 1
    print("ALL CHECKS OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
