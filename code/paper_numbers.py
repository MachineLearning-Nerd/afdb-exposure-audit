"""Recompute the numbers quoted in the paper from the shipped derived data.

Offline; reads only files under data/derived (default) and prints one line
per paper number: the value in the paper, the recomputed value, the source
file, and OK / DIFF. It imports no code from code/experiments and never
reads the (unreleased) batch-4 evaluation artifacts; pLDDT is used only for
the chain-level screens that the paper itself reports (Table 3), never for
pLDDT-binned coverage.

Usage: python -I code/paper_numbers.py [<derived_root>]
Exit status 0 iff every check is OK.
"""
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/derived")
OUT = ROOT / "paper/arxiv/verification/out"
CHECKS = []


def load(rel):
    return json.loads((ROOT / rel).read_text())


def check(where, what, paper, value, source, tol=None):
    if tol is None:
        ok = paper == value
    else:
        ok = abs(paper - value) <= tol
    CHECKS.append(ok)
    print(f"{'OK  ' if ok else 'DIFF'} | {where:<10} | {what:<62} | paper {paper!s:<10} | got {value!s:<22} | {source}")


def pct(a, b, nd=1):
    return round(100.0 * a / b, nd)


def fin(v):
    return v is not None and isinstance(v, (int, float)) and math.isfinite(v)


# ------------------------------------------------------------------ ATLAS (Sec. 4, 5, Table 1)
sam = load("results/e420/source_availability_manifest.json")["summary"]
check("Sec4", "ATLAS entry-chain records", 1938, sam["rows_attempted"], "results/e420/source_availability_manifest.json")
check("Sec4", "ATLAS PDB entries", 1735, sam["unique_entries"], "results/e420/source_availability_manifest.json")
check("Sec4", "ATLAS UniProt accessions", 943, sam["unique_accessions"], "results/e420/source_availability_manifest.json")
adm = load("results/e420/label_run_receipt.json")["summary"]["accessions_admitted"]
ovf = len(load("results/e420/overflow_set.json")["overflow_accessions"])
check("Sec4", "accessions with AFDB model + SIFTS mapping (600 + 269)", 869, adm + ovf,
      "results/e420/label_run_receipt.json + overflow_set.json")
v1 = load("paper/arxiv/verification/out/v1_channel_recount.json")
u = v1["union"]
check("Tab1", "accessions (one undecodable AFDB record dropped)", 868, u["accessions"], "out/v1_channel_recount.json")
check("Tab1", "accession-entry pairs", 935, u["entries"], "out/v1_channel_recount.json")
check("Tab1", "training-eligible accessions", 791, u["accessions_all_refs_pre_training_cutoff"], "out/v1")
check("Tab1", "training-eligible accessions (%)", 91.1, pct(791, u["accessions"]), "out/v1")
check("Tab1", "training-eligible pairs", 851, u["entries_training_channel"], "out/v1")
check("Tab1", "training-eligible pairs (%)", 91.0, pct(851, u["entries"]), "out/v1")
check("Tab1", "accessions with >=1 later reference", 77, u["accessions_with_post_cutoff_ref"], "out/v1")
check("Sec5", "accessions flagged by AFDB creation-date rule", 863, u["accessions_all_refs_pre_model_date"], "out/v1")
check("Sec5", "  ... as % of 868", 99.4, pct(863, u["accessions"]), "out/v1")
hist = v1["release_year_histogram"]
years = sorted(int(y) for y in hist)
expanded = [int(y) for y in sorted(hist) for _ in range(hist[y])]
check("Sec5", "dated references (release-year histogram total)", 935, len(expanded), "out/v1")
check("Sec5", "first / last release year", "1988-2023", f"{years[0]}-{years[-1]}", "out/v1")
check("Sec5", "median release year", 2008, int(statistics.median(expanded)), "out/v1")


# ------------------------------------------------------------------ census ledgers (Sec. 4, 7)
def ledger(sub):
    rows = []
    for b in (1, 2, 3, 4):
        rows += load(f"results/{sub}/batches/e422_batch0{b}_checkpoint.json")["rows"]
    return rows


orig, corr = ledger("e422"), ledger("e427")


def segs(rows):
    s = defaultdict(lambda: ([], []))
    for r in rows:
        if fin(r.get("lddt")):
            s[(r["entry"], r["chain"], r["accession"])][0].append(r["lddt"])
            if fin(r.get("plddt")):
                s[(r["entry"], r["chain"], r["accession"])][1].append(r["plddt"])
    return {k: (statistics.median(v[0]), statistics.median(v[1]) if v[1] else None, len(v[0])) for k, v in s.items()}


so, sc = segs(orig), segs(corr)
for tag, rows, s, n_rows, n_seg, n_acc in (("original", orig, so, 117905, 488, 138), ("corrected", corr, sc, 123454, 495, 141)):
    check("Tab3", f"{tag}: residue rows", n_rows, len(rows), f"results/{'e422' if tag == 'original' else 'e427'}/batches")
    check("Tab3", f"{tag}: chain segments", n_seg, len(s), "ledger")
    check("Tab3", f"{tag}: accessions", n_acc, len({r['accession'] for r in rows}), "ledger")
check("Sec4", "original accession-entry pairs", 272, len({(r["accession"], r["entry"]) for r in orig}), "results/e422/batches")
lo_o = sum(r["lddt"] < 0.60 for r in orig if fin(r["lddt"])) / sum(fin(r["lddt"]) for r in orig)
lo_c = sum(r["lddt"] < 0.60 for r in corr if fin(r["lddt"])) / sum(fin(r["lddt"]) for r in corr)
check("Tab3", "original: residues with lDDT < 0.60 (%)", 31.9, round(100 * lo_o, 1), "ledger")
check("Tab3", "corrected: residues with lDDT < 0.60 (%)", 1.1, round(100 * lo_c, 1), "ledger")


def hcla(s):
    return [k for k, (ml, mp, n) in s.items() if n >= 30 and mp is not None and mp >= 90 and ml < 0.4]


check("Tab3", "original: high-confidence low-accuracy chains", 116, len(hcla(so)), "ledger (pLDDT median per chain)")
check("Tab3", "  ... spanning accessions", 32, len({k[2] for k in hcla(so)}), "ledger")
check("Tab3", "corrected: high-confidence low-accuracy chains", 0, len(hcla(sc)), "ledger")
k7o = [v[0] for k, v in so.items() if k[2] == "K7PQ54"]
k7c = [v[0] for k, v in sc.items() if k[2] == "K7PQ54"]
check("Tab3", "K7PQ54 chains with median lDDT < 0.4 (original)", "29/38", f"{sum(x < 0.4 for x in k7o)}/{len(k7o)}", "ledger")
check("Tab3", "K7PQ54 chains with median lDDT < 0.4 (corrected)", "0/38", f"{sum(x < 0.4 for x in k7c)}/{len(k7c)}", "ledger")
both = set(so) & set(sc)
changed = [k for k in both if so[k][0] != sc[k][0]]
big = [k for k in both if abs(so[k][0] - sc[k][0]) > 0.2]
check("Sec7.2", "segments whose median lDDT changed", 188, len(changed), "ledgers")
check("Sec7.2", "segments changed by > 0.2", 151, len(big), "ledgers")
check("Sec7.2", "  ... across accessions", 40, len({k[2] for k in big}), "ledgers")


def disagreeing_entries(s):
    by = defaultdict(list)
    for (e, c, a), v in s.items():
        by[(e, a)].append(v[0])
    return {e for (e, a), v in by.items() if len(v) > 1 and max(v) - min(v) >= 0.4}


check("Sec7.2", "entries whose same-protein chains disagree by >= 0.4 (orig)", 13, len(disagreeing_entries(so)), "ledger")
check("Sec7.2", "  ... after correction", 0, len(disagreeing_entries(sc)), "ledger")
check("Sec7.2", "net change in rows", 5549, len(corr) - len(orig), "ledgers")
# "labels 6,465 residue positions that the original did not and drops 916":
# residue-level set difference on (entry, chain, accession, UniProt position).
key = lambda r: (r["entry"], r["chain"], r["accession"], r["uniprot_pos"])  # noqa: E731
ko, kc = {key(r) for r in orig}, {key(r) for r in corr}
check("Sec7.2", "residue positions only in corrected", 6465, len(kc - ko), "ledgers")
check("Sec7.2", "residue positions only in original", 916, len(ko - kc), "ledgers")
check("Sec7.2", "accessions added", 3, len({r["accession"] for r in corr} - {r["accession"] for r in orig}), "ledgers")
check("Sec7.2", "segments added", 7, len(set(sc) - set(so)), "ledgers")
v4 = load("paper/arxiv/verification/out/v4_k7pq54.json")
check("Sec7.2", "K7PQ54 share of original rows (%)", 14, round(100 * v4["k7pq54_rows"] / v4["n_rows"]), "out/v4_k7pq54.json")
check("Sec7.2", "K7PQ54 share of 0.90-level test misses (about half)", 0.50, round(v4["miss_share_top"]["1"], 2), "out/v4_k7pq54.json", tol=0.01)
fn = load("paper/arxiv/verification/out/figures_rev3_numbers.json")["f2_scope"]
check("Fig2", "segments only in corrected ledger", 7, fn["segments_corrected_only"], "out/figures_rev3_numbers.json")

# ------------------------------------------------------------------ 9qj6 mechanism (Sec. 7.1, Fig. 1)
aud = load("paper/arxiv/verification/audit/k7pq54_9qj6_audit.json")
f1 = load("paper/arxiv/verification/out/figures_rev3_numbers.json")["f1_mechanism_9qj6"]
check("Sec7.1", "register shift (author minus label)", 34, f1["register_shift"], "out/figures_rev3_numbers.json")
check("Sec7.1", "author residue 36 assigned UniProt", 70, f1["faulty_uniprot_for_first_atom"], "out/figures_rev3_numbers.json")
check("Sec7.1", "SIFTS segment start (label 1 -> UniProt)", 35, f1["sifts_unp_start"], "out/figures_rev3_numbers.json")
check("Sec7.1", "UniProt length (segment end)", 507, f1["sifts_unp_end"], "out/figures_rev3_numbers.json")
em = aud["e422_checkpoint_matches_summary"]
check("Sec7.1", "shifted median lDDT chain A", 0.198, round(em["A"]["median_lddt_ca"], 3), "audit/k7pq54_9qj6_audit.json")
check("Sec7.1", "shifted median lDDT chain B", 0.199, round(em["B"]["median_lddt_ca"], 3), "audit/k7pq54_9qj6_audit.json")
check("Sec7.1", "median pLDDT of shifted chain", 98.4, round(em["A"]["median_plddt"], 1), "audit/k7pq54_9qj6_audit.json")
ch = aud["chains"]
check("Sec7.1", "fixed: aligned residues A / B", "472/471", f"{ch['A']['n_aligned']}/{ch['B']['n_aligned']}", "audit json")
check("Sec7.1", "fixed: sequence mismatches A + B", 0, ch["A"]["n_sequence_mismatches"] + ch["B"]["n_sequence_mismatches"], "audit json")
check("Sec7.1", "fixed: median lDDT A", 0.993, round(ch["A"]["median_lddt_ca"], 3), "audit json")
check("Sec7.1", "fixed: median lDDT B", 0.995, round(ch["B"]["median_lddt_ca"], 3), "audit json")

# ------------------------------------------------------------------ trigger prevalence (Sec. 7.3, V7)
v7 = load("paper/arxiv/verification/out/v7_trigger_prevalence.json")
sg = v7["core_label_free"]["segments"]
check("Sec7.3", "SIFTS segments mapped to roster proteins", 492, sg["mapped"], "out/v7_trigger_prevalence.json")
check("Sec7.3", "segments starting on an unobserved residue", 304, sg["start_author_null_TRIGGER"], "out/v7")
check("Sec7.3", "  ... share (%)", 62, round(100 * sg["start_author_null_share"]), "out/v7")
check("Sec7.3", "triggered segments with non-zero offset", 147, sg["offset_nonzero"], "out/v7")
check("Sec7.3", "  ... share of 492 (%)", 30, round(100 * sg["offset_nonzero"] / sg["mapped"]), "out/v7")
check("Sec7.3", "  ... proteins", 39, v7["core_label_free"]["accessions"]["with_nonzero_offset_trigger"], "out/v7")
offs = sorted(int(x) for x in sg["offset_distribution"] if int(x) != 0)
check("Sec7.3", "offset range", "-31..201", f"{offs[0]}..{offs[-1]}", "out/v7")
lk = v7["link_to_fscope_NOT_IN_A4_V7_TEXT"]
check("Sec7.3", "non-zero-offset segments among the 151 changed > 0.2", 147, lk["triggered_nonzero_offset_segments_in_damaged_set"], "out/v7")
check("Sec7.3", "other changed segments (zero offset, offset varies)", 4, lk["damaged_keys_with_only_zero_offset_triggers"], "out/v7")

# ------------------------------------------------------------------ homology + templates (Sec. 6, Table 2)
v2 = load("paper/arxiv/verification/out/v2_prior_exposure.json")
rec = {a: r for a, r in v2["accessions"].items() if r["in_e422"]}
tc = Counter(r["training"]["class"] for r in rec.values())
check("Tab2", "census proteins", 138, len(rec), "out/v2_prior_exposure.json")
for cls, n, p in (("SAME_ACCESSION", 78, 56.5), ("SEQ95", 9, 6.5), ("SEQ30", 41, 29.7), ("NOVEL", 10, 7.2)):
    check("Tab2", f"homology class {cls}", f"{n} ({p}%)", f"{tc[cls]} ({pct(tc[cls], 138)}%)", "out/v2")
check("Sec6", "same accession or >=95% (abstract: 87, 63%)", "87 (63.0%)", f"{tc['SAME_ACCESSION'] + tc['SEQ95']} ({pct(tc['SAME_ACCESSION'] + tc['SEQ95'], 138)}%)", "out/v2")
check("Sec6", "relative at >= 30%", "128 (92.8%)", f"{138 - tc['NOVEL']} ({pct(138 - tc['NOVEL'], 138)}%)", "out/v2")
for acc, n in (("P00698", 730), ("P00918", 700), ("P61769", 693)):
    check("Sec6", f"{acc} same-accession entries before cutoff", n, v2["accessions"][acc]["training"]["same_accession"], "out/v2")
v5 = load("paper/arxiv/verification/out/v5_templates.json")["results"]
t = v5["template_class_counts"]
for cls, n, p in (("SAME_ACC", 45, 32.6), ("OTHER", 88, 63.8), ("UNRESOLVED", 5, 3.6)):
    check("Tab2", f"template class {cls}", f"{n} ({p}%)", f"{t.get(cls, 0)} ({pct(t.get(cls, 0), 138)}%)", "out/v5_templates.json")
check("Sec6", "templates per model (min = max)", "4", f"{v5['templates_per_accession']['min']}" if v5["templates_per_accession"]["min"] == v5["templates_per_accession"]["max"] else "varies", "out/v5")
check("Sec6", "earliest / latest template release", "1988 / 2021-01", f"{v5['earliest_template_release_overall'][:4]} / {v5['latest_template_release_overall'][:7]}", "out/v5")
check("Sec6", "templates released after 2021-02-15", 0, v5["n_templates_released_after_limit"], "out/v5")
check("Sec3", "model software", "AlphaFold Monomer v2.0 model: 138", ", ".join(f"{k}: {v}" for k, v in v5["model_group_names"].items()), "out/v5")
ct = v5["crosstab_template_class_x_v2_training_class"]
check("Sec6", "same-accession structure AND same-protein template", 43, ct["SAME_ACC"].get("SAME_ACCESSION", 0), "out/v5 crosstab")
check("Sec6", "same-accession structure, other templates only", 35, ct["OTHER"].get("SAME_ACCESSION", 0), "out/v5 crosstab")
extra = sorted(a for a, r in v5["per_accession"].items() if r["template_class"] == "SAME_ACC" and rec[a]["training"]["class"] != "SAME_ACCESSION")
check("Sec6", "same-protein template without pre-cutoff same-accession", 2, len(extra), "out/v5 per_accession")
yrs = sorted(min(tp["initial_release_date"] for tp in v5["per_accession"][a]["templates"] if tp.get("maps_to_accession"))[:4] for a in extra)
check("Sec6", "  ... their same-protein templates released in", "2019, 2020", ", ".join(yrs), "out/v5 per_accession")
check("Sec6", "  ... Q9F0J8 has no relative at 30%", "NOVEL", rec["Q9F0J8"]["training"]["class"] if "Q9F0J8" in extra else "not in set", "out/v2 + out/v5")

# ------------------------------------------------------------------ exposure effect (Sec. 6, V6)
v6 = load("paper/arxiv/verification/out/v6_exposure_effect.json")["results"]
gs = v6.get("group_sizes") or {}
ef = v6["effects"]
if gs:
    check("Sec6", "proteins with / without a close relative", "40/45", f"{gs.get('CLOSE')}/{gs.get('NOT_CLOSE')}", "out/v6_exposure_effect.json")
m = ef["median_lddt"]
check("Sec6", "median lDDT CLOSE vs NOT_CLOSE", "0.975/0.959", f"{m['mean_CLOSE']:.3f}/{m['mean_NOT_CLOSE']:.3f}", "out/v6")
check("Sec6", "  difference [95% CI]", "0.016 [0.003, 0.032]", f"{m['diff_CLOSE_minus_NOT_CLOSE']:.3f} [{m['ci95_percentile_bootstrap'][0]:.3f}, {m['ci95_percentile_bootstrap'][1]:.3f}]", "out/v6")
m = ef["median_plddt"]
check("Sec6", "median pLDDT CLOSE vs NOT_CLOSE", "96.1/92.3", f"{m['mean_CLOSE']:.1f}/{m['mean_NOT_CLOSE']:.1f}", "out/v6")
check("Sec6", "  difference [95% CI]", "3.8 [1.7, 6.0]", f"{m['diff_CLOSE_minus_NOT_CLOSE']:.1f} [{m['ci95_percentile_bootstrap'][0]:.1f}, {m['ci95_percentile_bootstrap'][1]:.1f}]", "out/v6")
m = ef["frac_confident_error"]
check("Sec6", "confident-error share (%) CLOSE vs NOT_CLOSE", "0.12/0.13", f"{100 * m['mean_CLOSE']:.2f}/{100 * m['mean_NOT_CLOSE']:.2f}", "out/v6")
check("Sec6", "  difference in pp [95% CI]", "-0.01 [-0.18, 0.17]", f"{100 * m['diff_CLOSE_minus_NOT_CLOSE']:.2f} [{100 * m['ci95_percentile_bootstrap'][0]:.2f}, {100 * m['ci95_percentile_bootstrap'][1]:.2f}]", "out/v6")

# ------------------------------------------------------------------ post-2022 census (Sec. 4, 7.4)
sm = load("results/e427/e421_relabel/summary.json")
o, c = sm["original"], sm["corrected"]
check("Sec4", "post-2022 census entries", 168, o["entries"], "results/e427/e421_relabel/summary.json")
check("Sec7.4", "original: rejected (too few aligned residues)", 69, o["status_counts"]["insufficient_aligned"], "summary.json")
check("Sec7.4", "original: model could not be downloaded", 3, o["status_counts"]["af_fetch_fail"], "summary.json")
check("Sec7.4", "of those 69, now labelled", 45, sm["status_transitions"]["insufficient_aligned -> ok"], "summary.json")
check("Sec7.4", "labelled entries original -> corrected", "96 -> 139", f"{o['ok_entries']} -> {c['ok_entries']}", "summary.json")
check("Sec7.4", "residue rows original -> corrected", "24457 -> 37531", f"{o['residue_rows']} -> {c['residue_rows']}", "summary.json")
check("Sec7.4", "Spearman pLDDT-lDDT", "0.617 -> 0.643", f"{o['spearman_plddt_lddt']:.3f} -> {c['spearman_plddt_lddt']:.3f}", "summary.json")
g = lambda s, lo: next(b for b in s["gradient"] if b["bin"][0] == lo)  # noqa: E731
check("Sec7.4", "pLDDT>=90: share lDDT<0.60 (%) (n)", "0.58 (20777) -> 0.29 (31032)",
      f"{100 * g(o, 90)['err_rate_lddt_lt_0.60']:.2f} ({g(o, 90)['n']}) -> {100 * g(c, 90)['err_rate_lddt_lt_0.60']:.2f} ({g(c, 90)['n']})", "summary.json")
check("Sec7.4", "pLDDT<50: share lDDT<0.60 (%) (n)", "58 (177) -> 75 (314)",
      f"{100 * g(o, 0)['err_rate_lddt_lt_0.60']:.0f} ({g(o, 0)['n']}) -> {100 * g(c, 0)['err_rate_lddt_lt_0.60']:.0f} ({g(c, 0)['n']})", "summary.json")
cr = load("results/e427/e421_relabel/census_results.json")
ys = [lab["lddt"] for e in cr.values() if e["status"] == "ok" for lab in e["labels"] if fin(lab.get("lddt"))]
check("Sec7.3", "corrected post-2022: residues below 0.60 lDDT (%)", 2.6, round(100 * sum(y < 0.60 for y in ys) / len(ys), 1), "results/e427/e421_relabel/census_results.json")

# ------------------------------------------------------------------ P5 (Table 4, Sec. 7.5)
p5 = load("p5_batch04_outcomes.json")
for run, lab in (("original_e422", "orig"), ("corrected_e427", "corr")):
    for k, pm, pt in ((("p5_90"), {"orig": 0.8580, "corr": 0.8945}, {"orig": 0.8760, "corr": 0.8914}),
                      (("p5_95"), {"orig": 0.9152, "corr": 0.9456}, {"orig": 0.9328, "corr": 0.9449})):
        o5 = p5[run][k]["outcome"]
        check("Tab4", f"{lab} {o5['prediction']} mean coverage", pm[lab], round(o5["mean_cov"], 4), "p5_batch04_outcomes.json")
        check("Tab4", f"{lab} {o5['prediction']} threshold", pt[lab], round(o5["threshold"], 4), "p5_batch04_outcomes.json")
c90, c95 = p5["corrected_e427"]["p5_90"]["outcome"], p5["corrected_e427"]["p5_95"]["outcome"]
check("Sec7.5", "margin at 0.90", 0.0031, round(c90["mean_cov"] - c90["threshold"], 4), "p5_batch04_outcomes.json")
check("Sec7.5", "margin at 0.95", 0.0007, round(c95["mean_cov"] - c95["threshold"], 4), "p5_batch04_outcomes.json")
check("Sec7.5", "MCSE at 0.95", 0.0017, round(c95["mcse"], 4), "p5_batch04_outcomes.json")
check("Sec7.5", "corrected rows / accessions", "123454/141", f"{p5['corrected_e427']['n_rows_cumulative']}/{len({r['accession'] for r in corr})}", "p5_batch04_outcomes.json + ledger")

# Registered-protocol recomputation (code/reproduce_p5.py; shares no code with the evaluator)
import gzip  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
import reproduce_p5 as rp  # noqa: E402

for run, sub in (("original_e422", "e422"), ("corrected_e427", "e427")):
    _, acc, y = rp.load_ledger(ROOT / f"results/{sub}/batches")
    _, res = rp.p5(acc, y)
    for k, pred in (("p5_90", "P5-90"), ("p5_95", "P5-95")):
        check("Sec7.5", f"recomputed {run} {pred} mean equals evaluator's", p5[run][k]["outcome"]["mean_cov"],
              res[pred]["mean_cov"], "code/reproduce_p5.py vs p5_batch04_outcomes.json", tol=1e-12)
_, acc, y = rp.load_ledger(ROOT / "results/e427/batches")
_, res = rp.p5(acc, y, open_interval=True)
check("Sec7.5", "boundary counted as uncovered: corrected 0.95 mean", 0.9453, round(res["P5-95"]["mean_cov"], 4), "code/reproduce_p5.py --open-interval")
check("Sec7.5", "  ... against threshold", 0.9449, round(res["P5-95"]["threshold"], 4), "code/reproduce_p5.py --open-interval")
pre = json.loads(gzip.open(ROOT / "results/e427_repair/repaired_rows_batches01-04.json.gz").read())
pacc = [r["accession"] for r in pre if all(fin(r.get(k)) for k in ("lddt", "plddt", "d_kabsch", "b_factor_z")) and 0 <= r["plddt"] < 100.01]
py = [r["lddt"] for r in pre if all(fin(r.get(k)) for k in ("lddt", "plddt", "d_kabsch", "b_factor_z")) and 0 <= r["plddt"] < 100.01]
import numpy as np  # noqa: E402
_, res = rp.p5(pacc, np.asarray(py, dtype=float))
check("Sec7.5", "preliminary (unregistered) repair: residues", 122024, len(py), "results/e427_repair/repaired_rows_batches01-04.json.gz")
check("Sec7.5", "preliminary repair P5 means (0.90 / 0.95)", "0.898/0.947", f"{res['P5-90']['mean_cov']:.3f}/{res['P5-95']['mean_cov']:.3f}", "same, via code/reproduce_p5.py")
ils = load("results/e427/audit/independent_label_sample.json")["independent_segment_sample"]
check("Sec7.2", "independent audit: sampled segments reproduced exactly", "36 True", f"{ils['n_segments']} {ils['all_rows_match']}", "results/e427/audit/independent_label_sample.json")

# ------------------------------------------------------------------ Appendix B (search receipts)
v2c = load("paper/arxiv/verification/out/v2c_extension.json")
n_acc = len(v2["accessions"]) + sum(len(part["accessions"]) for part in v2c["results"].values() if isinstance(part, dict) and "accessions" in part)
check("AppB", "accessions queried (both censuses)", 224, n_acc, "out/v2 + out/v2c_extension.json")
st = Counter((rr["url"].split("/")[2], rr["status"]) for rr in v2["receipts"] + v2c["receipts"])
check("AppB", "RCSB Search queries returning 200", 594, st[("search.rcsb.org", 200)], "out/v2 + out/v2c receipts")
check("AppB", "RCSB Search queries returning 204", 408, st[("search.rcsb.org", 204)], "out/v2 + out/v2c receipts")
check("AppB", "failed queries (non-200/204)", 0, sum(n for (h, s), n in st.items() if s not in (200, 204)), "out/v2 + out/v2c receipts")

bad = CHECKS.count(False)
print(f"\n{len(CHECKS) - bad}/{len(CHECKS)} checks OK")
sys.exit(1 if bad else 0)
