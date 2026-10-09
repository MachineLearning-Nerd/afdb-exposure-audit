"""Independent audit of e428 (reviewer phase). Offline: reads only cached files.

Run: python3 -I experiments/audit_e428.py  -> prints JSON to stdout.
Written without reading experiments/e428_exposure_audit.py or docs/e428_results.md.
"""
import csv, glob, hashlib, json, math, os, random, re, statistics, sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "e428")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CUT_AF2 = "2018-04-30T23:59:59Z"
out = {}

# ---------------------------------------------------------------- xlsx
X = os.path.join(D, "terwilliger", "x", "xl")
xlsx_sha = hashlib.sha256(open(os.path.join(D, "terwilliger", "AlphaFoldCrystal_M.xlsx"), "rb").read()).hexdigest()
out["xlsx_sha256_ok"] = xlsx_sha == "7b3e1473b99e1575019d8c095d292d6a96de702c36580241cf85932e7d69d143"
sst = []
for si in ET.parse(os.path.join(X, "sharedStrings.xml")).getroot().findall("m:si", NS):
    sst.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))


def sheet_path(name):
    wb = ET.parse(os.path.join(X, "workbook.xml")).getroot()
    rid = None
    for s in wb.find("m:sheets", NS):
        if s.get("name") == name:
            rid = s.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
    rels = ET.parse(os.path.join(X, "_rels", "workbook.xml.rels")).getroot()
    for r in rels:
        if r.get("Id") == rid:
            return os.path.join(X, r.get("Target"))


def cells(name):
    res = {}
    for c in ET.parse(sheet_path(name)).getroot().iter("{%s}c" % NS["m"]):
        v = c.find("m:v", NS)
        if v is None:
            is_ = c.find("m:is", NS)
            if is_ is None:
                continue
            val = "".join(t.text or "" for t in is_.iter("{%s}t" % NS["m"]))
        elif c.get("t") == "s":
            val = sst[int(v.text)]
        elif c.get("t") in ("str", "inlineStr", "b", "e"):
            val = v.text
        else:
            val = float(v.text)
        res[c.get("r")] = val
    return res


P = cells("PRED-RMSDLIST")
T = [(str(P["AP%d" % r]).strip().upper(), P["AR%d" % r]) for r in range(8, 110)]
out["T_n"] = len(T)
out["T_unique_ids"] = len({i for i, _ in T})
out["T_median_all"] = statistics.median(v for _, v in T)
S = cells("Summary")
plddt = {}
for k, v in S.items():
    m = re.fullmatch(r"B(\d+)", k)
    if m and isinstance(v, str) and re.fullmatch(r"[0-9][A-Za-z0-9]{3}", v.strip()):
        h = S.get("H" + m.group(1))
        if isinstance(h, float):
            plddt.setdefault(v.strip().upper(), []).append(h)
out["summary_ids_with_plddt"] = len(plddt)
out["summary_dup_ids"] = sorted(k for k, v in plddt.items() if len(v) > 1)

ts = json.load(open(os.path.join(ROOT, "results", "e428", "t_structures.json")))
tj = {r["pdb"].upper(): r for r in ts}
mism = []
for i, v in T:
    r = tj.get(i)
    if r is None or abs(r["rmsd"] - v) > 1e-9:
        mism.append((i, v, None if r is None else r["rmsd"]))
    p = plddt.get(i)
    if r is not None and (p is None or abs(p[0] - (r.get("plddt") or -1)) > 1e-6):
        mism.append(("plddt", i, p, r.get("plddt")))
out["T_json_vs_xlsx_mismatches"] = mism
out["T_ids_equal_set"] = set(tj) == {i for i, _ in T}

# ---------------------------------------------------------------- stats
CLOSE = {"SAME_ACCESSION", "SEQ95"}


def med_diff(a, b):
    return statistics.median(a) - statistics.median(b)


def boot_ci(a, b, n, seed):
    rng = random.Random(seed)
    ds = sorted(med_diff([rng.choice(a) for _ in a], [rng.choice(b) for _ in b]) for _ in range(n))
    # percentile, linear interpolation (numpy default)
    def q(p):
        x = p * (n - 1); lo = math.floor(x); hi = min(lo + 1, n - 1)
        return ds[lo] + (ds[hi] - ds[lo]) * (x - lo)
    return q(0.025), q(0.975)


def perm_p(a, b, n, seed):
    rng = random.Random(seed)
    obs = abs(med_diff(a, b)); pool = a + b; k = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(med_diff(pool[:len(a)], pool[len(a):])) >= obs - 1e-12:
            k += 1
    return (k + 1) / (n + 1)


def compare(val):
    a = [val[i] for i, r in tj.items() if r["class"] in CLOSE and val.get(i) is not None]
    b = [val[i] for i, r in tj.items() if r["class"] not in CLOSE | {"QUERY_FAILED"} and val.get(i) is not None]
    lo, hi = boot_ci(a, b, 10000, 12345)
    return {"n_close": len(a), "n_other": len(b), "median_close": statistics.median(a),
            "median_other": statistics.median(b), "diff": med_diff(a, b), "ci95": [lo, hi],
            "perm_p": perm_p(a, b, 10000, 54321), "supported_ci_upper_lt_0": hi < 0}


rm = {i: v for i, v in T}
out["P1_rmsd"] = compare(rm)
out["S2_plddt"] = compare({i: plddt[i][0] for i in tj if i in plddt})
out["S4_rmsd_median_by_class"] = {}
for c in ["SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL", "QUERY_FAILED"]:
    v = [rm[i] for i, r in tj.items() if r["class"] == c]
    out["S4_rmsd_median_by_class"][c] = (len(v), statistics.median(v) if v else None)
out["S4_plddt_median_by_class"] = {c: statistics.median([plddt[i][0] for i, r in tj.items() if r["class"] == c])
                                   for c in ["SAME_ACCESSION", "SEQ95", "SEQ30", "NOVEL"]}

# ---------------------------------------------------------------- cache
GET, SEQ, ACC = {}, {}, {}
for f in glob.glob(os.path.join(D, "http", "*.json")):
    d = json.load(open(f))
    assert hashlib.sha256((d["body"] or "").encode()).hexdigest() == d["sha256"]
    if d["method"] == "GET":
        GET[d["url"]] = d
        continue
    q = json.loads(d["request_body"])["query"]
    nodes = q["nodes"] if q.get("type") == "group" else [q]
    dated = any(n["parameters"].get("attribute") == "rcsb_accession_info.initial_release_date"
                and n["parameters"]["value"] == CUT_AF2 and n["parameters"]["operator"] == "less_or_equal" for n in nodes)
    ids = json.loads(d["body"])["result_set"] if d["status"] == 200 else ([] if d["status"] == 204 else None)
    for n in nodes:
        p = n["parameters"]
        if n["service"] == "sequence":
            assert p["evalue_cutoff"] == 0.1
            SEQ[(p["value"], p["identity_cutoff"], dated)] = ids
        elif "database_accession" in p.get("attribute", ""):
            ACC[(p["value"], dated)] = ids


def classify(seq, accs):
    s95, s30 = SEQ.get((seq, 0.95, True), "MISSING"), SEQ.get((seq, 0.3, True), "MISSING")
    ah = [ACC.get((a, True), "MISSING") for a in accs]
    if "MISSING" in (s95, s30) or "MISSING" in ah or None in (s95, s30) or None in ah:
        return "QUERY_FAILED/MISSING", None
    if any(ah):
        c = "SAME_ACCESSION"
    elif s95:
        c = "SEQ95"
    elif s30:
        c = "SEQ30"
    else:
        c = "NOVEL"
    return c, {"n_acc_hits": sum(len(x) for x in ah), "n95": len(s95), "n30": len(s30)}


def t_class(pdb):
    e = GET.get("https://data.rcsb.org/rest/v1/core/entry/%s" % pdb)
    if e is None or e["status"] != 200:
        return "QUERY_FAILED", {"entry_status": None if e is None else e["status"]}
    ids = sorted(int(x) for x in json.loads(e["body"])["rcsb_entry_container_identifiers"]["polymer_entity_ids"])
    for k in ids:
        pe = GET.get("https://data.rcsb.org/rest/v1/core/polymer_entity/%s/%d" % (pdb, k))
        if pe is None:
            return "ENTITY_NOT_CACHED", {"entity": k}
        b = json.loads(pe["body"])
        if b["entity_poly"]["rcsb_entity_polymer_type"] != "Protein":
            continue
        seq = b["entity_poly"]["pdbx_seq_one_letter_code_can"].replace("\n", "")
        accs = sorted({r["database_accession"] for r in
                       (b["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers") or [])
                       if r["database_name"] == "UniProt"})
        c, info = classify(seq, accs)
        return c, dict(info or {}, entity=k, accs=accs)
    return "NO_PROTEIN_ENTITY", {}


def sifts(pdb, chain=None):
    g = GET.get("https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/%s" % pdb.lower())
    if g is None:
        return None
    if g["status"] != 200:
        return []
    u = json.loads(g["body"])[pdb.lower()]["UniProt"]
    return sorted(a for a, v in u.items() if chain is None or any(m["chain_id"] == chain for m in v["mappings"]))


test = list(csv.DictReader(open(os.path.join(D, "alphaflow", "atlas_test.csv"))))
train = list(csv.DictReader(open(os.path.join(D, "alphaflow", "atlas_train.csv"))))
train_ids = {r["name"].split("_")[0].upper() for r in train}
train_acc = {}
n_train_nosifts = 0
for r in train:
    p, ch = r["name"].rsplit("_", 1)
    a = sifts(p, ch)
    if a is None or a == []:
        n_train_nosifts += 1
    for x in a or []:
        train_acc.setdefault(x, set()).add(p.upper())
out["train_chains"] = len(train)
out["train_chains_without_sifts_accession"] = n_train_nosifts


def a_class(row):
    p, ch = row["name"].rsplit("_", 1)
    accs = sifts(p, ch)
    if accs is None:
        return "SIFTS_NOT_CACHED", {}
    c, info = classify(row["seqres"], accs)
    n95 = SEQ.get((row["seqres"], 0.95, False)) or []
    n30 = SEQ.get((row["seqres"], 0.3, False)) or []
    tr = {"tr_same_acc": sorted({q for a in accs for q in train_acc.get(a, ())}),
          "tr95": sorted(set(n95) & train_ids), "tr30": sorted(set(n30) & train_ids),
          "undated_cached": (row["seqres"], 0.95, False) in SEQ and (row["seqres"], 0.3, False) in SEQ}
    return c, dict(info or {}, accs=accs, **tr)


# ---------------------------------------------------------------- class recheck
aj = {r["name"]: r for r in json.load(open(os.path.join(ROOT, "results", "e428", "a_chains.json")))}
rng = random.Random(7)
t_ids = sorted(tj); a_ids = sorted(aj)
t_s = sorted(rng.sample(t_ids, math.ceil(0.2 * len(t_ids))))
a_s = sorted(rng.sample(a_ids, math.ceil(0.2 * len(a_ids))))
out["sample_T"], out["sample_A"] = t_s, a_s
res_t = {i: t_class(i) for i in t_ids}
rows = {r["name"]: r for r in test}
res_a = {i: a_class(rows[i]) for i in a_ids}
out["sample_T_mismatch"] = [(i, res_t[i][0], tj[i]["class"]) for i in t_s if res_t[i][0] != tj[i]["class"]]
out["sample_A_mismatch"] = [(i, res_a[i][0], aj[i]["class"]) for i in a_s if res_a[i][0] != aj[i]["class"]]
out["full_T_mismatch"] = [(i, res_t[i][0], tj[i]["class"], res_t[i][1]) for i in t_ids if res_t[i][0] != tj[i]["class"]]
out["full_A_mismatch"] = [(i, res_a[i][0], aj[i]["class"]) for i in a_ids if res_a[i][0] != aj[i]["class"]]
out["T_entity_not_1"] = [(i, res_t[i][1].get("entity")) for i in t_ids if res_t[i][1].get("entity") not in (1, None)]
out["T_no_accession"] = [i for i in t_ids if res_t[i][1].get("accs") == []]
out["A_no_accession"] = [i for i in a_ids if res_a[i][1].get("accs") == []]


def cnt(c, res):
    from collections import Counter
    return dict(Counter(v[0] for v in res.values()))


out["T_class_counts_recomputed"] = cnt(None, res_t)
out["A_class_counts_recomputed"] = cnt(None, res_a)
out["A_released_le_2020_05_01"] = sum(r["release_date"] <= "2020-05-01" for r in test)
out["A_test_n"] = len(test)
out["A_train_overlap"] = {
    "same_accession": sum(bool(res_a[i][1]["tr_same_acc"]) for i in a_ids),
    "seq95": sum(bool(res_a[i][1]["tr95"]) for i in a_ids),
    "seq30": sum(bool(res_a[i][1]["tr30"]) for i in a_ids),
    "any": sum(bool(res_a[i][1]["tr_same_acc"] or res_a[i][1]["tr30"]) for i in a_ids),
    "all_undated_cached": all(res_a[i][1]["undated_cached"] for i in a_ids)}
out["A_train_overlap_vs_json_mismatch"] = [
    i for i in a_ids if (sorted(res_a[i][1]["tr95"]) != sorted(aj[i]["train_split_seq95"])
                         or sorted(res_a[i][1]["tr30"]) != sorted(aj[i]["train_split_seq30"])
                         or bool(res_a[i][1]["tr_same_acc"]) != bool(aj[i]["train_split_same_accession"]))]
out["A_fusion"] = {n: res_a[n][1]["accs"] for n in ["6xds_A", "6xrx_A", "7p41_D"]}
out["A_fusion_per_accession_training_era_hits"] = {a: len(ACC.get((a, True)) or []) for n in out["A_fusion"]
                                                   for a in out["A_fusion"][n]}
out["A_multi_accession_chains"] = {i: res_a[i][1]["accs"] for i in a_ids if len(res_a[i][1].get("accs", [])) > 1}
e = GET["https://data.rcsb.org/rest/v1/core/entry/7DRH"]
out["7DRH"] = {"status": e["status"], "body": e["body"][:120]}
json.dump(out, sys.stdout, indent=1, default=list)
